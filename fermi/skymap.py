import base64
from io import BytesIO
from pathlib import Path
from tempfile import TemporaryDirectory
import uuid
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib import animation, rc, _api
from matplotlib._animation_data import (
    DISPLAY_TEMPLATE, INCLUDED_FRAMES, JS_INCLUDE, STYLE_INCLUDE)
import cartopy.crs as ccrs
import fermi.tipsh as tipsh
import math
import numpy 

figsize = (16, 8)
title_fontsize = 24

def symLogNormShift(data, linthresh, linscale=1.0, vmin=None, vmax=None):
    min_data = numpy.min(data)
    max_data = numpy.max(data)
    mean_data = numpy.median(data)
    symLogNorm = colors.SymLogNorm(linthresh=linthresh,linscale=linscale,vmin=min_data-mean_data,vmax=max_data-mean_data)
    return colors.FuncNorm((lambda x: symLogNorm(x - mean_data), lambda x: symLogNorm.inverse(x) + mean_data), vmin=min_data, vmax=max_data)

def ticks(data, linthresh):
    min_data = numpy.min(data)
    max_data = numpy.max(data)
    n1 = math.floor(math.log10(linthresh))
    n2 = math.floor(math.log10(-min_data))
    n3 = math.floor(math.log10(max_data))
    return numpy.concatenate((-numpy.flip(numpy.logspace(n1, n2, n2-n1+1)), [0], numpy.logspace(n1, n3, n3-n1+1)))

def imshow_mollweide(data, cmap, norm):
    fig = plt.figure(figsize=figsize)
    ax = plt.subplot(projection=ccrs.Mollweide())
    #diff = plt.imshow(fits.getdata(data_file, ext=0), origin='lower')
    diff = ax.imshow(data, cmap=cmap, origin='lower', norm=norm, transform=ccrs.PlateCarree(), extent=(-180,180,-90,90))
    #ax.grid(color='black')
    #diff = plt.imshow(im0, cmap=plt.cm.viridis, origin='lower', norm=colors.LogNorm())
    plt.colorbar(diff)
    return fig, ax

def imshow_mollweide_multiple(data, keys, cmap, norm):
    N = len(keys)
    fig = plt.figure(figsize=(30, 15*N))
    for n in range(0, N):
        ax = plt.subplot(N, 1, n+1, projection=ccrs.Mollweide())
        diff = ax.imshow(data[keys[n]], cmap=cmap, origin='lower', norm=norm, transform=ccrs.PlateCarree(), extent=(-180,180,-90,90))
        #ax.grid(color='black')
        plt.colorbar(diff)
    #return fig

def autonorm(data, linthresh):
    min_data = numpy.min(data)
    max_data = numpy.max(data)
    if (min_data < 0 and max_data > 0):
        if (-min_data > max_data):
            vmin = min_data
            vmax = -min_data
        else:
            vmin = -max_data
            vmax = max_data
    elif max_data > 0:
        vmin = 0
        vmax = max_data
    else:
        vmin = min_data
        vmax = 0
    return colors.SymLogNorm(linthresh=linthresh, vmin=vmin, vmax=vmax)

def imshow_carree_impl(fig, ax, data, cmap, norm=None, interpolation='antialiased', extent=(-180,180,-90,90), linthresh=1.0, title=None, shrink=0.775, ticks=None, tick_labels=None):
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    fig.set_tight_layout({'pad': 0})
    norm = norm if norm else autonorm(data, linthresh)
    im = ax.imshow(data, cmap=cmap, origin='lower', norm=norm, extent=extent, interpolation=interpolation)
    ax.grid(color='black')
    if title:
        ax.set_title(title, fontsize=title_fontsize)
    cbar = fig.colorbar(im, ax=ax, shrink=shrink, ticks=ticks)
    cbar.ax.tick_params(labelsize=16)
    if tick_labels:
        cbar.ax.set_yticklabels(tick_labels)
    return im

def imshow_carree(data, cmap, norm=None, interpolation='antialiased', extent=(-180,180,-90,90), linthresh=1.0, title=None, ticks=None, tick_labels=None):
    fig = plt.figure(figsize=figsize)
    ax = plt.subplot()
    imshow_carree_impl(fig, ax, data, cmap, interpolation=interpolation, linthresh=linthresh, norm=norm, title=title, ticks=ticks, tick_labels=tick_labels)
    

def imshow_spherified_carree(data, cmap, norm, interpolation='antialiased'):
    rows, cols = data.shape
    fig = plt.figure(figsize=figsize)
    ax1 = plt.subplot(2, 1, 2)
    #diff = plt.imshow(fits.getdata(data_file, ext=0), origin='lower')
    im1 = ax1.imshow(tipsh.unspherify(data), cmap=cmap, origin='lower', norm=norm, extent=(-180,180,-90,90), interpolation=interpolation)
    ax1.grid(color='black')
    #diff = plt.imshow(im0, cmap=plt.cm.viridis, origin='lower', norm=colors.LogNorm())
    plt.colorbar(im1)

    ax2 = plt.subplot(2, 1, 1)
    #diff = plt.imshow(fits.getdata(data_file, ext=0), origin='lower')
    im2 = ax2.imshow(tipsh.unspherify_top(data), cmap=cmap, origin='lower', norm=norm, extent=(0, 360,90,-90), interpolation=interpolation)
    ax2.grid(color='black')
    #diff = plt.imshow(im0, cmap=plt.cm.viridis, origin='lower', norm=colors.LogNorm())
    plt.colorbar(im2)


def imshow(data, cmap, norm, interpolation='antialiased'):
    fig = plt.figure(figsize=figsize)
    ax = plt.subplot()
    diff = ax.imshow(data, cmap=cmap, origin='lower', norm=norm, interpolation=interpolation)
    ax.grid(color='black')
    plt.colorbar(diff)


def imshow_multiple(data, keys, cmap, linthresh=None, norms=None, interpolation='antialiased'):
    N = len(keys)
    fig = plt.figure(figsize=(30, 15*N))
    for n in range(0, N):
        ax = plt.subplot(N, 1, n+1)
        imshow_carree_impl(fig, ax, data[keys[n]], cmap, interpolation=interpolation, linthresh=linthresh[keys[n]])
        ax.set_title(keys[n], fontsize=title_fontsize)

# Taken directly from jakevdp's JSAnimation package at
# http://github.com/jakevdp/JSAnimation
def _included_frames(paths, frame_format):
    """paths should be a list of Paths"""
    return INCLUDED_FRAMES.format(Nframes=len(paths),
                                  frame_dir=paths[0].parent,
                                  frame_format=frame_format)


def _embedded_frames(frame_list, frame_format):
    """frame_list should be a list of base64-encoded png files"""
    if frame_format == 'svg':
        # Fix MIME type for svg
        frame_format = 'svg+xml'
    template = '  frames[{0}] = "data:image/{1};base64,{2}"\n'
    return "\n" + "".join(
        template.format(i, frame_format, frame_data.replace('\n', '\\\n'))
        for i, frame_data in enumerate(frame_list))


#@writers.register('html')
class HTMLWriter(animation.FileMovieWriter):
    """Writer for JavaScript-based HTML movies."""

    supported_formats = ['png', 'jpeg', 'tiff', 'svg']

    @classmethod
    def isAvailable(cls):
        return True

    def __init__(self, fps=30, codec=None, bitrate=None, extra_args=None,
                 metadata=None, embed_frames=False, default_mode='loop',
                 embed_limit=None):

        #if extra_args:
        #    _log.warning("HTMLWriter ignores 'extra_args'")
        extra_args = ()  # Don't lookup nonexistent rcParam[args_key].
        self.embed_frames = embed_frames
        self.default_mode = default_mode.lower()
        _api.check_in_list(['loop', 'once', 'reflect'],
                           default_mode=self.default_mode)

        # Save embed limit, which is given in MB
        if embed_limit is None:
            self._bytes_limit = plt.rcParams['animation.embed_limit']
        else:
            self._bytes_limit = embed_limit
        # Convert from MB to bytes
        self._bytes_limit *= 1024 * 1024

        super().__init__(fps, codec, bitrate, extra_args, metadata)

    def setup(self, fig, outfile, dpi, frame_dir=None):
        outfile = Path(outfile)
        _api.check_in_list(['.html', '.htm'], outfile_extension=outfile.suffix)

        self._saved_frames = []
        self._total_bytes = 0
        self._hit_limit = False

        if not self.embed_frames:
            if frame_dir is None:
                frame_dir = outfile.with_name(outfile.stem + '_frames')
            frame_dir.mkdir(parents=True, exist_ok=True)
            frame_prefix = frame_dir / 'frame'
        else:
            frame_prefix = None

        super().setup(fig, outfile, dpi, frame_prefix)
        self._clear_temp = False

    def grab_frame(self, **savefig_kwargs):
        if self.embed_frames:
            # Just stop processing if we hit the limit
            if self._hit_limit:
                return
            f = BytesIO()
            self.fig.savefig(f, format=self.frame_format,
                             dpi=self.dpi, **savefig_kwargs)
            imgdata64 = base64.encodebytes(f.getvalue()).decode('ascii')
            self._total_bytes += len(imgdata64)
            if self._total_bytes >= self._bytes_limit:
                #_log.warning(
                #    "Animation size has reached %s bytes, exceeding the limit "
                #    "of %s. If you're sure you want a larger animation "
                #    "embedded, set the animation.embed_limit rc parameter to "
                #    "a larger value (in MB). This and further frames will be "
                #    "dropped.", self._total_bytes, self._bytes_limit)
                self._hit_limit = True
            else:
                self._saved_frames.append(imgdata64)
        else:
            return super().grab_frame(**savefig_kwargs)

    def finish(self):
        # save the frames to an html file
        if self.embed_frames:
            fill_frames = _embedded_frames(self._saved_frames,
                                           self.frame_format)
            Nframes = len(self._saved_frames)
        else:
            # temp names is filled by FileMovieWriter
            fill_frames = _included_frames(self._temp_paths, self.frame_format)
            Nframes = len(self._temp_paths)
        mode_dict = dict(once_checked='',
                         loop_checked='',
                         reflect_checked='')
        mode_dict[self.default_mode + '_checked'] = 'checked'

        interval = 1000 // self.fps

        with open(self.outfile, 'w') as of:
            of.write(JS_INCLUDE + STYLE_INCLUDE)
            of.write(DISPLAY_TEMPLATE.format(id=uuid.uuid4().hex,
                                             Nframes=Nframes,
                                             fill_frames=fill_frames,
                                             interval=interval,
                                             **mode_dict))

        # duplicate the temporary file clean up logic from
        # FileMovieWriter.cleanup.  We can not call the inherited
        # versions of finish or cleanup because both assume that
        # there is a subprocess that we either need to call to merge
        # many frames together or that there is a subprocess call that
        # we need to clean up.
        if self._tmpdir:
            #_log.debug('MovieWriter: clearing temporary path=%s', self._tmpdir)
            self._tmpdir.cleanup()

def to_jshtml(self, fps=None, embed_frames=True, default_mode=None):
    """
    Generate HTML representation of the animation.
    Parameters
    ----------
    fps : int, optional
        Movie frame rate (per second). If not set, the frame rate from
        the animation's frame interval.
    embed_frames : bool, optional
    default_mode : str, optional
        What to do when the animation ends. Must be one of ``{'loop',
        'once', 'reflect'}``. Defaults to ``'loop'`` if ``self.repeat``
        is True, otherwise ``'once'``.
    """
    if fps is None and hasattr(self, '_interval'):
        # Convert interval in ms to frames per second
        fps = 1000 / self._interval

    # If we're not given a default mode, choose one base on the value of
    # the repeat attribute
    if default_mode is None:
        default_mode = 'loop' if self.repeat else 'once'

    if not hasattr(self, "_html_representation"):
        # Can't open a NamedTemporaryFile twice on Windows, so use a
        # TemporaryDirectory instead.
        with TemporaryDirectory() as tmpdir:
            path = Path(tmpdir, "temp.html")
            writer = HTMLWriter(fps=fps,
                                embed_frames=embed_frames,
                                default_mode=default_mode)
            self.save(str(path), writer=writer)
            self._html_representation = path.read_text()

    return self._html_representation

def animate_carree(data_fn, title_fn, N, cmap, linthresh=None, norms=None):
    fig = plt.figure(figsize=figsize)
    ax = plt.axes(xlim=(-180,180), ylim=(-90,90))
    #ax.grid()
    im = imshow_carree_impl(fig, ax, data_fn(0), cmap, norms(0) if norms else autonorm(data_fn(0), linthresh(0)))
    tx = ax.set_title(title_fn(0), fontsize=title_fontsize)
    def animate(i):
        im.set_data(data_fn(i))
        im.set_norm(norms(i) if norms else autonorm(data_fn(i), linthresh(i)))
        tx.set_text(title_fn(i))
        return [im]
    anim = animation.FuncAnimation(
                                fig, 
                                animate, 
                                frames = N,
                                interval = 1000 , # in ms
                                )
    plt.close()
    #plt.rcParams['animation.ffmpeg_path'] = '/usr/bin/ffmpeg'
    #FFwriter = animation.FFMpegWriter(fps=1, extra_args=['-vcodec', 'libx264'])
    #anim.save('unmodified_models_animation.mp4', writer=FFwriter)
    #return to_jshtml(anim)
    return anim