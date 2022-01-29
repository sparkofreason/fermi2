import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib import animation, rc
import cartopy.crs as ccrs
import fermi.tipsh as tipsh
import math
import numpy 

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
    fig = plt.figure(figsize=(30, 15))
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

def imshow_carree_impl(fig, ax, data, cmap, norm=None, interpolation='antialiased', extent=(-180,180,-90,90), linthresh=1.0):
    norm = norm if norm else autonorm(data, linthresh)
    im = ax.imshow(data, cmap=cmap, origin='lower', norm=norm, extent=extent, interpolation=interpolation)
    ax.grid(color='black')
    fig.colorbar(im, ax=ax)
    return im

def imshow_carree(data, cmap, norm=None, interpolation='antialiased', extent=(-180,180,-90,90), linthresh=1.0):
    fig = plt.figure(figsize=(30, 15))
    ax = plt.subplot()
    imshow_carree_impl(fig, ax, data, cmap, interpolation=interpolation, linthresh=linthresh)
    

def imshow_spherified_carree(data, cmap, norm, interpolation='antialiased'):
    rows, cols = data.shape
    fig = plt.figure(figsize=(30, 15))
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
    fig = plt.figure(figsize=(30, 15))
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
        ax.set_title(keys[n], fontsize=40)

def animate_carree(data, keys, cmap, linthresh=None, norms=None):
    N = len(keys)
    def animate(i):
        j, k = divmod(i, N)
        n = k if j % 2 == 0 else N - k - 1  
        #print(N, n, j, k)
        im.set_array(count_recs[energies[n]])
        im.set_norm(norms[energies[n]])
        ax.set_title(energies[n], fontsize=40)
        return [im]