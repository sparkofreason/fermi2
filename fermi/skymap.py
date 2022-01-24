import matplotlib.pyplot as plt
import matplotlib.colors as colors
import cartopy.crs as ccrs

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

def imshow_carree(data, cmap, norm, interpolation='antialiased', extent=(-180,180,-90,90)):
    fig = plt.figure(figsize=(30, 15))
    ax = plt.subplot()
    #diff = plt.imshow(fits.getdata(data_file, ext=0), origin='lower')
    diff = ax.imshow(data, cmap=cmap, origin='lower', norm=norm, extent=extent, interpolation=interpolation)
    ax.grid(color='black')
    #diff = plt.imshow(im0, cmap=plt.cm.viridis, origin='lower', norm=colors.LogNorm())
    plt.colorbar(diff)
    return fig, ax

def imshow_multiple(data, keys, cmap, norms):
    N = len(keys)
    fig = plt.figure(figsize=(30, 15*N))
    for n in range(0, N):
        ax = plt.subplot(N, 1, n+1)
        ax.grid()
        diff = ax.imshow(data[keys[n]], cmap=cmap, origin='lower', norm=norms[keys[n]], extent=(-180,180,-90,90))
        #ax.grid(color='black')
        plt.title(keys[n])
        plt.colorbar(diff)

