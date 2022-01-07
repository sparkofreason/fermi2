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