import fermi.tipsh as tipsh
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename
from mmappickle.dict import mmapdict
import datetime
import numpy
import time
from scipy.stats import poisson

energy = '5-10gev'
fits_root = 'fits_files'
data_file = get_pkg_data_filename(f"{fits_root}/lat_source_zmax90_{energy}_ccube.fits")
point_model_file = get_pkg_data_filename(f"{fits_root}/{energy}_source_point_model_map.fits")
diffuse_model_file = get_pkg_data_filename(f"{fits_root}/{energy}_source_diffuse_model_map.fits")
galactic_model_file = get_pkg_data_filename(f"{fits_root}/{energy}_source_galactic_model_map.fits")
galactic_v02_model_file = get_pkg_data_filename(f"{fits_root}/{energy}_source_galactic_v02_model_map.fits")
isotropic_model_file = get_pkg_data_filename(f"{fits_root}/{energy}_source_isotropic_model_map.fits")

count_cube = fits.getdata(data_file, ext=0)
count_data = count_cube[0]
point_model = fits.getdata(point_model_file, ext=0)
diffuse_model = fits.getdata(diffuse_model_file, ext=0)
galactic_model = fits.getdata(galactic_model_file, ext=0)
galactic_v02_model = fits.getdata(galactic_v02_model_file, ext=0)
isotropic_model = fits.getdata(isotropic_model_file, ext=0)

total_model = 1.0*point_model + 1.0*diffuse_model + 1.0*galactic_model + 1.0*isotropic_model
del count_cube, point_model, diffuse_model, galactic_model, isotropic_model

alpha = 1e-2
jmin = 0
fwer = 'uniform'
data = 'constant'
now = datetime.

if (data == 'simulated'):
    count_data = poisson.rvs(total_model, random_state=420)
elif (data == 'constant'):
    total_model = numpy.ones(total_model.shape)
    count_data = poisson.rvs(total_model, random_state=42)

tic = time.time()
tipsh.run_tipsh()
toc = time.time()
print("Elapsed", toc - tic)
#numpy.savez(filename, a=a, hs=hs, vs=vs, ds=ds, total_model=total_model)