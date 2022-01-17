import fermi.tipsh as tipsh
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename
import datetime
import numpy
import time

energy = '100-200mev'
fits_root = 'fits_files'
data_file = get_pkg_data_filename(f"{fits_root}/lat_source_zmax90_{energy}_ccube.fits")
point_model_file = get_pkg_data_filename(f"{fits_root}/{energy}_source_point_model_map.fits")
diffuse_model_file = get_pkg_data_filename(f"{fits_root}/{energy}_source_diffuse_model_map.fits")
galactic_model_file = get_pkg_data_filename(f"{fits_root}/{energy}_source_galactic_v02_model_map.fits")
isotropic_model_file = get_pkg_data_filename(f"{fits_root}/{energy}_source_isotropic_model_map.fits")

count_cube = fits.getdata(data_file, ext=0)
count_data = count_cube[0]
point_model = fits.getdata(point_model_file, ext=0)
diffuse_model = fits.getdata(diffuse_model_file, ext=0)
galactic_model = fits.getdata(galactic_model_file, ext=0)
isotropic_model = fits.getdata(isotropic_model_file, ext=0)

total_model = 1.0*point_model + 1.0*diffuse_model + 0.9*galactic_model + 1.2*isotropic_model
del count_cube, point_model, diffuse_model, galactic_model, isotropic_model

alpha = 0.01
jmin = 0
fwer = 'sidak'
mode = 'wavelet'
args = tipsh.skellam_inputs(tipsh.spherify(count_data), tipsh.spherify(total_model), alpha, jmin, fwer, mode)

tic = time.time()
a, hs, vs, ds = tipsh.haar_threshold_foo(args, 8)
toc = time.time()
print("Elapsed", toc - tic)

now = datetime.date.today()
filename = f"{now}_{energy}_{fwer}_{alpha}_{mode}"
print("Saving", filename)
numpy.savez(filename, a=a, hs=hs, vs=vs, ds=ds, total_model=total_model)