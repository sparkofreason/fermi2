import fermi.tipsh as tipsh
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename
from mmappickle.dict import mmapdict
import datetime
import numpy
import time

energies = ['100-200mev', '200-500mev', '500-1000mev', '1-2gev', '2-5gev', '5-10gev', '10-20gev', '20-50gev', '50-100gev', '100-200gev', '200-500gev', '500-1000gev']
fits_root = 'fits_files'
mode = 'wavelet'
now = datetime.date.today()
date = '2022-01-23'
alpha = 0.99
jmin = 0
fwer = 'sidak'

counts = mmapdict('counts.mmdpickle')
total_models = mmapdict('total_models_v02.mmdpickle')

now = datetime.date.today()
filename = f"v02_modified_{now}_{fwer}_{alpha}_{mode}.mmdpickle"
m = mmapdict(filename)
m.vacuum()

tic = time.time()
for energy in energies:
    print('----->', energy)
    args = tipsh.skellam_inputs(tipsh.spherify(counts[energy]), tipsh.spherify(total_models[energy]), alpha, jmin, fwer, mode)
    a, hs, vs, ds = tipsh.haar_threshold_pool(args, 8)
    count_rec = tipsh.inv_haar_sphere(a, hs, vs, ds)
    results = {'a': a, 'hs': hs, 'vs': vs, 'ds': ds, 'count_rec': count_rec}
    m[energy] = results
toc = time.time()
print("Elapsed", toc - tic)