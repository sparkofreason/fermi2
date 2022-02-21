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
alpha = 1e-6
jmin = 0
fwer = 'uniform'

counts = mmapdict('counts.mmdpickle')
total_models = mmapdict('total_models_v02.mmdpickle')

date = datetime.date.today()

tic = time.time()
for energy in energies:
    print('----->', energy)
    filename = f"v02_modified_{date}_{energy}_{fwer}_{alpha}.mmdpickle"
    tipsh.run_tipsh(counts[energy], total_models[energy], alpha, jmin, fwer, filename, poolWorkers=8)
toc = time.time()
print("Elapsed", toc - tic)
