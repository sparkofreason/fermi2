import fermi.tipsh as tipsh
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename
from mmappickle.dict import mmapdict
import datetime
import numpy
import time

data = '/home/dave/data/fermi-data'
alphas = [0.1, 0.05, 0.03, 0.01, 0.005, 0.003, 0.001, 5e-4, 3e-4, 1e-4, 5e-5, 3e-5, 1e-5, 5e-6, 3e-6, 1e-6, 5e-7, 3e-7, 1e-7, 5e-8, 3e-8, 1e-8]
fits_root = 'fits_files'
jmin = 0
fwer = 'uniform'
energy = '500-1000gev'

counts = mmapdict('counts.mmdpickle')
total_models = mmapdict('total_models_v02.mmdpickle')

date = datetime.date.today()

tic = time.time()
for alpha in alphas:
    print('----->', alpha)
    filename = f"{data}/alphas_v02_modified_{date}_{energy}_{fwer}_{alpha}.mmdpickle"
    tipsh.run_tipsh(counts[energy], total_models[energy], alpha, jmin, fwer, filename, poolWorkers=8)
toc = time.time()
print("Elapsed", toc - tic)