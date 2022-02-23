import fermi.tipsh as tipsh
from mmappickle.dict import mmapdict
import datetime
import time

data = '/home/dave/data/fermi-data'
energies = ['100-200mev', '200-500mev', '500-1000mev', '1-2gev', '2-5gev', '5-10gev', '10-20gev', '20-50gev', '50-100gev', '100-200gev', '200-500gev', '500-1000gev']
jmin = 0
fwer = 'uniform'
alpha = 1e-6
level_recs = False

counts = mmapdict('counts.mmdpickle', readonly=True)
total_models = mmapdict('total_models_v02.mmdpickle', readonly=True)

date = datetime.date.today()

if fwer == 'sidak':
    alpha_fn = lambda j: tipsh.sidak(alpha, j)
else:
    alpha_fn = lambda j: alpha

tic = time.time()
for energy in energies:
    print('----->', energy)
    pvalues_filename = f"{data}/v02_modified_{date}_{energy}_pvalues.mmdpickle"
    filename = f"{data}/v02_modified_{date}_{energy}_{fwer}_{alpha}_{level_recs}.mmdpickle"
    tipsh.run_threshold(pvalues_filename, alpha_fn, jmin, filename, level_recs=level_recs, poolWorkers=8)
toc = time.time()
print("Elapsed", toc - tic)
