import fermi.tipsh as tipsh
from mmappickle.dict import mmapdict
import datetime
import time

data = '/home/dave/data/fermi-data'
energies = ['100-200mev', '200-500mev', '500-1000mev', '1-2gev', '2-5gev', '5-10gev', '10-20gev', '20-50gev', '50-100gev', '100-200gev', '200-500gev', '500-1000gev']
jmin = 0

counts = mmapdict('counts.mmdpickle', readonly=True)
total_models = mmapdict('total_models_v02.mmdpickle', readonly=True)

date = datetime.date.today()

tic = time.time()
for energy in energies:
    print('----->', energy)
    filename = f"{data}/v02_modified_{date}_{energy}_pvalues.mmdpickle"
    tipsh.run_pvalues(counts[energy], total_models[energy], jmin, filename, poolWorkers=8)
toc = time.time()
print("Elapsed", toc - tic)
