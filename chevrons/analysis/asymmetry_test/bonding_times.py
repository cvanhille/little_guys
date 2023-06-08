import numpy as np
import pandas as pd
from ovito.io import *
from ovito.modifiers import *
import sys
import glob
from tqdm.auto import tqdm

def times(path):
    file = '%s/output.xyz'%(path)
    bonds = '%s/bonds.dump'%(path)

    p = import_file(file)
    mod = LoadTrajectoryModifier()
    p.modifiers.append(mod)
    mod.source.load(bonds, columns = ['Particle Identifiers.1', 'Particle Identifiers.2', 'Bond Type', 'Energy.1', 'Force.1', 'Length'], multiple_frames = True)

    ff = p.source.num_frames
    nbonds = []
    frames = []
    # for i, f in enumerate(tqdm(range(ff+1))):
    for f in range(ff+1):
        d = p.compute(f)
        top = np.array(d.particles.bonds['Topology'])
        nbonds.append(len(top))
        frames.append(f)
        
    nbonds = np.array(nbonds)
    frames = np.array(frames)

    nbonds -= nbonds[0]

    if len(frames[nbonds >= 1]) == 0:
        tb1 = np.nan
        tb2 = np.nan
        tb3 = np.nan
    else:
        tb1 = np.min(frames[nbonds >= 1])
        if len(frames[nbonds >= 2]) == 0:
            tb2 = np.nan
            tb3 = np.nan
        else:
            tb2 = np.min(frames[nbonds >= 2])
            if len(frames[nbonds >= 3]) == 0:
                tb3 = np.nan
            else:
                tb3 = np.min(frames[nbonds >= 3])

    return tb1,tb2,tb3

if len(sys.argv) < 2:
    print()
    print("Error! Wrong number of arguments. Please provide: 1) general path to sample")
    print()
    exit()

gpath = sys.argv[1]

files = glob.glob('%s/sd*/output.xyz'%(gpath))
print()
ck = input("%d simulations to analyse! Continue? [y/n] "%(len(files)))
print()
if not (ck == 'y' or ck == 'yes' or ck == 'Y'):
    print("Ok! Aborting...")
    print()
    exit()

seeds = []
time1 = []
time2 = []
time3 = []
for i, file in enumerate(tqdm(files)):
    path = file.split('/output.xyz')[0]
    seed = int(path.split('/sd')[1])
    t1,t2,t3 = times(path)
    seeds.append(seed)
    time1.append(t1)
    time2.append(t2)
    time3.append(t3)

data = pd.DataFrame(index = np.arange(len(seeds)), columns = ['seed', 'time1', 'time2', 'time3'])
data['seed'] = seeds
data['time1'] = time1
data['time2'] = time2
data['time3'] = time3
data.to_csv('%s/bonding_times.txt'%(gpath))
