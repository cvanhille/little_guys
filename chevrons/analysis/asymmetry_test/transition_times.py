import numpy as np
import pandas as pd
from ovito.io import *
import sys
import glob
from tqdm.auto import tqdm

def angles(pos,typ,dtyp):
    v1 = pos[typ == 3+dtyp][0] - pos[typ == 4+dtyp][0]
    v2 = pos[typ == 5+dtyp][0] - pos[typ == 4+dtyp][0]
    v1 /= np.linalg.norm(v1)
    v2 /= np.linalg.norm(v2)
    dp = np.dot(v1,v2)
    ang = np.arccos(dp)
    atip = ang/np.pi*180
    v1 = pos[typ == 2+dtyp][0] - pos[typ == 1+dtyp][0]
    v2 = pos[typ == 6+dtyp][0] - pos[typ == 1+dtyp][0]
    v1 /= np.linalg.norm(v1)
    v2 /= np.linalg.norm(v2)
    dp = np.dot(v1,v2)
    ang = np.arccos(dp)
    abac = ang/np.pi*180
    return atip, abac

def ats(p,dtyp):
    tips = []
    bacs = []
    frms = []
    dists = []
    for fr in range(p.source.num_frames+1):
        d = p.compute(fr)
        pos = np.array(d.particles['Position'])[:,:2]
        typ = np.array(d.particles['Particle Type'])
        atip, abac = angles(pos,typ,dtyp)
        d1 = np.linalg.norm(pos[typ == 1][0] - pos[typ == 10][0])
        d2 = np.linalg.norm(pos[typ == 2][0] - pos[typ == 9][0])
        d3 = np.linalg.norm(pos[typ == 6][0] - pos[typ == 11][0])
        dists.append([d1,d2,d3])
        tips.append(atip)
        bacs.append(abac)
        frms.append(fr)
    frms = np.array(frms)
    tips = np.array(tips)
    bacs = np.array(bacs)
    dists = np.array(dists)
    return frms, tips, bacs, dists

if len(sys.argv) < 3:
    print()
    print("Error! Wrong number of arguments. Please provide: 1) general path to sample 2) key for type of binding (LB/RB) 3) binding mde (dist/nodist)")
    print()
    exit()

gpath = sys.argv[1]
key = sys.argv[2]
mode = sys.argv[3]

if key == 'LB' or key == 'L' or key == 'left':
    dtyp = 0
elif key == 'RB' or key == 'R' or key == 'right':
    dtyp = 6

files = glob.glob('%s/sd*/output.xyz'%(gpath))
print()
ck = input("%d simulations to analyse! Continue? [y/n] "%(len(files)))
print()
if not (ck == 'y' or ck == 'yes' or ck == 'Y'):
    print("Ok! Aborting...")
    print()
    exit()

seeds = []
times = []
for i, file in enumerate(tqdm(files)):
    seed = int(file.split('/output.xyz')[0].split('/sd')[1])
    p = import_file(file)
    frms, tips, bacs, dists = ats(p,dtyp)
    shape = np.abs(bacs-120)/60
    shape[shape > 0.8] = 1.0
    shape[shape < 0.2] = 0.0
    shape[shape%1.0 > 0] = np.nan
    sw = 120+shape*60

    if mode == 'dist':
        tchevs = frms[shape == 0]
        tclose = frms[np.mean(dists,1) < 1.02]
        tvals = np.intersect1d(tchevs,tclose)
        if len(tvals) == 0:
            time = np.nan
        else:
            time = np.min(tvals)
    else:
        time = np.min(frms[shape == 0])

    if np.isnan(time):
        continue

    seeds.append(seed)
    times.append(time)

data = pd.DataFrame(index = np.arange(len(seeds)), columns = ['seed', 'time'])
data['seed'] = seeds
data['time'] = times
data.to_csv('%s/binding_times.txt'%(gpath))
