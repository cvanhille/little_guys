import numpy as np
import pandas as pd
from ovito.io import *
import sys
import glob
from tqdm.auto import tqdm

def angles(pos,typ):
    v1 = pos[typ == 3][0] - pos[typ == 4][0]
    v2 = pos[typ == 5][0] - pos[typ == 4][0]
    # PBCs
    v1[v1 >= 2.5] -= 5.0
    v1[v1 < -2.5] += 5.0
    v2[v2 >= 2.5] -= 5.0
    v2[v2 < -2.5] += 5.0
    v1 /= np.linalg.norm(v1)
    v2 /= np.linalg.norm(v2)
    dp = np.dot(v1,v2)
    ang = np.arccos(dp)
    atip = ang/np.pi*180
    v1 = pos[typ == 2][0] - pos[typ == 1][0]
    v2 = pos[typ == 6][0] - pos[typ == 1][0]
    # PBCs
    v1[v1 >= 2.5] -= 5.0
    v1[v1 < -2.5] += 5.0
    v2[v2 >= 2.5] -= 5.0
    v2[v2 < -2.5] += 5.0
    v1 /= np.linalg.norm(v1)
    v2 /= np.linalg.norm(v2)
    dp = np.dot(v1,v2)
    ang = np.arccos(dp)
    abac = ang/np.pi*180
    return atip, abac

def ats(p):
    tips = []
    bacs = []
    frms = []
    for fr in range(p.source.num_frames+1):
        d = p.compute(fr)
        pos = np.array(d.particles['Position'])[:,:2]
        typ = np.array(d.particles['Particle Type'])
        atip, abac = angles(pos,typ)
        tips.append(atip)
        bacs.append(abac)
        frms.append(fr)
    frms = np.array(frms)
    tips = np.array(tips)
    bacs = np.array(bacs)

    diff = np.abs(bacs-tips)
    
    shape = np.abs(bacs-120)/60
    shape[shape > 0.8] = 1.0
    shape[shape < 0.2] = 0.0
    shape[shape%1.0 > 0] = np.nan
    sw = 120+shape*60
    
    cfrs = frms[shape >= 0.0]
    cshp = shape[shape >= 0.0]
    
    dfrs = cfrs[1:]
    dshp = cshp[1:]-cshp[:-1]
    
    cdfrs = dfrs[np.abs(dshp)>0]
    cdshp = dshp[np.abs(dshp)>0]

    times = cdfrs[1:]-cdfrs[:-1]
    times = np.concatenate(([cdfrs[0]],times))
    shapes = np.zeros(len(times))
    shapes[1:] += cdshp[:-1]
    if cdshp[0] == 1:
        shapes[0] = -1
    else:
        shapes[0] = 1
    
    # Cat == -1
    # Triangle == 1
    ctimes = times[shapes == -1]
    ttimes = times[shapes == 1]
    
    return ctimes, ttimes

if len(sys.argv) < 2:
    print()
    print("Error! Wrong number of arguments. Please provide: 1) general path to sample")
    print()
    exit()

gpath = sys.argv[1]

files = glob.glob('%s/sd*/output.xyz'%(gpath)) 

seeds = []
ctimes = []
ttimes = []
for i, file in enumerate(tqdm(files)):

    p = import_file(file)

    ctimesS, ttimesS = ats(p)
    seed = int(file.split('/')[-2].split('sd')[1])
    sds = seed*np.ones(len(ctimesS))

    ctimes = np.concatenate((ctimes, ctimesS))
    ttimes = np.concatenate((ttimes, ttimesS))
    seeds = np.concatenate((seeds, sds))

data = pd.DataFrame(index = np.arange(len(seeds)), columns = ['seed','ctime','ttime'])
data['seed'] = seeds
data['ctime'] = ctimes
data['ttime'] = ttimes
data.to_csv('%s/fluctuations.txt'%(gpath))