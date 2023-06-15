import numpy as np
import pandas as pd
from ovito.io import *
from ovito.modifiers import *
import sys
import glob
from tqdm.auto import tqdm

def check_frame_5X(parts, top):

    # Bonds that can form

    # RB:
    # 1-10 -- straight
    # 1-11 -- long diagonal bottom
    # 1-9 -- long diagonal top
    # 6-10 -- short diagonal bottom
    # 2-10 -- short diagonal top

    RBds = []
    RBbs = []

    r = parts[parts[:,0] == 1][:,3:] - parts[parts[:,0] == 10][:,3:]
    r[r >= 10] -= 20
    r[r < -10] += 20
    r = np.linalg.norm(r)
    RBds.append(r)

    r = parts[parts[:,0] == 1][:,3:] - parts[parts[:,0] == 11][:,3:]
    r[r >= 10] -= 20
    r[r < -10] += 20
    r = np.linalg.norm(r)
    RBds.append(r)

    r = parts[parts[:,0] == 1][:,3:] - parts[parts[:,0] == 9][:,3:]
    r[r >= 10] -= 20
    r[r < -10] += 20
    r = np.linalg.norm(r)
    RBds.append(r)

    r = parts[parts[:,0] == 6][:,3:] - parts[parts[:,0] == 10][:,3:]
    r[r >= 10] -= 20
    r[r < -10] += 20
    r = np.linalg.norm(r)
    RBds.append(r)

    r = parts[parts[:,0] == 2][:,3:] - parts[parts[:,0] == 10][:,3:]
    r[r >= 10] -= 20
    r[r < -10] += 20
    r = np.linalg.norm(r)
    RBds.append(r)

    RBds = np.array(RBds)

    b = 0
    ts = top[top[:,0] == 1]
    tss = ts[ts[:,1] == 10]
    if len(tss) > 0:
        b = 1
    ts = top[top[:,1] == 1]
    tss = ts[ts[:,0] == 10]
    if len(tss) > 0:
        b = 1
    RBbs.append(b)

    b = 0
    ts = top[top[:,0] == 1]
    tss = ts[ts[:,1] == 11]
    if len(tss) > 0:
        b = 1
    ts = top[top[:,1] == 1]
    tss = ts[ts[:,0] == 11]
    if len(tss) > 0:
        b = 1
    RBbs.append(b)

    b = 0
    ts = top[top[:,0] == 1]
    tss = ts[ts[:,1] == 9]
    if len(tss) > 0:
        b = 1
    ts = top[top[:,1] == 1]
    tss = ts[ts[:,0] == 9]
    if len(tss) > 0:
        b = 1
    RBbs.append(b)

    b = 0
    ts = top[top[:,0] == 6]
    tss = ts[ts[:,1] == 10]
    if len(tss) > 0:
        b = 1
    ts = top[top[:,1] == 6]
    tss = ts[ts[:,0] == 10]
    if len(tss) > 0:
        b = 1
    RBbs.append(b)

    b = 0
    ts = top[top[:,0] == 2]
    tss = ts[ts[:,1] == 10]
    if len(tss) > 0:
        b = 1
    ts = top[top[:,1] == 2]
    tss = ts[ts[:,0] == 10]
    if len(tss) > 0:
        b = 1
    RBbs.append(b)

    RBbs = np.array(RBbs)

    # LB:
    # 7-4 -- straight
    # 7-5 -- long diagonal bottom
    # 7-3 -- long diagonal top
    # 12-4 -- short diagonal bottom
    # 8-4 -- short diagonal top
    
    LBds = []
    LBbs = []

    r = parts[parts[:,0] == 7][:,3:] - parts[parts[:,0] == 4][:,3:]
    r[r >= 10] -= 20
    r[r < -10] += 20
    r = np.linalg.norm(r)
    LBds.append(r)

    r = parts[parts[:,0] == 7][:,3:] - parts[parts[:,0] == 5][:,3:]
    r[r >= 10] -= 20
    r[r < -10] += 20
    r = np.linalg.norm(r)
    LBds.append(r)

    r = parts[parts[:,0] == 7][:,3:] - parts[parts[:,0] == 3][:,3:]
    r[r >= 10] -= 20
    r[r < -10] += 20
    r = np.linalg.norm(r)
    LBds.append(r)

    r = parts[parts[:,0] == 12][:,3:] - parts[parts[:,0] == 4][:,3:]
    r[r >= 10] -= 20
    r[r < -10] += 20
    r = np.linalg.norm(r)
    LBds.append(r)

    r = parts[parts[:,0] == 8][:,3:] - parts[parts[:,0] == 4][:,3:]
    r[r >= 10] -= 20
    r[r < -10] += 20
    r = np.linalg.norm(r)
    LBds.append(r)

    LBds = np.array(LBds)

    b = 0
    ts = top[top[:,0] == 7]
    tss = ts[ts[:,1] == 4]
    if len(tss) > 0:
        b = 1
    ts = top[top[:,1] == 7]
    tss = ts[ts[:,0] == 4]
    if len(tss) > 0:
        b = 1
    LBbs.append(b)

    b = 0
    ts = top[top[:,0] == 7]
    tss = ts[ts[:,1] == 5]
    if len(tss) > 0:
        b = 1
    ts = top[top[:,1] == 7]
    tss = ts[ts[:,0] == 5]
    if len(tss) > 0:
        b = 1
    LBbs.append(b)

    b = 0
    ts = top[top[:,0] == 7]
    tss = ts[ts[:,1] == 3]
    if len(tss) > 0:
        b = 1
    ts = top[top[:,1] == 7]
    tss = ts[ts[:,0] == 3]
    if len(tss) > 0:
        b = 1
    LBbs.append(b)

    b = 0
    ts = top[top[:,0] == 12]
    tss = ts[ts[:,1] == 4]
    if len(tss) > 0:
        b = 1
    ts = top[top[:,1] == 12]
    tss = ts[ts[:,0] == 4]
    if len(tss) > 0:
        b = 1
    LBbs.append(b)

    b = 0
    ts = top[top[:,0] == 8]
    tss = ts[ts[:,1] == 4]
    if len(tss) > 0:
        b = 1
    ts = top[top[:,1] == 8]
    tss = ts[ts[:,0] == 4]
    if len(tss) > 0:
        b = 1
    LBbs.append(b)

    LBbs = np.array(LBbs)
    
    return LBds, LBbs, RBds, RBbs

def check_frame_3X(parts, top):

    # Bonds that can form

    # RB:
    # 1-10 -- straight
    # 6-10 -- short diagonal bottom
    # 2-10 -- short diagonal top

    RBds = []
    RBbs = []

    r = parts[parts[:,0] == 1][:,3:] - parts[parts[:,0] == 10][:,3:]
    r[r >= 10] -= 20
    r[r < -10] += 20
    r = np.linalg.norm(r)
    RBds.append(r)

    r = parts[parts[:,0] == 6][:,3:] - parts[parts[:,0] == 10][:,3:]
    r[r >= 10] -= 20
    r[r < -10] += 20
    r = np.linalg.norm(r)
    RBds.append(r)

    r = parts[parts[:,0] == 2][:,3:] - parts[parts[:,0] == 10][:,3:]
    r[r >= 10] -= 20
    r[r < -10] += 20
    r = np.linalg.norm(r)
    RBds.append(r)

    RBds = np.array(RBds)

    b = 0
    ts = top[top[:,0] == 1]
    tss = ts[ts[:,1] == 10]
    if len(tss) > 0:
        b = 1
    ts = top[top[:,1] == 1]
    tss = ts[ts[:,0] == 10]
    if len(tss) > 0:
        b = 1
    RBbs.append(b)

    b = 0
    ts = top[top[:,0] == 6]
    tss = ts[ts[:,1] == 10]
    if len(tss) > 0:
        b = 1
    ts = top[top[:,1] == 6]
    tss = ts[ts[:,0] == 10]
    if len(tss) > 0:
        b = 1
    RBbs.append(b)

    b = 0
    ts = top[top[:,0] == 2]
    tss = ts[ts[:,1] == 10]
    if len(tss) > 0:
        b = 1
    ts = top[top[:,1] == 2]
    tss = ts[ts[:,0] == 10]
    if len(tss) > 0:
        b = 1
    RBbs.append(b)

    RBbs = np.array(RBbs)

    # LB:
    # 7-4 -- straight
    # 12-4 -- short diagonal bottom
    # 8-4 -- short diagonal top
    
    LBds = []
    LBbs = []

    r = parts[parts[:,0] == 7][:,3:] - parts[parts[:,0] == 4][:,3:]
    r[r >= 10] -= 20
    r[r < -10] += 20
    r = np.linalg.norm(r)
    LBds.append(r)

    r = parts[parts[:,0] == 12][:,3:] - parts[parts[:,0] == 4][:,3:]
    r[r >= 10] -= 20
    r[r < -10] += 20
    r = np.linalg.norm(r)
    LBds.append(r)

    r = parts[parts[:,0] == 8][:,3:] - parts[parts[:,0] == 4][:,3:]
    r[r >= 10] -= 20
    r[r < -10] += 20
    r = np.linalg.norm(r)
    LBds.append(r)

    LBds = np.array(LBds)

    b = 0
    ts = top[top[:,0] == 7]
    tss = ts[ts[:,1] == 4]
    if len(tss) > 0:
        b = 1
    ts = top[top[:,1] == 7]
    tss = ts[ts[:,0] == 4]
    if len(tss) > 0:
        b = 1
    LBbs.append(b)

    b = 0
    ts = top[top[:,0] == 12]
    tss = ts[ts[:,1] == 4]
    if len(tss) > 0:
        b = 1
    ts = top[top[:,1] == 12]
    tss = ts[ts[:,0] == 4]
    if len(tss) > 0:
        b = 1
    LBbs.append(b)

    b = 0
    ts = top[top[:,0] == 8]
    tss = ts[ts[:,1] == 4]
    if len(tss) > 0:
        b = 1
    ts = top[top[:,1] == 8]
    tss = ts[ts[:,0] == 4]
    if len(tss) > 0:
        b = 1
    LBbs.append(b)

    LBbs = np.array(LBbs)
    
    return LBds, LBbs, RBds, RBbs

def check_frame_3P(parts, top):

    # Bonds that can form

    # RB:
    # 1-10 -- straight middle
    # 2-9 -- straight top
    # 6-11 -- straight bottom

    RBds = []
    RBbs = []

    r = parts[parts[:,0] == 1][:,3:] - parts[parts[:,0] == 10][:,3:]
    r[r >= 10] -= 20
    r[r < -10] += 20
    r = np.linalg.norm(r)
    RBds.append(r)

    r = parts[parts[:,0] == 2][:,3:] - parts[parts[:,0] == 9][:,3:]
    r[r >= 10] -= 20
    r[r < -10] += 20
    r = np.linalg.norm(r)
    RBds.append(r)

    r = parts[parts[:,0] == 6][:,3:] - parts[parts[:,0] == 11][:,3:]
    r[r >= 10] -= 20
    r[r < -10] += 20
    r = np.linalg.norm(r)
    RBds.append(r)

    RBds = np.array(RBds)

    b = 0
    ts = top[top[:,0] == 1]
    tss = ts[ts[:,1] == 10]
    if len(tss) > 0:
        b = 1
    ts = top[top[:,1] == 1]
    tss = ts[ts[:,0] == 10]
    if len(tss) > 0:
        b = 1
    RBbs.append(b)

    b = 0
    ts = top[top[:,0] == 2]
    tss = ts[ts[:,1] == 9]
    if len(tss) > 0:
        b = 1
    ts = top[top[:,1] == 2]
    tss = ts[ts[:,0] == 9]
    if len(tss) > 0:
        b = 1
    RBbs.append(b)

    b = 0
    ts = top[top[:,0] == 6]
    tss = ts[ts[:,1] == 11]
    if len(tss) > 0:
        b = 1
    ts = top[top[:,1] == 6]
    tss = ts[ts[:,0] == 11]
    if len(tss) > 0:
        b = 1
    RBbs.append(b)

    RBbs = np.array(RBbs)

    # LB:
    # 7-4 -- straight middle
    # 8-3 -- straight diagonal top
    # 12-5 -- straight diagonal bottom
    
    LBds = []
    LBbs = []

    r = parts[parts[:,0] == 7][:,3:] - parts[parts[:,0] == 4][:,3:]
    r[r >= 10] -= 20
    r[r < -10] += 20
    r = np.linalg.norm(r)
    LBds.append(r)

    r = parts[parts[:,0] == 8][:,3:] - parts[parts[:,0] == 3][:,3:]
    r[r >= 10] -= 20
    r[r < -10] += 20
    r = np.linalg.norm(r)
    LBds.append(r)

    r = parts[parts[:,0] == 12][:,3:] - parts[parts[:,0] == 5][:,3:]
    r[r >= 10] -= 20
    r[r < -10] += 20
    r = np.linalg.norm(r)
    LBds.append(r)

    LBds = np.array(LBds)

    b = 0
    ts = top[top[:,0] == 7]
    tss = ts[ts[:,1] == 4]
    if len(tss) > 0:
        b = 1
    ts = top[top[:,1] == 7]
    tss = ts[ts[:,0] == 4]
    if len(tss) > 0:
        b = 1
    LBbs.append(b)

    b = 0
    ts = top[top[:,0] == 8]
    tss = ts[ts[:,1] == 3]
    if len(tss) > 0:
        b = 1
    ts = top[top[:,1] == 8]
    tss = ts[ts[:,0] == 3]
    if len(tss) > 0:
        b = 1
    LBbs.append(b)

    b = 0
    ts = top[top[:,0] == 12]
    tss = ts[ts[:,1] == 5]
    if len(tss) > 0:
        b = 1
    ts = top[top[:,1] == 12]
    tss = ts[ts[:,0] == 5]
    if len(tss) > 0:
        b = 1
    LBbs.append(b)

    LBbs = np.array(LBbs)
    
    return LBds, LBbs, RBds, RBbs

def check_frame_CR(parts, top):

    # Bonds that can form

    # RB:
    # 1-12 -- straight middle
    # 2-11 -- straight top
    # 6-13 -- straight bottom

    RBds = []
    RBbs = []

    r = parts[parts[:,0] == 1][:,3:] - parts[parts[:,0] == 12][:,3:]
    r[r >= 10] -= 20
    r[r < -10] += 20
    r = np.linalg.norm(r)
    RBds.append(r)

    r = parts[parts[:,0] == 2][:,3:] - parts[parts[:,0] == 11][:,3:]
    r[r >= 10] -= 20
    r[r < -10] += 20
    r = np.linalg.norm(r)
    RBds.append(r)

    r = parts[parts[:,0] == 6][:,3:] - parts[parts[:,0] == 13][:,3:]
    r[r >= 10] -= 20
    r[r < -10] += 20
    r = np.linalg.norm(r)
    RBds.append(r)

    RBds = np.array(RBds)

    b = 0
    ts = top[top[:,0] == 1]
    tss = ts[ts[:,1] == 12]
    if len(tss) > 0:
        b = 1
    ts = top[top[:,1] == 1]
    tss = ts[ts[:,0] == 12]
    if len(tss) > 0:
        b = 1
    RBbs.append(b)

    b = 0
    ts = top[top[:,0] == 2]
    tss = ts[ts[:,1] == 11]
    if len(tss) > 0:
        b = 1
    ts = top[top[:,1] == 2]
    tss = ts[ts[:,0] == 11]
    if len(tss) > 0:
        b = 1
    RBbs.append(b)

    b = 0
    ts = top[top[:,0] == 6]
    tss = ts[ts[:,1] == 13]
    if len(tss) > 0:
        b = 1
    ts = top[top[:,1] == 6]
    tss = ts[ts[:,0] == 13]
    if len(tss) > 0:
        b = 1
    RBbs.append(b)

    RBbs = np.array(RBbs)

    # LB:
    # 9-4 -- straight middle
    # 10-3 -- straight diagonal top
    # 14-5 -- straight diagonal bottom
    
    LBds = []
    LBbs = []

    r = parts[parts[:,0] == 9][:,3:] - parts[parts[:,0] == 4][:,3:]
    r[r >= 10] -= 20
    r[r < -10] += 20
    r = np.linalg.norm(r)
    LBds.append(r)

    r = parts[parts[:,0] == 10][:,3:] - parts[parts[:,0] == 3][:,3:]
    r[r >= 10] -= 20
    r[r < -10] += 20
    r = np.linalg.norm(r)
    LBds.append(r)

    r = parts[parts[:,0] == 14][:,3:] - parts[parts[:,0] == 5][:,3:]
    r[r >= 10] -= 20
    r[r < -10] += 20
    r = np.linalg.norm(r)
    LBds.append(r)

    LBds = np.array(LBds)

    b = 0
    ts = top[top[:,0] == 9]
    tss = ts[ts[:,1] == 4]
    if len(tss) > 0:
        b = 1
    ts = top[top[:,1] == 9]
    tss = ts[ts[:,0] == 4]
    if len(tss) > 0:
        b = 1
    LBbs.append(b)

    b = 0
    ts = top[top[:,0] == 10]
    tss = ts[ts[:,1] == 3]
    if len(tss) > 0:
        b = 1
    ts = top[top[:,1] == 10]
    tss = ts[ts[:,0] == 3]
    if len(tss) > 0:
        b = 1
    LBbs.append(b)

    b = 0
    ts = top[top[:,0] == 14]
    tss = ts[ts[:,1] == 5]
    if len(tss) > 0:
        b = 1
    ts = top[top[:,1] == 14]
    tss = ts[ts[:,0] == 5]
    if len(tss) > 0:
        b = 1
    LBbs.append(b)

    LBbs = np.array(LBbs)
    
    return LBds, LBbs, RBds, RBbs

def times(path,mode):
    file = '%s/output.xyz'%(path)
    bonds = '%s/bonds.dump'%(path)
    p = import_file(file)
    mod = LoadTrajectoryModifier()
    p.modifiers.append(mod)
    mod.source.load(bonds, columns = ['Particle Identifiers.1', 'Particle Identifiers.2', 'Bond Type', 'Energy.1', 'Force.1', 'Length'], multiple_frames = True)
    # Run through frames and compute bonds and distances
    # bondsL = []
    # bondsR = []
    # distsL = []
    # distsR = []
    frames = []
    nbndsL = []
    nbndsR = []
    # for i, f in enumerate(tqdm(range(p.source.num_frames+1))):
    for f in range(p.source.num_frames+1):
        d = p.compute(f)
        pos = np.array(d.particles['Position'])
        ids = np.array(d.particles['Particle Identifier'])
        tps = np.array(d.particles['Particle Type'])
        mls = np.array(d.particles['Molecule Identifier'])
        top = np.array(d.particles.bonds['Topology'])
        parts = np.empty((len(ids),5))
        parts[:,0] = ids
        parts[:,1] = tps
        parts[:,2] = mls
        parts[:,3] = pos[:,0]
        parts[:,4] = pos[:,1]
        top = ids[top]
        # Analyse frame
        if mode == '5X':
            LBds, LBbs, RBds, RBbs = check_frame_5X(parts, top)
        elif mode == '3X':
            LBds, LBbs, RBds, RBbs = check_frame_3X(parts, top)
        elif mode == '3P' or mode == '3L':
            LBds, LBbs, RBds, RBbs = check_frame_3P(parts, top)
        elif mode == 'CR':
            LBds, LBbs, RBds, RBbs = check_frame_CR(parts, top)
        # bondsL.append(LBbs)
        # bondsR.append(RBbs)
        # distsL.append(LBds)
        # distsR.append(RBds)
        frames.append(f)
        nbndsL.append(len(LBbs[LBbs>0]))
        nbndsR.append(len(RBbs[RBbs>0]))        
    # bondsL = np.array(bondsL)
    # bondsR = np.array(bondsR)
    # distsL = np.array(distsL)
    # distsR = np.array(distsR)
    frames = np.array(frames)
    nbndsL = np.array(nbndsL)
    nbndsR = np.array(nbndsR)
    # Find times for left binding
    Lts = []
    if len(frames[nbndsL>=1]) == 0:
        Lts.append(np.nan)
    else:
        Lts.append(np.min(frames[nbndsL>=1]))
    if len(frames[nbndsL>=2]) == 0:
        Lts.append(np.nan)
    else:
        Lts.append(np.min(frames[nbndsL>=2]))
    if len(frames[nbndsL>=3]) == 0:
        Lts.append(np.nan)
    else:
        Lts.append(np.min(frames[nbndsL>=3]))
    if len(frames[nbndsL>=4]) == 0:
        Lts.append(np.nan)
    else:
        Lts.append(np.min(frames[nbndsL>=4]))
    if len(frames[nbndsL>=5]) == 0:
        Lts.append(np.nan)
    else:
        Lts.append(np.min(frames[nbndsL>=5]))
    # Find times for right binding
    Rts = []
    if len(frames[nbndsR>=1]) == 0:
        Rts.append(np.nan)
    else:
        Rts.append(np.min(frames[nbndsR>=1]))
    if len(frames[nbndsR>=2]) == 0:
        Rts.append(np.nan)
    else:
        Rts.append(np.min(frames[nbndsR>=2]))
    if len(frames[nbndsR>=3]) == 0:
        Rts.append(np.nan)
    else:
        Rts.append(np.min(frames[nbndsR>=3]))
    if len(frames[nbndsR>=4]) == 0:
        Rts.append(np.nan)
    else:
        Rts.append(np.min(frames[nbndsR>=4]))
    if len(frames[nbndsR>=5]) == 0:
        Rts.append(np.nan)
    else:
        Rts.append(np.min(frames[nbndsR>=5]))
    return np.array(Lts), np.array(Rts)

if len(sys.argv) < 3:
    print()
    print("Error! Wrong number of arguments. Please provide: 1) general path to sample 2) bond binding mode [5X,3X,3P,3L,CR] 3) [y/n] to do the analysis no matter the number of seeds --optional")
    print()
    exit()

gpath = sys.argv[1]
mode = sys.argv[2]
if not (mode == '5X' or mode == '3X' or mode == '3P' or mode == '3L' or mode == 'CR'):
    print()
    print("Error! Wrong choice of bond binding mode. Pick from [5X,3X,3P,3L,CR]")
    print()
    exit()

files = glob.glob('%s/sd*/output.xyz'%(gpath))

if len(sys.argv) == 4:
    ck = sys.argv[3]
else:
    print()
    ck = input("%d simulations to analyse! Continue? [y/n] "%(len(files)))
    print()

if not (ck == 'y' or ck == 'yes' or ck == 'Y'):
    print("Ok! Aborting...")
    print()
    exit()

seeds = []
timeL1 = []
timeL2 = []
timeL3 = []
timeL4 = []
timeL5 = []
timeR1 = []
timeR2 = []
timeR3 = []
timeR4 = []
timeR5 = []
for i, file in enumerate(tqdm(files)):
    path = file.split('/output.xyz')[0]
    seed = int(path.split('/sd')[1])
    LTs, RTs = times(path,mode)
    seeds.append(seed)
    timeL1.append(LTs[0])
    timeL2.append(LTs[1])
    timeL3.append(LTs[2])
    timeL4.append(LTs[3])
    timeL5.append(LTs[4])
    timeR1.append(RTs[0])
    timeR2.append(RTs[1])
    timeR3.append(RTs[2])
    timeR4.append(RTs[3])
    timeR5.append(RTs[4])

data = pd.DataFrame(index = np.arange(len(seeds)), columns = ['seed', 'L1', 'L2', 'L3', 'L4', 'L5', 'R1', 'R2', 'R3', 'R4', 'R5'])
data['seed'] = seeds
data['L1'] = timeL1
data['L2'] = timeL2
data['L3'] = timeL3
data['L4'] = timeL4
data['L5'] = timeL5
data['R1'] = timeR1
data['R2'] = timeR2
data['R3'] = timeR3
data['R4'] = timeR4
data['R5'] = timeR5
data.to_csv('%s/bonding_times.txt'%(gpath))
