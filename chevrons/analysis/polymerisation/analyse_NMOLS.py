import numpy as np
import pandas as pd
from ovito.io import *
from ovito.modifiers import *
import sys
import glob
import os
from tqdm.auto import tqdm
import matplotlib.pyplot as plt
import time

def chains_in_frame(p, f):
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
    tbs = tps[top]
    mbs = mls[top]
    top = ids[top]
    L = d.cell[0,0]
    L2 = L/2.0

    imbs = []
    for i in range(len(tbs)):
        mb = mbs[i]
        tb = tbs[i]
        if mb[0] == mb[1]:
            continue
        if tb[0] == 1 and tb[1] == 4:
            imbs.append([tb[0],tb[1],mb[0],mb[1],1,0])
        elif tb[0] == 4 and tb[1] == 1:
            imbs.append([tb[1],tb[0],mb[1],mb[0],-1,0])
        elif tb[0] == 2 and tb[1] == 3:
            imbs.append([tb[0],tb[1],mb[0],mb[1],1,1])
        elif tb[0] == 3 and tb[1] == 2:
            imbs.append([tb[1],tb[0],mb[1],mb[0],-1,1])
        elif tb[0] == 6 and tb[1] == 5:
            imbs.append([tb[0],tb[1],mb[0],mb[1],1,-1])
        elif tb[0] == 5 and tb[1] == 6:
            imbs.append([tb[1],tb[0],mb[1],mb[0],-1,-1])
    imbs = np.array(imbs)
    
    if len(imbs) == 0:
        return [], []

    m0 = imbs[:,2]
    m1 = imbs[:,3]
    um0 = np.unique(m0)
    cm0 = np.zeros(len(um0))
    um1 = um0.copy()
    for i in range(len(um0)):
        cm0[i] = len(m0[m0 == um0[i]])
        um1[i] = np.min(m1[m0 == um0[i]])

    fum0 = um0[cm0 >= 3]
    fum1 = um1[cm0 >= 3]
    fcm0 = cm0[cm0 >= 3]

    chains = []
    lengths = []
    skips = []
    for i in range(len(fum0)):
        if skips.count(i) > 0:
            continue
        ch = np.array([fum0[i],fum1[i]])
        w = np.where(fum0 == fum1[i])[0]
        while len(w) > 0:
            ni = w[0]
            ch = np.concatenate((ch,[fum1[ni]]))
            skips.append(ni)
            w = np.where(fum0 == fum1[ni])[0]
        w = np.where(fum1 == fum0[i])[0]
        while len(w) > 0:
            ni = w[0]
            ch = np.concatenate(([fum0[ni]],ch))
            skips.append(ni)
            w = np.where(fum1 == fum0[ni])[0]
        chains.append(ch)
        lengths.append(len(ch))

    lengths = np.array(lengths)
    
    return chains, lengths

def analyse(path, gpath):

    tic = time.perf_counter()

    spath = path.split(gpath)[1]
    print(spath)

    file = '%s/output.xyz'%(path)
    bonds = '%s/bonds.dump'%(path)
    p = import_file(file)
    mod = LoadTrajectoryModifier()
    p.modifiers.append(mod)
    mod.source.load(bonds, columns = ['Particle Identifiers.1', 'Particle Identifiers.2', 'Bond Type'], multiple_frames = True)
    maxf = p.source.num_frames

    pols = []
    mons = []
    mpls = []
    frms = []
    plen = []

    heads = 0
    tails = 0
    merges = 0
    nucs = 0

    hts = [0]
    tts = [0]
    mts = [0]
    nts = [0]
    mhts = [0]
    mtts = [0]
    mmts = [0]
    mnts = [0]

    fl = open('%s/polymerisation_analysis_events/%s.txt'%(gpath,spath), 'w')

    for pf in range(maxf):
        f = pf+1

        chains, lengths = chains_in_frame(p, f)
        pchains, plengths = chains_in_frame(p, pf)

        f = f*10

        frms.append(f)
        pols.append(len(lengths))
        pms = np.sum(lengths)
        mons.append(80-pms)
        mpls.append(pms)
        if len(lengths) > 0:
            plen.append(lengths.mean())
        else:
            plen.append(0)

        if len(lengths) > 0:
            
            # Nucleation
            idx = np.where(lengths == 2)[0]
            pidx = np.where(plengths == 2)[0]
            if len(idx) > 0:
                for i in idx:
                    c = chains[i]
                    exists = 0
                    if len(pidx) > 0:
                        for j in pidx:
                            pc = pchains[j]
                            if pc[0] == c[0] and pc[1] == c[1]:
                                exists = 1
                    else:
                        exists = 0
                    if exists == 0:
                        nucs += 1
                        nts.append(f)
                        mnts.append(80-pms)
                        fl.write('Nucleation at time %.1f - ... -> ['%(f))
                        for ci in c:
                            fl.write('%d '%(ci))
                        fl.write(']\n')
            
            # Polymerisation and merging
            idx = np.where(lengths > 2)[0]
            if len(idx) > 0:
                for i in idx:
                    c = chains[i]
                    interacts = []
                    types = []
                    for pc in pchains:
                        if pc[0] == c[0]:
                            if pc[-1] != c[-1]:
                                interacts.append(pc)
                                types.append(-1)
                            else:
                                interacts.append(0)
                                types.append(0)
                        elif pc[-1] == c[-1]:
                                interacts.append(pc)
                                types.append(1)
                    if len(interacts) == 1:
                        if types[0] == 1:
                            heads += 1
                            hts.append(f)
                            mhts.append(80-pms)
                            fl.write('Head growth at time %.1f - '%(f))
                            for ci in interacts[0]:
                                fl.write('%d '%(ci))
                            fl.write('] -> [')
                            for ci in c:
                                fl.write('%d '%(ci))
                            fl.write(']\n')
                        elif types[0] == -1:
                            tails += 1
                            tts.append(f)
                            mtts.append(80-pms)
                            fl.write('Tail growth at time %.1f - '%(f))
                            for ci in interacts[0]:
                                fl.write('%d '%(ci))
                            fl.write('] -> [')
                            for ci in c:
                                fl.write('%d '%(ci))
                            fl.write(']\n')
                    elif len(interacts) == 2:
                        merges += 1
                        mts.append(f)
                        mmts.append(80-pms)
                        fl.write('Merging at time %.1f - '%(f))
                        for ci in interacts[0]:
                            fl.write('%d '%(ci))
                        fl.write('] + [')
                        for ci in interacts[1]:
                            fl.write('%d '%(ci))
                        fl.write('] -> [')
                        for ci in c:
                            fl.write('%d '%(ci))
                        fl.write(']\n')
                    elif len(interacts) > 2:
                        fl.write('Too many chains merging! - time %.1f\n'%(f))
                    elif len(interacts) == 0:
                        nucs += 1
                        nts.append(f)
                        mnts.append(80-pms)
                        fl.write('Big nucleation (size %d) at time %.1f - ... -> ['%(len(c),f))
                        for ci in c:
                            fl.write('%d '%(ci))
                        fl.write(']\n')
    
    fl.close()


    pols = np.array(pols)
    mons = np.array(mons)
    mpls = np.array(mpls)
    frms = np.array(frms)
    plen = np.array(plen)

    data = pd.DataFrame(index = frms, columns = ['solmons','polmons','polymers','avgsize'])
    data['solmons'] = mons
    data['polmons'] = mpls
    data['polymers'] = pols
    data['avgsize'] = plen
    data.to_csv('%s/polymerisation_analysis_timeseries/%s.txt'%(gpath,spath))

    hts = np.array(hts)
    tts = np.array(tts)
    mts = np.array(mts)
    nts = np.array(nts)

    mhts = np.array(mhts)
    mtts = np.array(mtts)
    mmts = np.array(mmts)
    mnts = np.array(mnts)

    data = pd.DataFrame(index = np.arange(np.max([len(hts),len(tts),len(mts),len(nts)])), columns = ['head_times','tail_times','merge_times','nuc_times','head_mols','tail_mols','merge_mols','nuc_mols'])
    data['head_times'].iloc[:len(hts)] = hts
    data['tail_times'].iloc[:len(tts)] = tts
    data['merge_times'].iloc[:len(mts)] = mts
    data['nuc_times'].iloc[:len(nts)] = nts
    data['head_mols'].iloc[:len(mhts)] = mhts
    data['tail_mols'].iloc[:len(mtts)] = mtts
    data['merge_mols'].iloc[:len(mmts)] = mmts
    data['nuc_mols'].iloc[:len(mnts)] = mnts
    data.to_csv('%s/polymerisation_analysis_results/%s.txt'%(gpath,spath))

    print()
    print("%d head polymerisation events in %d tau -- %.4f events/tau"%(heads, f, heads/f))
    print("%d tail polymerisation events in %d tau -- %.4f events/tau"%(tails, f, tails/f))
    print("%d merging events in %d tau -- %.4f events/tau"%(merges, f, merges/f))
    print("%d nucleation events in %d tau -- %.4f events/tau -- %.2f tau on average"%(nucs, f, nucs/f, np.nanmean(nts[1:])))
    print("   ", heads/tails, "H/T ratio")
    print("   ", heads/nucs, "H/N ratio")
    print("   ", heads/merges, "H/M ratio")
    print("   ", tails/merges, "T/M ratio")

    toc = time.perf_counter()
    print()
    print(toc-tic, 'seconds elapsed')
    
    return

if len(sys.argv) < 2:
    print()
    print("Error! Wrong number of arguments. Please provide: 1) general path to sample")
    print()
    exit()

gpath = sys.argv[1]

files = glob.glob('%s/sd*/output.xyz'%(gpath))

r = os.system('mkdir %s/polymerisation_analysis_events'%(gpath))
r = os.system('mkdir %s/polymerisation_analysis_timeseries'%(gpath))
r = os.system('mkdir %s/polymerisation_analysis_results'%(gpath))

for i, file in enumerate(tqdm(files)):
    path = file.split('/output.xyz')[0]
    analyse(path, gpath)
