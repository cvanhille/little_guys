import numpy as np
import sys
import argparse
import os
import glob
from tqdm.auto import tqdm

parser = argparse.ArgumentParser(description='')
parser.add_argument('-p', '--path', help='path to simulation set - REQUIRED', required=True, type=str)
parser.add_argument('-Kbond', '--Kbond', help='bond constant [kT/sigma2]', required=False, type=float, default=500.0)
parser.add_argument('-eps', '--eps', help='binding constant [kT]', required=False, type=float, default=0.0)
parser.add_argument('-epsA', '--epsA', help='binding constant for intra line pair [kT]', required=False, type=float, default=0.0)
parser.add_argument('-epsB', '--epsB', help='binding constant for intra perp pair [kT]', required=False, type=float, default=0.0)
parser.add_argument('-coff', '--coff', help='cutoff distance [sigma]', required=False, type=float, default=1.2)
parser.add_argument('-runtime', '--runtime', help='simulation run time [simulation time units]', required=False, type=float, default=2000.0)
parser.add_argument('-frate', '--frate', help='frame rate [simulation time units]', required=False, type=float, default=1.0)
parser.add_argument('-N','--N', help='number of repetitions', required=False, type=int, default=1000)
parser.add_argument('-config','--config', help='configuration to simulate - REQUIRED - recommend 2MOLS_C-B_R', required=True, type=str)
parser.add_argument('-bonding','--bonding', help='bonding mode', required=False, action = 'store_true')
parser.add_argument('-bonds','--bonds', help='number of binding bonds [5X, 3X, 3P, 3L]', required=False, choices=['5X','3X','3P','3L'], default='3P')
parser.add_argument('-repside', '--repside', help='EV size for 2-2, 3-3, 5-5 and 6-6 interactions [sigma]', required=False, type=float, default=1.0)
parser.add_argument('-capside','--capside', help='cap the side of the molecule', required=False, action = 'store_true')
parser.add_argument('-phi', '--phi', help='area fraction', required=False, type=float, default=0.0)
parser.add_argument('-L', '--L', help='box size', required=False, type=float, default=0.0)

args = parser.parse_args()
gpath = args.path
Kbond = float(args.Kbond)
eps = float(args.eps)
epsA = float(args.epsA)
epsB = float(args.epsB)
coff = float(args.coff)
runtime = float(args.runtime)
frate = float(args.frate)
N = int(args.N)
config = args.config
bonding = args.bonding
nbonds = args.bonds
repside = float(args.repside)
capside = args.capside
phi = float(args.phi)
L = float(args.L)

if config == 'NMOLS_FREEZE' or config == 'NMOLS':
	if phi == 0 or L == 0:
		print()
		print("Error! Missing L and phi arguments! Aborting!!")
		print()
		exit()

seeds = np.random.randint(0,9000, size = N) + 1000

while len(np.unique(seeds)) < N:
	seeds0 = np.unique(seeds)
	seeds1 = np.random.randint(0,9000, size = N-len(seeds0)) + 1000
	seeds = np.concatenate((seeds0,seeds1))

if os.access(gpath, os.F_OK):
	print()
	print("Careful! This path already exists!")
	print("Make sure this is the set you want to create. If so, remove the existing folder by running 'rm -r %s'"%(gpath))
	print()
	exit()

r = os.system('mkdir %s'%(gpath))

for i, seed in enumerate(tqdm(seeds)):
	command = 'python3 make_5B_files.py -p %s/sd%d -config %s -runtime %f -frate %f -Kbond %f -eps %f -epsA %f -epsB %f -coff %f -repside %f -phi %f -L %f -sd %d'%(gpath,seed,config,runtime,frate,Kbond,eps,epsA,epsB,coff,repside,phi,L,seed)
	if bonding:
		command = '%s -bonding -bonds %s'%(command,nbonds)
	if capside:
		command = '%s -capside'%(command)
	# print(command)
	r = os.system(command)

print("Done!")
