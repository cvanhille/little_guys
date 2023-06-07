import numpy as np
import sys
import argparse
import os
import glob
from tqdm.auto import tqdm

parser = argparse.ArgumentParser(description='')
parser.add_argument('-p', '--path', help='path to simulation set - REQUIRED', required=True, type=str)
parser.add_argument('-Kbond', '--Kbond', help='bond constant [kT/sigma2]', required=False, type=float, default=500.0)
parser.add_argument('-eps', '--eps', help='binding constant [kT]', required=False, type=float, default=8.0)
parser.add_argument('-epsA', '--epsA', help='binding constant for intra line pair [kT]', required=False, type=float, default=0.0)
parser.add_argument('-epsB', '--epsB', help='binding constant for intra perp pair [kT]', required=False, type=float, default=0.0)
parser.add_argument('-runtime', '--runtime', help='simulation run time [simulation time units]', required=False, type=float, default=10000.0)
parser.add_argument('-frate', '--frate', help='frame rate [simulation time units]', required=False, type=float, default=1.0)
parser.add_argument('-N','--N', help='number of repetitions', required=False, type=int, default=1000)
parser.add_argument('-config','--config', help='configuration to simulate - REQUIRED', required=True, type=str)

args = parser.parse_args()
gpath = args.path
Kbond = float(args.Kbond)
eps = float(args.eps)
epsA = float(args.epsA)
epsB = float(args.epsB)
runtime = float(args.runtime)
frate = float(args.frate)
N = int(args.N)
config = args.config

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
	command = 'python3 make_input_files.py -p %s/sd%d -config %s -runtime %f -frate %f -Kbond %f -eps %f -epsA %f -epsB %f -sd %d'%(gpath,seed,config,runtime,frate,Kbond,eps,epsA,epsB,seed)
	r = os.system(command)

print("Done!")
