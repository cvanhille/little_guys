import numpy as np
import sys
import argparse
import os
import glob

parser = argparse.ArgumentParser(description='')
parser.add_argument('-p', '--path', help='path to simulation folder - REQUIRED', required=True, type=str)
parser.add_argument('-Kbond', '--Kbond', help='bond constant [kT/sigma2]', required=False, type=float, default=500.0)
parser.add_argument('-eps', '--eps', help='binding constant [kT]', required=False, type=float, default=8.0)
parser.add_argument('-coff', '--coff', help='cutoff distance [sigma]', required=False, type=float, default=1.2)
parser.add_argument('-tstep', '--tstep', help='simulation timestep [simulation time units]', required=False, type=float, default=0.01)
parser.add_argument('-runtime', '--runtime', help='simulation run time [simulation time units]', required=False, type=float, default=100.0)
parser.add_argument('-frate', '--frate', help='frame rate [simulation time units]', required=False, type=float, default=0.1)
parser.add_argument('-sd','--seed', help='random number generator seed', required=False, type=int, default=1234)
parser.add_argument('-config','--config', help='configuration to simulate - REQUIRED', required=True, type=str)

args = parser.parse_args()
gpath = args.path
Kbond = float(args.Kbond)
eps = float(args.eps)
coff = float(args.coff)
tstep = float(args.tstep)
runtime = float(args.runtime)
frate = float(args.frate)
seed = int(args.seed)
runsteps = int(runtime/tstep)
dump = int(frate/tstep)
config = args.config

if not os.access('%s/in.local'%(gpath), os.F_OK):
	r = os.system('mkdir %s'%(gpath))
else:
	print()
	print("Careful! This simulation path already exists and contains simulation files. Double check that you want to write this file.")
	print("If you do, remove the path with 'rm -r %s' and rerun the generator"%(gpath))
	print()
	exit()

if os.access('config_files/%s.config'%(config), os.F_OK):
	r = os.system('cp config_files/%s.config %s/configuration.dat'%(config,gpath))
else:
	print()
	print("Error! Wrong choice of configuration, the file does not exist. Please try again and make sure you choose an existing file.")
	print("If the file you need does not exist you may need to write it yourself...")
	print("Available configurations:")
	cfs = glob.glob('config_files/*.config')
	for f in cfs:
		print('\t%s'%(f.split('config_files/')[1].split('.config')[0]))
	print()
	exit()

f = open('%s/info.txt'%(gpath), 'w')
f.write("path:\t\t\t%s\n"%(gpath))
f.write("Kbond [kT/sigma2]:\t%.1f\n"%(Kbond))
f.write("eps [kT]:\t\t%.1f\n"%(eps))
f.write("coff [sigma]:\t\t%f\n"%(coff))
f.write("tstep [tau]:\t\t%.5f\n"%(tstep))
f.write("runtime [tau]:\t\t%.1f - %d steps\n"%(runtime,runsteps))
f.write("frate [tau]:\t\t%.1f - %d steps\n"%(frate,dump))
f.write("seed:\t\t\t%d\n"%(seed))
f.write("configuration:\t\t%s\n"%(config))
f.close()

f = open('%s/in.local'%(gpath), 'w')
f.write('''# Input script for LAMMPS simulation of individual chevron-like molecules fluctuating between two configurations

units               lj
atom_style          molecular
dimension           2 
boundary            p p p
log                 log.txt
read_data           configuration.dat

variable            Kbond equal %.1f                           	# bond constant [kT/sigma2]
variable            eps equal %.1f                           	# binding constant [kT]
variable            coff equal %f                           	# cutoff distance [sigma]
variable            tstep equal %f                             	# simulation timestep size [seconds]
variable            seed equal %d                               	# random number generator seed
variable            run_steps equal %d          					# simulation run time [simulation steps]
variable            dump_time equal %d        						# dumping interval [simulation steps]
'''%(Kbond,eps,coff,tstep,seed,runsteps,dump))
f.write('''
special_bonds       lj 1.0 1.0 1.0
bond_style          harmonic
bond_coeff          1 ${Kbond} 1.0

pair_style          hybrid/overlay zero 1.50 cosine/squared 1.50
pair_coeff          * * cosine/squared 0.00 1.00 1.10
pair_coeff          * * zero 1.50
pair_coeff          * * cosine/squared 1.00 1.00 1.00 wca
''')
if eps > 0:
	f.write("pair_coeff          2 9 cosine/squared ${eps} 1.00 ${coff} wca 			# top interaction (horizontal) 2-3\n")
	f.write("pair_coeff          6 11 cosine/squared ${eps} 1.00 ${coff} wca 			# bottom interaction (vertical) 6-5\n")
f.write('''
fix                 fLang all langevin 1.0 1.0 1.0 ${seed}
fix                 fNVE all nve

dump                1 all custom ${dump_time} output.xyz id mol type x y
dump_modify         1 format line "%d %d %d %.2f %.2f"

thermo              ${dump_time}
thermo_style        custom step temp pe ke etotal epair ebond press vol density atoms

compute             cBonds all property/local batom1 batom2 btype
compute             cBondDxys all bond/local engpot force dist

dump                2 all local ${dump_time} bonds.dump c_cBonds[*] c_cBondDxys[*]
dump_modify         2 format line "%f %f %f %.2f %.2f %.2f"

fix                 twodim all enforce2d

timestep            ${tstep}
run                 ${run_steps}
''')
f.close()