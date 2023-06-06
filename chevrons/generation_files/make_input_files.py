import numpy as np
import sys
import argparse
import os
import glob

parser = argparse.ArgumentParser(description='')
parser.add_argument('-p', '--path', help='path to simulation folder - REQUIRED', required=True, type=str)
parser.add_argument('-Kbond', '--Kbond', help='bond constant [kT/sigma2]', required=False, type=float, default=500.0)
parser.add_argument('-eps', '--eps', help='binding constant [kT]', required=False, type=float, default=8.0)
parser.add_argument('-epsA', '--epsA', help='binding constant for intra line pair [kT]', required=False, type=float, default=0.0)
parser.add_argument('-epsB', '--epsB', help='binding constant for intra perp pair [kT]', required=False, type=float, default=0.0)
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
epsA = float(args.epsA)
epsB = float(args.epsB)
coff = float(args.coff)
tstep = float(args.tstep)
runtime = float(args.runtime)
frate = float(args.frate)
seed = int(args.seed)
runsteps = int(runtime/tstep)
dump = int(frate/tstep)
config = args.config

np.random.seed(seed)

if not os.access('%s/in.local'%(gpath), os.F_OK):
	r = os.system('mkdir %s'%(gpath))
else:
	print()
	print("Careful! This simulation path already exists and contains simulation files. Double check that you want to write this file.")
	print("If you do, remove the path with 'rm -r %s' and rerun the generator"%(gpath))
	print()
	exit()

if config.split('_')[0] == 'FAR':
	print()
	print("Generating new config file for FAR setup of binding interactions...")
	if config.split('_')[1] == 'LB':
		print("Left-binding setup")
		f = open('%s/configuration.dat'%(gpath), 'w')
		f.write("# Configuration file for the initial conditions of the chevron-like molecules simulations\n\n")
		f.write("12 atoms\n")
		f.write("17 bonds\n")
		f.write("12 atom types\n")
		f.write("1 bond types\n")
		f.write("-10.0 10.0 xlo xhi\n")
		f.write("-10.0 10.0 ylo yhi\n")
		f.write("-0.25 0.25 zlo zhi\n\n")
		f.write("Masses\n\n")
		f.write("1 1\n")
		f.write("2 1\n")
		f.write("3 1\n")
		f.write("4 1\n")
		f.write("5 1\n")
		f.write("6 1\n")
		f.write("7 1\n")
		f.write("8 1\n")
		f.write("9 1\n")
		f.write("10 1\n")
		f.write("11 1\n")
		f.write("12 1\n\n")
		f.write("Atoms\n\n")
		pos0x = -10*np.random.random()
		pos0y = 20*(0.5-np.random.random())
		pos1x = 10*np.random.random()
		pos1y = 20*(0.5-np.random.random())
		# f.write("1 1 1 %f %f 0.0000000\n"%(pos0x,pos0y))
		# f.write("2 1 2 %f %f 0.0000000\n"%(pos0x+0.5,pos0y+0.8660254))
		# f.write("3 1 3 %f %f 0.0000000\n"%(pos0x-0.5,pos0y+0.8660254))
		# f.write("4 1 4 %f %f 0.0000000\n"%(pos0x-1.0,pos0y))
		# f.write("5 1 5 %f %f 0.0000000\n"%(pos0x-0.5,pos0y-0.8660254))
		# f.write("6 1 6 %f %f 0.0000000\n"%(pos0x+0.5,pos0y-0.8660254))
		# f.write("7 2 7 %f %f 0.0000000\n"%(pos1x,pos1y))
		# f.write("8 2 8 %f %f 0.0000000\n"%(pos1x,pos1y+1.0))
		# f.write("9 2 9 %f %f 0.0000000\n"%(pos1x-0.8660254,pos1y+0.5))
		# f.write("10 2 10 %f %f 0.0000000\n"%(pos1x-1.73205081,pos1y))
		# f.write("11 2 11 %f %f 0.0000000\n"%(pos1x-0.8660254,pos1y-0.5))
		# f.write("12 2 12 %f %f 0.0000000\n\n"%(pos1x,pos1y-1.0))
		f.write("1 1 1 %f %f 0.0000000\n"%(pos0x,pos0y))
		f.write("2 1 2 %f %f 0.0000000\n"%(pos0x,pos0y+1.0))
		f.write("3 1 3 %f %f 0.0000000\n"%(pos0x-0.8660254,pos0y+0.5))
		f.write("4 1 4 %f %f 0.0000000\n"%(pos0x-1.73205081,pos0y))
		f.write("5 1 5 %f %f 0.0000000\n"%(pos0x-0.8660254,pos0y-0.5))
		f.write("6 1 6 %f %f 0.0000000\n"%(pos0x,pos0y-1.0))
		f.write("7 2 7 %f %f 0.0000000\n"%(pos1x,pos1y))
		f.write("8 2 8 %f %f 0.0000000\n"%(pos1x+0.5,pos1y+0.8660254))
		f.write("9 2 9 %f %f 0.0000000\n"%(pos1x-0.5,pos1y+0.8660254))
		f.write("10 2 10 %f %f 0.0000000\n"%(pos1x-1.0,pos1y))
		f.write("11 2 11 %f %f 0.0000000\n"%(pos1x-0.5,pos1y-0.8660254))
		f.write("12 2 12 %f %f 0.0000000\n\n"%(pos1x+0.5,pos1y-0.8660254))
		f.write("Bonds\n\n")
		f.write("1 1 1 2\n")
		f.write("2 1 1 3\n")
		f.write("3 1 1 5\n")
		f.write("4 1 1 6\n")
		f.write("5 1 2 3\n")
		f.write("6 1 3 4\n")
		f.write("7 1 4 5\n")
		f.write("8 1 5 6\n")
		f.write("9 1 7 8\n")
		f.write("10 1 7 9\n")
		f.write("11 1 7 11\n")
		f.write("12 1 7 12\n")
		f.write("13 1 8 9\n")
		f.write("14 1 9 10\n")
		f.write("15 1 10 11\n")
		f.write("16 1 11 12\n")
		f.write("17 1 7 10\n")
		f.close()
	if config.split('_')[1] == 'RB':
		print("Right-binding setup")
		f = open('%s/configuration.dat'%(gpath), 'w')
		f.write("# Configuration file for the initial conditions of the chevron-like molecules simulations\n\n")
		f.write("12 atoms\n")
		f.write("17 bonds\n")
		f.write("12 atom types\n")
		f.write("1 bond types\n")
		f.write("-20.0 20.0 xlo xhi\n")
		f.write("-20.0 20.0 ylo yhi\n")
		f.write("-0.25 0.25 zlo zhi\n\n")
		f.write("Masses\n\n")
		f.write("1 1\n")
		f.write("2 1\n")
		f.write("3 1\n")
		f.write("4 1\n")
		f.write("5 1\n")
		f.write("6 1\n")
		f.write("7 1\n")
		f.write("8 1\n")
		f.write("9 1\n")
		f.write("10 1\n")
		f.write("11 1\n")
		f.write("12 1\n\n")
		f.write("Atoms\n\n")
		pos0x = -10*np.random.random()
		pos0y = 20*(0.5-np.random.random())
		pos1x = 10*np.random.random()
		pos1y = 20*(0.5-np.random.random())
		if np.sqrt((pos0x-pos1x)**2+(pos0y-pos1y)**2) < 2.0:
			pos0x -= 4.0
		# f.write("1 1 1 %f %f 0.0000000\n"%(pos0x,pos0y))
		# f.write("2 1 2 %f %f 0.0000000\n"%(pos0x,pos0y+1.0))
		# f.write("3 1 3 %f %f 0.0000000\n"%(pos0x-0.8660254,pos0y+0.5))
		# f.write("4 1 4 %f %f 0.0000000\n"%(pos0x-1.73205081,pos0y))
		# f.write("5 1 5 %f %f 0.0000000\n"%(pos0x-0.8660254,pos0y-0.5))
		# f.write("6 1 6 %f %f 0.0000000\n"%(pos0x,pos0y-1.0))
		# f.write("7 2 7 %f %f 0.0000000\n"%(pos1x,pos1y))
		# f.write("8 2 8 %f %f 0.0000000\n"%(pos1x+0.5,pos1y+0.8660254))
		# f.write("9 2 9 %f %f 0.0000000\n"%(pos1x-0.5,pos1y+0.8660254))
		# f.write("10 2 10 %f %f 0.0000000\n"%(pos1x-1.0,pos1y))
		# f.write("11 2 11 %f %f 0.0000000\n"%(pos1x-0.5,pos1y-0.8660254))
		# f.write("12 2 12 %f %f 0.0000000\n\n"%(pos1x+0.5,pos1y-0.8660254))
		f.write("1 1 1 %f %f 0.0000000\n"%(pos0x,pos0y))
		f.write("2 1 2 %f %f 0.0000000\n"%(pos0x+0.5,pos0y+0.8660254))
		f.write("3 1 3 %f %f 0.0000000\n"%(pos0x-0.5,pos0y+0.8660254))
		f.write("4 1 4 %f %f 0.0000000\n"%(pos0x-1.0,pos0y))
		f.write("5 1 5 %f %f 0.0000000\n"%(pos0x-0.5,pos0y-0.8660254))
		f.write("6 1 6 %f %f 0.0000000\n"%(pos0x+0.5,pos0y-0.8660254))
		f.write("7 2 7 %f %f 0.0000000\n"%(pos1x,pos1y))
		f.write("8 2 8 %f %f 0.0000000\n"%(pos1x,pos1y+1.0))
		f.write("9 2 9 %f %f 0.0000000\n"%(pos1x-0.8660254,pos1y+0.5))
		f.write("10 2 10 %f %f 0.0000000\n"%(pos1x-1.73205081,pos1y))
		f.write("11 2 11 %f %f 0.0000000\n"%(pos1x-0.8660254,pos1y-0.5))
		f.write("12 2 12 %f %f 0.0000000\n\n"%(pos1x,pos1y-1.0))
		f.write("Bonds\n\n")
		f.write("1 1 1 2\n")
		f.write("2 1 1 3\n")
		f.write("3 1 1 5\n")
		f.write("4 1 1 6\n")
		f.write("5 1 2 3\n")
		f.write("6 1 3 4\n")
		f.write("7 1 4 5\n")
		f.write("8 1 5 6\n")
		f.write("9 1 7 8\n")
		f.write("10 1 7 9\n")
		f.write("11 1 7 11\n")
		f.write("12 1 7 12\n")
		f.write("13 1 8 9\n")
		f.write("14 1 9 10\n")
		f.write("15 1 10 11\n")
		f.write("16 1 11 12\n")
		f.write("17 1 1 4\n")
		f.close()
	print()
else:
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
variable            epsA equal %.1f                           	# binding constant A pair [kT]
variable            epsB equal %.1f                           	# binding constant B pair [kT]
variable            coff equal %f                           	# cutoff distance [sigma]
variable            tstep equal %f                             	# simulation timestep size [seconds]
variable            seed equal %d                               	# random number generator seed
variable            run_steps equal %d          					# simulation run time [simulation steps]
variable            dump_time equal %d        						# dumping interval [simulation steps]
'''%(Kbond,eps,epsA,epsB,coff,tstep,seed,runsteps,dump))
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
	f.write("pair_coeff          2 9 cosine/squared ${eps} 1.00 ${coff} wca 			# top interaction 2-3\n")
	f.write("pair_coeff          6 11 cosine/squared ${eps} 1.00 ${coff} wca 			# bottom interaction 6-5\n")
	if config.split('_')[-1] == 'FR' or config.split('_')[0] == 'FAR':
		f.write("pair_coeff          1 10 cosine/squared ${eps} 1.00 ${coff} wca 			# middle interaction 1-4\n")
	if epsA > 0:
		if config.split('_')[1] == 'LB':
			f.write("pair_coeff          1 4 cosine/squared ${epsA} 1.00 ${coff} wca 			# A pair intra 1-4\n")
		elif config.split('_')[1] == 'RB':
			f.write("pair_coeff          7 10 cosine/squared ${epsA} 1.00 ${coff} wca 			# A pair intra 1-4\n")
	if epsB > 0:
		if config.split('_')[1] == 'LB':
			f.write("pair_coeff          3 5 cosine/squared ${epsB} 1.00 ${coff} wca 			# B pair intra 3-5\n")
		elif config.split('_')[1] == 'RB':
			f.write("pair_coeff          9 11 cosine/squared ${epsB} 1.00 ${coff} wca 			# B pair intra 3-5\n")
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