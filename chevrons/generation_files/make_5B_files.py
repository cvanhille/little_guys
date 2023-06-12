import numpy as np
import sys
import argparse
import os
import glob

def shape(X0,Y0,key,L):
	Xs = [X0]
	Ys = [Y0]
	ids = np.array([1,2,3,4,5,6])
	types = np.array([1,2,3,4,5,6])
	if key == 'RNDM' or key == 'R' or key == 'RANDOM':
		rs = np.random.random()
	elif key == 'CAT' or key == 'CHEVRON' or key == 'C':
		rs = 0.0
	elif key == 'TRIANGLE' or key == 'TRI' or key == 'T':
		rs = 1.0
	rd = np.random.randint(6)
	if rs < 0.5:
		# Cat
		# Bead 2
		XN = X0+np.cos(rd*np.pi/3)
		YN = Y0+np.sin(rd*np.pi/3)
		if XN >= L/2.0:
			XN -= L
		if XN < -L/2.0:
			XN += L
		if YN >= L/2.0:
			YN -= L
		if YN < -L/2.0:
			YN += L
		Ys.append(YN)
		Xs.append(XN)
		# Bead 3
		XN = X0+np.cos((rd+1)*np.pi/3)
		YN = Y0+np.sin((rd+1)*np.pi/3)
		if XN >= L/2.0:
			XN -= L
		if XN < -L/2.0:
			XN += L
		if YN >= L/2.0:
			YN -= L
		if YN < -L/2.0:
			YN += L
		Ys.append(YN)
		Xs.append(XN)
		# Bead 4
		XN = X0+np.cos((rd+2)*np.pi/3)
		YN = Y0+np.sin((rd+2)*np.pi/3)
		if XN >= L/2.0:
			XN -= L
		if XN < -L/2.0:
			XN += L
		if YN >= L/2.0:
			YN -= L
		if YN < -L/2.0:
			YN += L
		Ys.append(YN)
		Xs.append(XN)
		# Bead 5
		XN = X0+np.cos((rd+3)*np.pi/3)
		YN = Y0+np.sin((rd+3)*np.pi/3)
		if XN >= L/2.0:
			XN -= L
		if XN < -L/2.0:
			XN += L
		if YN >= L/2.0:
			YN -= L
		if YN < -L/2.0:
			YN += L
		Ys.append(YN)
		Xs.append(XN)
		# Bead 6
		XN = X0+np.cos((rd+4)*np.pi/3)
		YN = Y0+np.sin((rd+4)*np.pi/3)
		if XN >= L/2.0:
			XN -= L
		if XN < -L/2.0:
			XN += L
		if YN >= L/2.0:
			YN -= L
		if YN < -L/2.0:
			YN += L
		Ys.append(YN)
		Xs.append(XN)
	else:
		# Triangle
		# Bead 2
		XN = X0+np.cos(rd*np.pi/3)
		YN = Y0+np.sin(rd*np.pi/3)
		if XN >= L/2.0:
			XN -= L
		if XN < -L/2.0:
			XN += L
		if YN >= L/2.0:
			YN -= L
		if YN < -L/2.0:
			YN += L
		Ys.append(YN)
		Xs.append(XN)
		# Bead 3
		XN = X0+np.cos((rd+1)*np.pi/3)
		YN = Y0+np.sin((rd+1)*np.pi/3)
		if XN >= L/2.0:
			XN -= L
		if XN < -L/2.0:
			XN += L
		if YN >= L/2.0:
			YN -= L
		if YN < -L/2.0:
			YN += L
		Ys.append(YN)
		Xs.append(XN)
		# Bead 4
		XN = X0+np.cos((rd+1)*np.pi/3)+np.cos((rd+2)*np.pi/3)
		YN = Y0+np.sin((rd+1)*np.pi/3)+np.sin((rd+2)*np.pi/3)
		if XN >= L/2.0:
			XN -= L
		if XN < -L/2.0:
			XN += L
		if YN >= L/2.0:
			YN -= L
		if YN < -L/2.0:
			YN += L
		Ys.append(YN)
		Xs.append(XN)
		# Bead 5
		XN = X0+np.cos((rd+2)*np.pi/3)
		YN = Y0+np.sin((rd+2)*np.pi/3)
		if XN >= L/2.0:
			XN -= L
		if XN < -L/2.0:
			XN += L
		if YN >= L/2.0:
			YN -= L
		if YN < -L/2.0:
			YN += L
		Ys.append(YN)
		Xs.append(XN)
		# Bead 6
		XN = X0+np.cos((rd+3)*np.pi/3)
		YN = Y0+np.sin((rd+3)*np.pi/3)
		if XN >= L/2.0:
			XN -= L
		if XN < -L/2.0:
			XN += L
		if YN >= L/2.0:
			YN -= L
		if YN < -L/2.0:
			YN += L
		Ys.append(YN)
		Xs.append(XN)
	Xs = np.array(Xs)
	Ys = np.array(Ys)
	return Xs, Ys, ids, types

def setup_2mols(L, key1, key2):
	x0 = -L/2.0
	y0 = -L/2.0
	dx = 1.0
	dy = 2*np.cos(np.pi/6)
	ddx = np.sin(np.pi/6)
	ddy = np.cos(np.pi/6)
	nrows = int(np.floor(L/dy))
	ncols = int(np.floor(L/dx))
	if L-ncols*dx < ddx:
		ncols -= 1
	xs = []
	ys = []
	for nr in range(nrows):
		for nc in range(ncols):
			xs.append(x0+nc*dx)
			ys.append(y0+nr*dy)
			xs.append(x0+nc*dx+ddx)
			ys.append(y0+nr*dy+ddy)
	xs = np.array(xs)
	ys = np.array(ys)
	xs += (L/2.0-xs.max())/2.0
	ys += (L/2.0-ys.max())/2.0
	# Molecule 1
	idx = np.random.randint(len(xs))
	X0 = xs[idx]
	Y0 = ys[idx]
	Xs, Ys, ids, types = shape(X0,Y0,key1,L)
	# Empty safe space
	dxs = xs-X0
	dys = ys-Y0
	dxs[dxs < -L/2.0] += L
	dxs[dxs >= L/2.0] -= L
	dys[dys < -L/2.0] += L
	dys[dys >= L/2.0] -= L
	dr2 = dxs*dxs+dys*dys
	nxs = xs[dr2 > 9.0]
	nys = ys[dr2 > 9.0]
	# Molecule 2
	idx = np.random.randint(len(nxs))
	nX0 = nxs[idx]
	nY0 = nys[idx]
	nXs, nYs, nids, ntypes = shape(nX0,nY0,key2,L)
	# Assemble
	fids = np.concatenate((ids,nids+ids.max()))
	ftypes = np.concatenate((types,ntypes))
	fX = np.concatenate((Xs,nXs))
	fY = np.concatenate((Ys,nYs))
	return fids, ftypes, fX, fY

parser = argparse.ArgumentParser(description='')
parser.add_argument('-p', '--path', help='path to simulation folder - REQUIRED', required=True, type=str)
parser.add_argument('-Kbond', '--Kbond', help='bond constant [kT/sigma2]', required=False, type=float, default=500.0)
parser.add_argument('-eps', '--eps', help='binding constant [kT]', required=False, type=float, default=0.0)
parser.add_argument('-epsA', '--epsA', help='binding constant for intra line pair [kT]', required=False, type=float, default=0.0)
parser.add_argument('-epsB', '--epsB', help='binding constant for intra perp pair [kT]', required=False, type=float, default=0.0)
parser.add_argument('-coff', '--coff', help='cutoff distance [sigma]', required=False, type=float, default=1.2)
parser.add_argument('-tstep', '--tstep', help='simulation timestep [simulation time units]', required=False, type=float, default=0.01)
parser.add_argument('-runtime', '--runtime', help='simulation run time [simulation time units]', required=False, type=float, default=100.0)
parser.add_argument('-frate', '--frate', help='frame rate [simulation time units]', required=False, type=float, default=0.1)
parser.add_argument('-sd','--seed', help='random number generator seed', required=False, type=int, default=1234)
parser.add_argument('-config','--config', help='configuration to simulate - REQUIRED', required=True, type=str)
parser.add_argument('-bonding','--bonding', help='bonding mode', required=False, action = 'store_true')

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
bonding = args.bonding

np.random.seed(seed)

if not os.access('%s/in.local'%(gpath), os.F_OK):
	r = os.system('mkdir %s'%(gpath))
else:
	print()
	print("Careful! This simulation path already exists and contains simulation files. Double check that you want to write this file.")
	print("If you do, remove the path with 'rm -r %s' and rerun the generator"%(gpath))
	print()
	exit()

if config.split('_')[0] == '2MOLS':
	ids, types, X, Y = setup_2mols(20, config.split('_')[1].split('-')[0], config.split('_')[2])
	if not len(ids) == 12:
		print()
		print("Error! Wrong number of atoms generated by the setup function... Check code please!")
		print()
		exit()
	mols = np.ones(len(X))
	mols[6:] = 2
	f = open('%s/configuration.dat'%(gpath), 'w')
	f.write("# Configuration file for the initial conditions of the chevron-like molecules simulations\n\n")
	f.write("12 atoms\n")
	if config.split('_')[1].split('-')[1] == 'B':
		f.write("17 bonds\n")
	elif config.split('_')[1].split('-')[1] == 'F':
		f.write("16 bonds\n")
	f.write("6 atom types\n")
	f.write("2 bond types\n")
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
	f.write("Atoms\n\n")
	for i in range(len(X)):
		f.write("%d %d %d %f %f 0.0000000\n"%(ids[i],mols[i],types[i],X[i],Y[i]))
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
	if config.split('_')[1].split('-')[1] == 'B':
		f.write("17 1 1 4\n")
	f.close()
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
f.write("epsA [kT]:\t\t%.1f\n"%(epsA))
f.write("epsB [kT]:\t\t%.1f\n"%(epsB))
f.write("coff [sigma]:\t\t%f\n"%(coff))
f.write("tstep [tau]:\t\t%.5f\n"%(tstep))
f.write("runtime [tau]:\t\t%.1f - %d steps\n"%(runtime,runsteps))
f.write("frate [tau]:\t\t%.1f - %d steps\n"%(frate,dump))
f.write("seed:\t\t\t%d\n"%(seed))
f.write("configuration:\t\t%s\n"%(config))
if bonding:
	f.write("bonding mode:\t\tTrue\n")
else:
	f.write("bonding mode:\t\tFalse\n")
f.close()

f = open('%s/in.local'%(gpath), 'w')
f.write('''# Input script for LAMMPS simulation of individual chevron-like molecules fluctuating between two configurations

units               lj
atom_style          molecular
dimension           2 
boundary            p p p
log                 log.txt
read_data           configuration.dat''')
if bonding:
	f.write("  extra/bond/per/atom 5  extra/special/per/atom 20  extra/angle/per/atom 3\n\n")
f.write('''variable            Kbond equal %.1f                           	# bond constant [kT/sigma2]
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
bond_coeff          2 ${Kbond} 1.7320508075688774

pair_style          hybrid/overlay zero 2.00 cosine/squared 2.00
pair_coeff          * * cosine/squared 0.00 1.00 1.10
pair_coeff          * * zero 2.00
pair_coeff          * * cosine/squared 1.00 1.00 1.00 wca
''')
if eps > 0:
	f.write("pair_coeff          2 9 cosine/squared ${eps} 1.00 ${coff} wca 			# top interaction 2-3\n")
	f.write("pair_coeff          6 11 cosine/squared ${eps} 1.00 ${coff} wca 			# bottom interaction 6-5\n")
	if config.split('_')[-1] == 'FR' or config.split('_')[0] == 'FAR':
		f.write("pair_coeff          1 10 cosine/squared ${eps} 1.00 ${coff} wca 			# middle interaction 1-4\n")
if epsA > 0:
	if (len(config.split('_')) > 1 and config.split('_')[1] == 'LB') or config == 'TETRAMER' or config == 'DECAMER':
		f.write("pair_coeff          1 4 cosine/squared ${epsA} 1.00 ${coff} wca 			# A pair intra 1-4\n")
	elif (len(config.split('_')) > 1 and config.split('_')[1] == 'RB'):
		f.write("pair_coeff          7 10 cosine/squared ${epsA} 1.00 ${coff} wca 			# A pair intra 1-4\n")
	elif (len(config.split('_')) > 1 and config.split('_')[1] == 'LBF') or (len(config.split('_')) > 1 and config.split('_')[1] == 'RBF'):
		f.write("pair_coeff          1 4 cosine/squared ${epsA} 1.00 ${coff} wca 			# A pair intra 1-4\n")
		f.write("pair_coeff          7 10 cosine/squared ${epsA} 1.00 ${coff} wca 			# A pair intra 1-4\n")
if epsB > 0:
	if (len(config.split('_')) > 1 and config.split('_')[1] == 'LB') or config == 'TETRAMER' or config == 'DECAMER':
		f.write("pair_coeff          3 5 cosine/squared ${epsB} 1.00 ${coff} wca 			# B pair intra 3-5\n")
	elif (len(config.split('_')) > 1 and config.split('_')[1] == 'RB'):
		f.write("pair_coeff          9 11 cosine/squared ${epsB} 1.00 ${coff} wca 			# B pair intra 3-5\n")
	elif (len(config.split('_')) > 1 and config.split('_')[1] == 'LBF') or (len(config.split('_')) > 1 and config.split('_')[1] == 'RBF'):
		f.write("pair_coeff          3 5 cosine/squared ${epsB} 1.00 ${coff} wca 			# B pair intra 3-5\n")
		f.write("pair_coeff          9 11 cosine/squared ${epsB} 1.00 ${coff} wca 			# B pair intra 3-5\n")
if bonding:
	f.write("\n")
	f.write("fix					fBindC all bond/create 1 1 4 1.05 1 inter_mol\n")
	f.write("fix					fBindTS all bond/create 1 2 4 1.05 1 inter_mol\n")
	f.write("fix					fBindBS all bond/create 1 6 4 1.05 1 inter_mol\n")
	f.write("fix					fBindTL all bond/create 1 1 3 1.78 2 inter_mol\n")
	f.write("fix					fBindBL all bond/create 1 1 5 1.78 2 inter_mol\n")
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
