import numpy as np
import sys
import argparse
import os
import glob

def shape(X0,Y0,key,L,capside,orkey):
	Xs = [X0]
	Ys = [Y0]
	if capside:
		ids = np.array([1,2,3,4,5,6,7,8])
		types = np.array([1,2,3,4,5,6,7,7])
	else:
		ids = np.array([1,2,3,4,5,6])
		types = np.array([1,2,3,4,5,6])
	if key == 'RNDM' or key == 'R' or key == 'RANDOM':
		rs = np.random.random()
	elif key == 'CAT' or key == 'CHEVRON' or key == 'C':
		rs = 0.0
	elif key == 'TRIANGLE' or key == 'TRI' or key == 'T':
		rs = 1.0
	if orkey == 'RNDM' or orkey == 'R' or orkey == 'RANDOM':
		rd = np.random.randint(6)
	else:
		rd = int(orkey)
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
		if capside:
			# Bead 7
			XN = X0+np.cos(rd*np.pi/3)+np.cos((rd+1)*np.pi/3)
			YN = Y0+np.sin(rd*np.pi/3)+np.sin((rd+1)*np.pi/3)
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
			# Bead 8
			XN = X0+np.cos((rd+4)*np.pi/3)+np.cos((rd+3)*np.pi/3)
			YN = Y0+np.sin((rd+4)*np.pi/3)+np.sin((rd+3)*np.pi/3)
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
		if capside:
			# Bead 7
			XN = X0+np.cos(rd*np.pi/3)+np.cos((rd+1)*np.pi/3)
			YN = Y0+np.sin(rd*np.pi/3)+np.sin((rd+1)*np.pi/3)
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
			# Bead 8
			XN = X0+np.cos((rd+3)*np.pi/3)+np.cos((rd+2)*np.pi/3)
			YN = Y0+np.sin((rd+3)*np.pi/3)+np.sin((rd+2)*np.pi/3)
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

def setup_2mols(L, key1, key2, capside):
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
	Xs, Ys, ids, types = shape(X0,Y0,key1,L,capside,'R')
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
	nXs, nYs, nids, ntypes = shape(nX0,nY0,key2,L,capside,'R')
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
parser.add_argument('-bcoff', '--bcoff', help='bound label cutoff distance [sigma]', required=False, type=float, default=1.1)
parser.add_argument('-tstep', '--tstep', help='simulation timestep [simulation time units]', required=False, type=float, default=0.01)
parser.add_argument('-runtime', '--runtime', help='simulation run time [simulation time units]', required=False, type=float, default=10000.0)
parser.add_argument('-frate', '--frate', help='frame rate [simulation time units]', required=False, type=float, default=100.0)
parser.add_argument('-sd','--seed', help='random number generator seed', required=False, type=int, default=1234)
parser.add_argument('-L', '--L', help='box size', required=False, type=float, default=20.0)

args = parser.parse_args()
gpath = args.path
Kbond = float(args.Kbond)
eps = float(args.eps)
epsA = float(args.epsA)
epsB = float(args.epsB)
coff = float(args.coff)
bcoff = float(args.bcoff)
tstep = float(args.tstep)
runtime = float(args.runtime)
frate = float(args.frate)
seed = int(args.seed)
runsteps = int(runtime/tstep)
dump = int(frate/tstep)
L = float(args.L)

np.random.seed(seed)

if not os.access('%s/in0.local'%(gpath), os.F_OK):
	r = os.system('mkdir %s'%(gpath))
else:
	print()
	print("Careful! This simulation path already exists and contains simulation files. Double check that you want to write this file.")
	print("If you do, remove the path with 'rm -r %s' and rerun the generator"%(gpath))
	print()
	exit()

# Setup initial conditions: 2 molecules - one is randomly chevron or triangle and one is a chevron with one more bond to stabilise it
ids, types, X, Y = setup_2mols(L, 'C', 'R', True)
if not len(ids) == 16:
	print()
	print("Error! Wrong number of atoms generated by the setup function... Check code please!")
	print()
	exit()
mols = np.ones(len(X))
mols[8:] = 2
f = open('%s/configuration.dat'%(gpath), 'w')
f.write("# Configuration file for the initial conditions of the chevron-like molecules simulations\n\n")
f.write("16 atoms\n")
f.write("25 bonds\n")
f.write("7 atom types\n")
f.write("2 bond types\n")
f.write("-%.1f %.1f xlo xhi\n"%(-L/2.0,L/2.0))
f.write("-%.1f %.1f ylo yhi\n"%(-L/2.0,L/2.0))
f.write("-0.25 0.25 zlo zhi\n\n")
f.write("Masses\n\n")
f.write("1 1\n")
f.write("2 1\n")
f.write("3 1\n")
f.write("4 1\n")
f.write("5 1\n")
f.write("6 1\n")
f.write("7 1\n\n")
f.write("Atoms\n\n")
for i in range(len(X)):
	f.write("%d %d %d %f %f 0.0\n"%(ids[i],mols[i],types[i],X[i],Y[i]))
f.write('\n')
f.write("Bonds\n\n")
f.write("1 1 1 2\n")
f.write("2 1 1 3\n")
f.write("3 1 1 5\n")
f.write("4 1 1 6\n")
f.write("5 1 2 3\n")
f.write("6 1 3 4\n")
f.write("7 1 4 5\n")
f.write("8 1 5 6\n")
f.write("9 1 2 7\n")
f.write("10 1 3 7\n")
f.write("11 1 5 8\n")
f.write("12 1 6 8\n")
f.write("13 1 9 10\n")
f.write("14 1 9 11\n")
f.write("15 1 9 13\n")
f.write("16 1 9 14\n")
f.write("17 1 10 11\n")
f.write("18 1 11 12\n")
f.write("19 1 12 13\n")
f.write("20 1 13 14\n")
f.write("21 1 10 15\n")
f.write("22 1 11 15\n")
f.write("23 1 13 16\n")
f.write("24 1 14 16\n")
f.write("25 1 1 4\n")
f.close()

f = open('%s/info.txt'%(gpath), 'w')
f.write("path:\t\t\t%s\n"%(gpath))
f.write("L [sigma]:\t\t%.2f\n"%(L))
f.write("Kbond [kT/sigma2]:\t%.1f\n"%(Kbond))
f.write("eps [kT]:\t\t%.1f\n"%(eps))
f.write("epsA [kT]:\t\t%.1f\n"%(epsA))
f.write("epsB [kT]:\t\t%.1f\n"%(epsB))
f.write("coff [sigma]:\t\t%f\n"%(coff))
f.write("bcoff [sigma]:\t\t%f\n"%(bcoff))
f.write("tstep [tau]:\t\t%.5f\n"%(tstep))
f.write("runtime [tau]:\t\t%.1f - %d steps\n"%(runtime,runsteps))
f.write("frate [tau]:\t\t%.1f - %d steps\n"%(frate,dump))
f.write("seed:\t\t\t%d\n"%(seed))
f.write("configuration:\t\t2MOLS_C-B_R\n")
f.write("Caps on side:\t\tTrue\n")
f.close()

f = open('%s/in0.local'%(gpath), 'w')
f.write('''# Input script for LAMMPS simulation of binding two chevron/triangle molecules

units               lj
atom_style          molecular
dimension           2 
boundary            p p p
log                 log0.txt
read_data           configuration.dat

group               moving id > 0

variable            Kbond equal %.1f                           	# bond constant [kT/sigma2]
variable            eps equal %.1f                           	# binding constant [kT]
variable            epsA equal %.1f                           	# binding constant A pair [kT]
variable            epsB equal %.1f                           	# binding constant B pair [kT]
variable            coff equal %f                           	# cutoff distance [sigma]
variable            bcoff equal %f                           	# cutoff distance [sigma]
variable            rdist equal 1.0                           	# cutoff distance [sigma]
variable            tstep equal %f                             	# simulation timestep size [seconds]
variable            seed equal %d                               	# random number generator seed
variable            run_steps equal %d          					# simulation run time [simulation steps]
variable            dump_steps equal %d        						# dumping interval [simulation steps]
'''%(Kbond,eps,epsA,epsB,coff,bcoff,tstep,seed,runsteps,dump))
f.write('''
special_bonds       lj 0.0 1.0 1.0
bond_style          harmonic
bond_coeff          1 ${Kbond} 1.0
bond_coeff          2 ${Kbond} 1.7320508075688774

pair_style          hybrid/overlay zero 2.00 cosine/squared 2.00
pair_coeff          * * cosine/squared 0.00 1.00 1.10
pair_coeff          * * zero 2.00
pair_coeff          * * cosine/squared 1.00 1.00 1.00 wca
pair_coeff          2 2 cosine/squared 1.00 ${rdist} ${rdist} wca
pair_coeff          3 3 cosine/squared 1.00 ${rdist} ${rdist} wca
pair_coeff          5 5 cosine/squared 1.00 ${rdist} ${rdist} wca
pair_coeff          6 6 cosine/squared 1.00 ${rdist} ${rdist} wca
pair_coeff          2 3 cosine/squared ${eps} 1.00 ${coff} wca 			# top interaction 2-3
pair_coeff          6 5 cosine/squared ${eps} 1.00 ${coff} wca 			# bottom interaction 6-5
pair_coeff          3 5 cosine/squared ${epsB} 1.00 ${coff} wca 			# B pair intra 3-5

compute             cPE all pe/atom pair
compute             cN2 all coord/atom cutoff ${bcoff} 2
compute             cN3 all coord/atom cutoff ${bcoff} 3
compute             cN5 all coord/atom cutoff ${bcoff} 5
compute             cN6 all coord/atom cutoff ${bcoff} 6
variable            vN2 atom c_cN2
variable            vN3 atom c_cN3
variable            vN5 atom c_cN5
variable            vN6 atom c_cN6
variable            vTI equal v_vN6[13]+v_vN2[11]
variable            vBI equal v_vN5[14]+v_vN3[10]

fix                 fLang moving langevin 1.0 1.0 1.0 ${seed}
fix                 fNVE moving nve

dump                1 all custom ${run_steps} output_0.xyz id mol type x y c_cPE v_vN2 v_vN3 v_vN5 v_vN6
dump_modify         1 format line "%d %d %d %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f" first yes

thermo              ${run_steps}
thermo_style        custom step temp pe ke etotal epair ebond press atoms v_vTI v_vBI

compute             cBonds all property/local batom1 batom2 btype
compute             cBondDxys all bond/local engpot force dist

dump                2 all local ${run_steps} bonds_1.dump c_cBonds[*]
dump_modify         2 format line "%.0f %.0f %.0f" first yes

fix                 twodim all enforce2d

fix                 fHaltT all halt 1 v_vTI == 1 error continue
fix                 fHaltB all halt 1 v_vBI == 1 error continue

timestep            ${tstep}
run                 ${run_steps}

dump                3 all custom ${run_steps} output_1.xyz id mol type x y c_cPE v_vN2 v_vN3 v_vN5 v_vN6
dump_modify         3 format line "%d %d %d %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f" first yes

dump                4 all local ${run_steps} bonds_2.dump c_cBonds[*]
dump_modify         4 format line "%.0f %.0f %.0f" first yes

timestep            ${tstep}
run                 ${run_steps}
''')
f.close()

f = open('%s/in1.local'%(gpath), 'w')
f.write('''# Input script for LAMMPS simulation of binding two chevron/triangle molecules

units               lj
atom_style          molecular
dimension           2 
boundary            p p p
log                 log1.txt
read_data           configuration.dat

group               moving id > 0

variable            Kbond equal %.1f                           	# bond constant [kT/sigma2]
variable            eps equal %.1f                           	# binding constant [kT]
variable            epsA equal %.1f                           	# binding constant A pair [kT]
variable            epsB equal %.1f                           	# binding constant B pair [kT]
variable            coff equal %f                           	# cutoff distance [sigma]
variable            bcoff equal %f                           	# cutoff distance [sigma]
variable            rdist equal 1.0                           	# cutoff distance [sigma]
variable            tstep equal %f                             	# simulation timestep size [seconds]
variable            seed equal %d                               	# random number generator seed
variable            run_steps equal %d          					# simulation run time [simulation steps]
variable            dump_steps equal %d        						# dumping interval [simulation steps]
'''%(Kbond,eps,epsA,epsB,coff,bcoff,tstep,seed,runsteps,dump))
f.write('''
special_bonds       lj 0.0 1.0 1.0
bond_style          harmonic
bond_coeff          1 ${Kbond} 1.0
bond_coeff          2 ${Kbond} 1.7320508075688774

pair_style          hybrid/overlay zero 2.00 cosine/squared 2.00
pair_coeff          * * cosine/squared 0.00 1.00 1.10
pair_coeff          * * zero 2.00
pair_coeff          * * cosine/squared 1.00 1.00 1.00 wca
pair_coeff          2 2 cosine/squared 1.00 ${rdist} ${rdist} wca
pair_coeff          3 3 cosine/squared 1.00 ${rdist} ${rdist} wca
pair_coeff          5 5 cosine/squared 1.00 ${rdist} ${rdist} wca
pair_coeff          6 6 cosine/squared 1.00 ${rdist} ${rdist} wca
pair_coeff          2 3 cosine/squared ${eps} 1.00 ${coff} wca 			# top interaction 2-3
pair_coeff          6 5 cosine/squared ${eps} 1.00 ${coff} wca 			# bottom interaction 6-5
pair_coeff          3 5 cosine/squared ${epsB} 1.00 ${coff} wca 			# B pair intra 3-5

compute             cPE all pe/atom pair
compute             cN2 all coord/atom cutoff ${bcoff} 2
compute             cN3 all coord/atom cutoff ${bcoff} 3
compute             cN5 all coord/atom cutoff ${bcoff} 5
compute             cN6 all coord/atom cutoff ${bcoff} 6
variable            vN2 atom c_cN2
variable            vN3 atom c_cN3
variable            vN5 atom c_cN5
variable            vN6 atom c_cN6
variable            vTI equal v_vN6[13]+v_vN2[11]
variable            vBI equal v_vN5[14]+v_vN3[10]

fix                 fLang moving langevin 1.0 1.0 1.0 ${seed}
fix                 fNVE moving nve

dump                1 all custom ${run_steps} output_0.xyz id mol type x y c_cPE v_vN2 v_vN3 v_vN5 v_vN6
dump_modify         1 format line "%d %d %d %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f" first yes

thermo              ${run_steps}
thermo_style        custom step temp pe ke etotal epair ebond press atoms v_vTI v_vBI

compute             cBonds all property/local batom1 batom2 btype
compute             cBondDxys all bond/local engpot force dist

dump                2 all local ${run_steps} bonds_0.dump c_cBonds[*]
dump_modify         2 format line "%.0f %.0f %.0f" first yes

fix                 twodim all enforce2d

fix                 fHaltT all halt 1 v_vTI == 2 error continue
fix                 fHaltB all halt 1 v_vBI == 2 error continue

timestep            ${tstep}
run                 ${run_steps}

dump                3 all custom ${run_steps} output_2.xyz id mol type x y c_cPE v_vN2 v_vN3 v_vN5 v_vN6
dump_modify         3 format line "%d %d %d %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f" first yes

dump                4 all local ${run_steps} bonds_2.dump c_cBonds[*]
dump_modify         4 format line "%.0f %.0f %.0f" first yes

timestep            ${tstep}
run                 ${run_steps}
''')
f.close()
