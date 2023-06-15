import sys
import glob
import os

wd = os.getcwd()
command = '%s/LAMMPS/lammps/src/lmp_mpi -in in.local'%(wd.split('/chevrons/')[0])

gpath = input("What is the general path to this simulation set? (please include full tree!) ")
print("OK! Will move to %s"%(gpath))

r = os.chdir(gpath)
print("Moved to %s"%(gpath))

if os.access('submission_files/job_array.sh', os.F_OK):
	print("Submission files here already exist! Either remove them if you want new ones or make sure you are aiming for the right directory...")
	print("To remove these submission files please run 'rm -r %s/submission_files'"%(gpath))
	exit()

r = os.system("mkdir submission_files")

paths_seed = input("Please provide a seed for glob.glob() to find all simulations to run. Only give the local tree (from here) ")
print("OK! Will sample %s"%(paths_seed))

paths = glob.glob(paths_seed)
check = input("%d simulations to run here... Are you sure? [y/n] "%(len(paths)))

wall_hours = int(input("What wall time hours would you like? (in hours please, only integer numbers) "))
wall_mins = int(input("What wall time minutes would you like? (in minutes please, only integer numbers) "))
print("OK! Will set a wall time of %d hours and %d minutes!"%(wall_hours, wall_mins))

memory = int(input("How much memory would you like? (in GB please, only integer numbers) "))
print("OK! Will request %d GB!"%(memory))

gname = input("What general name would you like to give this set of simulations? ")
print("OK! Naming these simulations %s"%(gname))

if check == 'n':
	print("Submission files writing cancelled. Aborting...")
	exit()
elif check == 'y':
	print("Proceeding...")
else:
	print("Wrong answer! Please try again :(")
	print("Aborting...")
	exit()

r = os.chdir('submission_files')

spath = os.getcwd()

f = open("paths.txt", "w")
g = open("names.txt", "w")
for p in paths:
	f.write("%s/%s\n"%(gpath,p))
	ps = p.split('/')
	n = ps[0]
	for i in range(len(ps)-1):
		n = '%s_%s'%(n,ps[i+1])
	g.write("%s\n"%(n))
f.close()
g.close()

f = open("job_array.sh", "w")
f.write('''#!/bin/bash
#
''')
f.write("#SBATCH --job-name=%s\n"%(gname))
f.write("#SBATCH --output=%s"%(gname))
f.write("_%j.out\n")
f.write("#SBATCH --array=1-%d:1\n"%(len(paths)))
f.write('''#
#Maximum runtime is limited to 10 days, ie. 240 hours
''')
f.write("#SBATCH --time=%d:%d:00\n#\n#SBATCH --mem=%dG\n"%(wall_hours,wall_mins,memory))
f.write('''#
#Do not requeue the job in the case it fails.
#SBATCH --no-requeue
#
#Do not export the local environment to the compute nodes
#SBATCH --export=NONE
unset SLURM_EXPORT_ENV
#
#for single-CPU jobs make sure that they use a single thread
export OMP_NUM_THREADS=1
#
#run the respective binary through SLURM's srun
path=$(sed "${SLURM_ARRAY_TASK_ID}q;d" paths.txt)
job_name=$(sed "${SLURM_ARRAY_TASK_ID}q;d" names.txt)
cd ${path}
echo init
echo task job id:       ${SLURM_JOB_ID}
echo array job id:      ${SLURM_ARRAY_JOB_ID}
echo array task index:  ${SLURM_ARRAY_TASK_ID}
echo simulation folder: ${path}
echo name:              ${job_name}
echo hostname:          ${SLURM_JOB_NODELIST}
echo submit host:       ${SLURM_SUBMIT_HOST}
''')
f.write("ln %s/%s_${SLURM_JOB_ID}.out %s/${job_name}.out\n"%(spath,gname,spath))
f.write("srun --cpu_bind=verbose %s\n"%(command))
f.close()

