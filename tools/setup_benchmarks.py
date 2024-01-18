#!/usr/bin/env python
import os, sys, re, shutil

basedir = {}
options = {}
print('####### Starting setup_benchmarks.py #######')
if (len(sys.argv) < 2):
	print('''
	Usage: setup_benchmarks <path/to/folder/with/benchmark/folders> [OPTIONS]
	Desc:  Setup benchmarks and job files for each parsec version detected in the current directory. Job files are done on a per-version basis to reduce the total number of jobs. This allows you to go to a $PARSEC/parsec/benchmarks folder and use one command to submit all the jobs
		--parsec  designates what folder to use as the directory with all the parsec branches in it : default is os.getcwd()
		--jobs    designates what job system to create jobs for - default is SLURM (stampede)
		--clean   removes all benchmarks prior to starting the copy
''')

	sys.exit(1)


basedir['benchmarks'] = sys.argv[1]
basedir['parsecs'] = sys.argv[sys.argv.index('--parsecs')+1] if ('--parsecs' in sys.argv) else os.getcwd()
options['jobs'] = sys.argv[sys.argv.index('--jobs')+1] if ('--jobs' in sys.argv) else 'BASH'
options['clean_directory'] = ('--clean' in sys.argv)

print basedir

''' get all benchmark directories '''
benchmarks = [k for k in os.listdir(basedir['benchmarks']) if (os.path.isdir(os.path.join(basedir['benchmarks'],k)) and k[0] != '.')]

''' get all the parsec versions '''
print('ParsecDir , Status')
for l in os.listdir(basedir['parsecs']):
	print(l, os.path.isdir(os.path.join(basedir['parsecs'],l,'parsec','benchmarks')))

parsecs = [os.path.join(basedir['parsecs'],l) for l in os.listdir(basedir['parsecs']) if (os.path.isdir(os.path.join(basedir['parsecs'],l,'parsec','benchmarks')) and l[0] != '.')]

if len(parsecs) == 0:
	e = Exception('Error: No parsecs found!')
	raise e

''' put all the benchmarks into the parsecs, cleaning out the parsec.out files ''' 
for parsec in parsecs:
	print(parsec+'....')
	if options['clean_directory'] : 
		print('...removing prior benchmarks')
		shutil.rmtree(os.path.join(parsec,'parsec/benchmarks'))
		os.mkdir(os.path.join(parsec,'parsec/benchmarks'))

	versionjobs = open(parsec+'/parsec/benchmarks/versionjobs.sbatch','w')
	for benchmark in benchmarks:
		try: shutil.rmtree(os.path.join(parsec,'parsec/benchmarks/'+benchmark))
		except: pass
		print('...copying benchmark %s' % (benchmark))
		shutil.copytree(os.path.join(basedir['benchmarks'],benchmark),os.path.join(parsec,'parsec/benchmarks/'+benchmark),ignore=shutil.ignore_patterns('.svn','.*','parsec.out'))
		F = open(parsec+'/parsec/benchmarks/'+benchmark+'/job.sbatch','w')
		if(options['jobs'] == 'SLURM'):
			outstr = '''#!/bin/bash
#SBATCH -J benchmark-%s   # Job Name
#SBATCH -o %s.o%%j    # Output and error file name (%%j expands to jobID)
#SBATCH -n 16           # Total number of mpi tasks requested
#SBATCH -p development  # Queue (partition) name -- normal, development, etc.
#SBATCH -t 00:30:00
#SBATCH -A TG-DMR090026
 ibrun ../../src/parsec.mpi
''' % (benchmark,benchmark)
		elif (options['jobs'] == 'BASH'): 
			outstr = '''#!/bin/bash
timeout 2h mpiexec -np 2 ../../src/parsec.mpi'''
		elif (options['jobs'] == 'SGE-LONESTAR'): 
			outstr = '''#!/bin/bash
#$ -V                           # Inherit the submission environment
#$ -cwd                         # Start job in submission directory
#$ -A Nanomaterials
#$ -N bench-%s                    # job Name
#$ -o bench-%s.o$JOB_ID               # Name of the output file (eg. myMPI.oJobID)
#$ -j y                         # combine stderr &
#$ -pe  12way 12                # Requests 12cores/node
#$ -q   normal                  # Queue name
#$ -l h_rt=00:45:00             # Run time (hh:mm:ss)
set -x
ibrun tacc_affinity ../../src/parsec.mpi
''' % (benchmark,benchmark)

		F.write(outstr)
		F.close()
		versionjobs.write('cd %s ; qsub job.sbatch ; cd ..\n' % (benchmark))
	versionjobs.close()


print basedir
print benchmarks
print parsecs

