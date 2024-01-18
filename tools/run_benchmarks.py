#!/usr/bin/env python
import os,sys,re

accepted_options = ['--parsecs','--split','--timelimit','--exec','--title','--preamble','--ext']
basedir = {}
options = {}
if ('--help' in sys.argv) or ('-h' in sys.argv) or ('help' in sys.argv): 
	print('''
	Usage: run_benchmarks.py [OPTIONS]
	Desc:  Makes separate job scripts for running shorter benchmarks sequentially on a single node
		--parsecs  designates what folder to use as the directory with all the parsec branches in it : default is os.getcwd()
		--split=[1] how many files to create
		--timelimit how long to restrain the individual commands to - default is 15 minutes "15m"
		--exec      the verb used to execute parsec program - default ibrun
		--title     titles the jobs
		--preamble  includes the submission queue for my use
		--ext       default is mpi, but you might be using ser if you use valgrind
	
	Example using Valgrind:
		./run_benchmarks.py --parsecs comborun/ \
		--split 5 --timelimit 2h \
		--exec "valgrind --track-origins=yes --leak-check=full --read-var-info=yes" \
		--preamble lonestar --ext ser

''')

for k in [l for l in sys.argv if (l[0:2] == '--') ] :
	if k not in accepted_options:
		print('Error: %s not recognized option. Use help to see currently allowed options' % (k))
		sys.exit(1)

basedir['parsecs'] = sys.argv[sys.argv.index('--parsecs')+1] if ('--parsecs' in sys.argv) else os.getcwd()
options['split'] = int(sys.argv[sys.argv.index('--split')+1]) if ('--split' in sys.argv) else 0
options['timelimit'] = ('timeout %s' % sys.argv[sys.argv.index('--timelimit')+1]) if ('--timelimit' in sys.argv) else ''
options['exec'] = (sys.argv[sys.argv.index('--exec')+1]) if ('--exec' in sys.argv) else 'ibrun' 
options['title'] = (sys.argv[sys.argv.index('--title')+1]) if ('--title' in sys.argv) else 'untitled' 
options['preamble'] = (sys.argv[sys.argv.index('--preamble')+1]) if ('--preamble' in sys.argv) else 'lonestar'
options['ext'] = (sys.argv[sys.argv.index('--ext')+1]) if ('--ext' in sys.argv) else 'mpi'

parsec_version  = [l for l in os.listdir(basedir['parsecs']) if (os.path.isdir(os.path.join(basedir['parsecs'],l,'parsec','benchmarks')) and l[0] != '.')]
print parsec_version
jobs = []
for dir in parsec_version: 
	subdir = os.path.join(basedir['parsecs'],dir,'parsec','benchmarks')
	print os.listdir(subdir)
	benchmarks = [l for l in os.listdir(subdir) if (l[0] != '.' and os.path.isdir(os.path.join(subdir,l)))]
	for benchmark in benchmarks:
		jobs.append('cd '+os.path.join(os.getcwd(),subdir,benchmark)+' ; %s %s ../../src/parsec.%s &> benchmark.out' % (options['timelimit'],options['exec'],options['ext']))
options['split'] = len(jobs)+1 if options['split']==0 else options['split']
for i in range(0,(len(jobs)+options['split'])/options['split']):
	start = i*options['split']
	end   = min((i+1)*options['split'],len(jobs))
	F = open('jobscript_%05d.sh' % (i),'w')
	if options['preamble']=='stampede':
		F.write('''#!/bin/bash
#SBATCH -J %s-bench   # Job Name
#SBATCH -o %s.o%%j    # Output and error file name (%%j expands to jobID)
#SBATCH -n 16           # Total number of mpi tasks requested
#SBATCH -p development  # Queue (partition) name -- normal, development, etc.
#SBATCH -t 02:00:00
#SBATCH -A TG-DMR090026\n''' % (options['title'],options['title']))
	if options['preamble']=='lonestar':
		F.write('''#!/bin/bash
#$ -N %s-bench   # Job Name
#$ -V            #
#$ -cwd
#$ -o %s.o$JOB_ID    # Output and error file name (%%j expands to jobID)
#$ -pe 12way 12           # Total number of mpi tasks requested
#$ -q normal  # Queue (partition) name -- normal, development, etc.
#$ -l h_rt=10:00:00
#$ -A Nanomaterials\n''' % (options['title'],options['title']))
	F.write(('%s\n'*len(jobs[start:end])) % (tuple(jobs[start:end])))
	F.close()

print jobs[0:2],jobs[-3:]
