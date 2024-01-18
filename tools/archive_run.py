#!/usr/bin/env python
import tarfile, sys, os

overwrite = True if sys.argv.count('--force') else False
tarname = sys.argv[sys.argv.index('--name')+1] if sys.argv.count('--name') else 'default.tar.gz'
parsecdir = sys.argv[sys.argv.index('--parsecs')+1] if sys.argv.count('--parsecs') else os.getcwd()

if os.path.exists(tarname) and not overwrite:
	msg = '''%s already exists.\nIf you desire to overwrite pre-existing archive, use --force''' % (tarname)
	raise ValueError(msg)
try:
	archive = tarfile.open(tarname,'w:gz')
except tarfile.CompressionError:
	print '''Can't use zlib through python in this configuration - disabling compression'''
	archive = tarfile.open(tarname,'w:')



parsecs = [l for l in os.listdir(parsecdir) if (os.path.isdir(os.path.join(parsecdir,l,'parsec','benchmarks')) and l[0] != '.')]
if not len(parsecs):
	msg = '''No parsecs found at %s\nUse --parsecs <path> to give the base directory containing the parsec subdirectories''' % (parsecdir)
	raise ValueError(msg)

for parsec in parsecs:
	benchdir = os.path.join(parsecdir,parsec,'parsec','benchmarks')
	benchmarks = [l for l in os.listdir(benchdir) if (os.path.isdir(os.path.join(benchdir,l)) and l[0] != '.')]
	for benchmark in benchmarks:
		archive.add(os.path.join(benchdir,benchmark,'parsec.out'))
		outfiles = [file for file in os.listdir(os.path.join(benchdir,benchmark)) if (file[0:4] == 'out.')]
		for outfile in outfiles: 
			print os.path.join(benchdir,benchmark,outfile)
			archive.add(os.path.join(benchdir,benchmark,outfile))
		print(os.path.join(benchdir,benchmark,'parsec.out'))
		archive.add(os.path.join(benchdir,benchmark,'parsec.out'))
		print(os.path.join(benchdir,benchmark,'benchmark.out'))
		archive.add(os.path.join(benchdir,benchmark,'benchmark.out'))

if os.path.exists('Records.txt'): 
	print('Records.txt')
	archive.add('Records.txt')

archive.close()
