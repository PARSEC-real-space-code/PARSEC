#!/usr/bin/env python
import os, sys, re, shutil as sh
import subprocess as sp

accepted_options = ['--parsecs','--makefile','--cleanall','--clean-only']
basedir = {}
options = {}

if ('--help' in sys.argv) or ('-h' in sys.argv) or ('help' in sys.argv): 
	print('''
	Usage: make_parsecs [OPTIONS]
	Desc:  Compiles all the parsec versions either as is or using a new makefile include.
		--parsecs  designates what folder to use as the directory with all the parsec branches in it : default is os.getcwd()
		--makefile designates what submake file to use (examples located in the parsec+/parsec/src/config directory) - otherwise keeps the include to whatever it was originally.
		--cleanall forces a make cleanall operation before a make (if makefile provided)
		--clean-only forces a make cleanall for each parsec and then exits
''')

	sys.exit(0)

for k in [l for l in sys.argv if (l[0:2] == '--') ] :
	if k not in accepted_options:
		print('Error: %s not recognized option. Use help to see currently allowed options' % (k))
		sys.exit(1)

basedir['parsecs'] = sys.argv[sys.argv.index('--parsecs')+1] if ('--parsecs' in sys.argv) else os.getcwd()
basedir['makefile'] = sys.argv[sys.argv.index('--makefile')+1] if ('--makefile' in sys.argv) else ''

CLEAN_BEFORE_MAKE = True
CLEAN_ONLY = False
print basedir

returncode = {}

''' get all the parsec versions '''
parsecs = [os.path.join(basedir['parsecs'],k,'parsec','src') for k in os.listdir(basedir['parsecs']) if (os.path.isdir(os.path.join(basedir['parsecs'],k)) and k[0] != '.' and os.path.isdir(os.path.join(basedir['parsecs'],k+'/parsec')))]

if not len(parsecs):
	print "Zero parsec directories found in %s\n\tlooking for parsecs directory for quick bootstrap\n" % (basedir['parsecs'])
	parsecs = [os.path.join('parsecs',k,'parsec','src') for k in os.listdir('parsecs') if (os.path.isdir(os.path.join('parsecs',k)) and k[0] != '.' and os.path.isdir(os.path.join('parsecs',k+'/parsec')))]
	if not len(parsecs):
		print "Zero parsec directories found in 'parsecs': exiting now..."
		sys.exit(1)

''' check to see if makefile is a file '''
NEW_MAKEFILE = os.path.isfile(basedir['makefile'])
make_src = basedir['makefile']

for parsec in parsecs:
	print(parsec+'....')
	''' copy makefile into various parsecs if it exists '''
	make_dst = os.path.join(parsec,'config',os.path.basename(basedir['makefile']))
	print make_dst
	if NEW_MAKEFILE:
		sh.copy2(make_src,make_dst)
		# read in old Makefile
		F = open(os.path.join(parsec,'Makefile'),'r')
		make_data = F.read()
		F.close()
		# write new makefile
		F = open(os.path.join(parsec,'Makefile'),'w')
		(new_data,n)  = re.subn(r"include\s+config/\S+","include config/%s" %  (os.path.basename(basedir['makefile'])),make_data)
		F.write(new_data)
		F.close()
		if n != 1: raise ValueError("n = %d\n%500s\n" % (n,new_data))

	''' invoke make '''
	taskmake_str = 'make -C %s cleanall ; ' % (parsec) if (('--cleanall' in sys.argv) or ('--clean-only' in sys.argv)) else ''
	taskmake_str += ('make -j8 -C %s' % (parsec) if not ('--cleanonly' in sys.argv) else '')
	print taskmake_str

	try: 
		taskmake = sp.Popen(taskmake_str,shell=True,stdout=sp.PIPE,stdin=sp.PIPE)
		taskmake.communicate()[0]
		returncode[ parsec.split('/')[-3] ] = taskmake.returncode
	except: 
		print taskmake_str
		raise


print returncode

