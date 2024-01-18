#!/usr/bin/env python

import os, sys, re, copy 
sys.path.append('../charles/scripts')
sys.path.append('../scripts')
import parsec

ALLOWED_SEARCHTERMS = parsec.OUTPUT_FIELDS
def parse_benchmarks(ex_basedir = {},**kwargs):
	''' check that key_ is allowed '''
	basedir = {}
	basedir['parsecs'] = ex_basedir.get('parsecs','./')

	key_ = kwargs.get('key_')

	if not key_ in parsec.OUTPUT_FIELDS: raise KeyError('%s is not stored currently\nAllowed Keys:\n' % (key_)+str(parsec.OUTPUT_FIELDS))
	
	parsec_version  = [l for l in os.listdir(basedir['parsecs']) if (os.path.isdir(os.path.join(basedir['parsecs'],l,'parsec','benchmarks')) and l[0] != '.')]
	print parsec_version
	compileddata = {}
	for l in parsec_version: compileddata[l] = {}
	for parsec_ in parsec_version:
		benchmarkdir = os.path.join(basedir['parsecs'],parsec_,'parsec','benchmarks')
		benchmarks = [l for l in os.listdir(benchmarkdir) if os.path.isdir(os.path.join(benchmarkdir,l))]
		for benchmark in benchmarks:
			'''parse parsec.out'''
			try:
				file = open(os.path.join(benchmarkdir,benchmark)+'/parsec.out','r')
				thisrun = parsec.Parsec_Interface(os.environ['WORK']+'/parsec.mpi')
				thisrun.parse_output(file=file,searchterm_=key_)
				compileddata[parsec_][benchmark] = eval('thisrun.%s' % key_)
			except IOError:
				''' parsec file wasn't there '''
			except IndexError: 
				'''Probably couldn't find the values required - probably means the case started but didn't finish'''
				compileddata[parsec_][benchmark] = -1
				file.close()
			
	return compileddata

if __name__ == '__main__':
	basedir = {}
	basedir['parsecs'] = sys.argv[sys.argv.index('--parsecs')+1] if ('--parsecs' in sys.argv) else os.getcwd()
	compileddata = parse_benchmarks(key_ = "Stats['SREs']", ex_basedir=basedir)

