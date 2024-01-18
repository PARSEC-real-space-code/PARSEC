#!/usr/bin/env python
from parse_benchmarks import *
import numpy
import sys
searchterms = sys.argv[sys.argv.index('-s')+1] if ('-s' in sys.argv) else ALLOWED_SEARCHTERMS[0:ALLOWED_SEARCHTERMS.index("Stats['SREs']")]
searchterms = sys.argv[sys.argv.index('--search')+1] if ('--search' in sys.argv) else searchterms
searchterms = searchterms.split(',') if type(searchterms)==type('') else searchterms
NUMPRINT_BENCHMARKS = int(sys.argv[sys.argv.index('-c')+1]) if ('-c' in sys.argv) else 4
NUMPRINT_BENCHMARKS = int(sys.argv[sys.argv.index('--columns')+1]) if ('--columns' in sys.argv) else NUMPRINT_BENCHMARKS
isSortedByTitle = '--sortbytitle' in sys.argv
basedir = {}
basedir['parsecs'] = sys.argv[sys.argv.index('--parsecs')+1] if ('--parsecs' in sys.argv) else os.getcwd()

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
# turn off updating

compileddata = parse_benchmarks(key_ = "Stats['SREs']", ex_basedir=basedir)

theParsecKeys    = compileddata.keys() # ex. '1.4'
theBenchmarkKeys = set([]) # ex. '10_6.00_nosym'
for l in compileddata:
	theBenchmarkKeys = theBenchmarkKeys.union(set(tuple(compileddata[l].keys()))) 
theBenchmarkKeys = list(theBenchmarkKeys)


for experiment in theBenchmarkKeys:
	curfig = plt.figure()
	curlegend = []
	for aparsec in theParsecKeys:
		try:
			if type(compileddata[aparsec][experiment]) == type(dict()):
				for scf_iter,datum in compileddata[aparsec][experiment].items():
					x = datum[0]
					y = datum[1]
					plt.semilogy(x,y,'-o')
					curlegend.append('%s - %d' % (aparsec,scf_iter))
			else:
				print "Skipping: %s, %s..." % (aparsec, experiment)
		except:
			pass

	plt.legend(curlegend)
	plt.title('Experiment %s' % (experiment))
	plt.xlabel('Iteration')
	plt.ylabel('SRE')
	plt.savefig('%s.eps' % (experiment))



