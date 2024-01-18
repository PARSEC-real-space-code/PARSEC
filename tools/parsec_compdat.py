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

# print out all the vector based ones
for searchterm in searchterms:
	compileddata = parse_benchmarks(key_ = searchterm, ex_basedir=basedir)

	theParsecKeys    = compileddata.keys()
	if type(compileddata[compileddata.keys()[0]].values()) != type(list()): continue
	theBenchmarkKeys = set([])
	for l in compileddata:
		theBenchmarkKeys = theBenchmarkKeys.union(set(tuple(compileddata[l].keys()))) 
	theBenchmarkKeys = list(theBenchmarkKeys)
	# counts of passes (per version)
	passes = [len(compileddata[l].values())-compileddata[l].values().count(-1.00) for l in compileddata]
	theParsecKeys = sorted( theParsecKeys, key = lambda t: passes[theParsecKeys.index(t)] ) 
	# counts number of failures (per benchmark)
	failures = [ [compileddata[key].get(l,-1.00) for key in theParsecKeys].count(-1.00) for l in theBenchmarkKeys ]
	theBenchmarkKeys = sorted( theBenchmarkKeys, key = lambda t: failures[theBenchmarkKeys.index(t)])
	if isSortedByTitle:
		theBenchmarkKeys = sorted( theBenchmarkKeys, key = lambda t: t)
	BenchmarkKeyWidth = [max(15,len(k)+1) for k in theBenchmarkKeys]
	ParsecKeyWidth = max([len(k) for k in theParsecKeys])
	''' find MAKEFILE flags '''
	F = open('Records.txt','a')
	F.write(searchterm)
	for start in  range(0,len(theBenchmarkKeys),NUMPRINT_BENCHMARKS):
		stop = start + min(NUMPRINT_BENCHMARKS,len(theBenchmarkKeys)-start)
		fmt_data = []
		for l in zip(BenchmarkKeyWidth[start:stop],theBenchmarkKeys[start:stop]):
			fmt_data.append(l[0])
			fmt_data.append(l[1])

		fmt = '\n'+('%-'+'%d' % (ParsecKeyWidth+1)+'s ') % 'Benchmarks' + '%*s'*(stop-start)
		print fmt % (tuple(fmt_data))
		F.write(fmt % (tuple(fmt_data))+'\n')
		for l in theParsecKeys:
			dat = [l]
			fmt = '%-'+'%d' % (ParsecKeyWidth+1)+'s '+'%*.6f'*(stop-start)

			# next line is overly complicated but guarantees the same ordering in each one
			i = 0
			for k in theBenchmarkKeys[start:stop]: 
				dat.extend( (fmt_data[2*i],compileddata[l].get(k,numpy.Inf)) )
				i+=1

			print  fmt % (tuple(dat))
			F.write(fmt % (tuple(dat))+'\n')

	F.write('\n')
	F.close()

