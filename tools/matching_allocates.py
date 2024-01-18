#!/usr/bin/env python
import re, os, copy

allocates = {}
deallocates = {}

dir = 'parsecs/1.4/parsec/src'

files = os.listdir(dir)
files = [os.path.join(dir,file) for file in files if not os.path.isdir(os.path.join(dir,file))]
for file in files:
	#find allocates
	F = open(file,'r')
	data = F.read()
	results = re.findall(r'[^e]allocate\(([\S ,]*)\)',data,re.IGNORECASE)
	allocates[os.path.basename(file)] = copy.deepcopy(results)
	results = re.findall(r'deallocate\((\S+)',data,re.IGNORECASE)
	deallocates[os.path.basename(file)] = copy.deepcopy(results)

	F.close()

