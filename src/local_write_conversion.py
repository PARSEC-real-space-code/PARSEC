#!/usr/bin/env python
'''This changes all source files that have write(unit#,
   to either a new number, or comments them out entirely.
   What will probably need to happen at some future date are
   comment guards and a debug version of parsec / or verbose.
   We don't just want if checks because those require CPU thought.
   I'm going to keep this around for automatically changing your source files,
   Since this knowledge will be rather slow in trickling around.
   -Charles'''

import re, os, sys
searchunit = 6
replaceunit= 9
USE_IFDEFGUARDS = False
EXECUTE = False if not '-confirm' in sys.argv else True
OVERRIDE = False if not '-override' in sys.argv else True
files = [file for file in os.listdir('./') if (os.path.splitext(file)[1] in ('.f90','.F90','.f90p','.f90z','.f'))]
if USE_IFDEFGUARDS: verbstr = ('\n#if defined(VERBOSE)','\n#endif')
else:               verbstr = ('','')

count = []
for file in files:
	F = open(file,'r')
	data = F.read()
	F.close()
	
	# finds the open({searchunit},) and replaces it as well
	findstr = r'\n(\s*\S*\s*)open(\s*)\((\s*[unit]*\s*[=]*\s*)%d(\s*)([ ,\)])(.*)' % (replaceunit)
	numOccurrence = len(re.findall(findstr,data))
	if numOccurrence > 0: 
		count.append((file,numOccurrence))

if len(count) > 0 and not OVERRIDE:
	for c in count: print c
	e = Exception('replaceunit (%d) already exists in this fortran project!\nChoose a different unit, or change that existing one first!' % (replaceunit))
	raise e

for file in files:
	F = open(file,'r')
	data = F.read()
	F.close()
	lines = data.split('\n')
	i = 1
	filen = 0
	for line in lines:
		# find write(unit=num or write(num 
		for str in ['write','open','myflush','close']:
			findstr = r'(\s*\S*\s*)%s(\s*)\((\s*[unit]*\s*[=]*\s*)%d(\s*)([ ,\)])(.*)' % (str,searchunit)
			repstr = r'%s\1%s\2(\g<3>%d\4\5\6%s' % (verbstr[0],str,replaceunit,verbstr[1])
			(line,n) = re.subn(findstr, repstr, line)
			if n > 0: 
				lines[i-1] = line
				filen += n
				print i,lines[i-1]
		i+=1

	# finds the open({searchunit},) and replaces it as well
	#findstr = r'\n(\s*\S*\s*)open(\s*)\((\s*[unit]*\s*[=]*\s*)%d(\s*)([ ,\)])(.*)' % (searchunit)
	#repstr = r'%s\n\1open\2(%d\4\5\6%s' % (verbstr[0],replaceunit,verbstr[1])
	#data = re.sub(findstr, repstr, data)

	if (EXECUTE and filen > 0):
		print 'Writing to file.........%s' % file
		F = open(file,'w')
		F.write('%s\n'*(len(lines)-1) % (tuple(lines[0:-1])))
		F.write('%s' % lines[-1])
		F.close()

	
