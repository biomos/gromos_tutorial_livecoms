#!/usr/bin/env python

import os
import sys

workdir = os.getcwd()

usage0 = '*******************\n'
usage1 = 'Rep_ana program for mpi release (becuase its not working)\n'
usage2 = 'usage rep_change [replica info] @type [int]\n'
usage3 = '@type flag define replica info format :\n 0 = new format (def)\n 1 = old format\n'
usage = usage0+usage1+usage2+usage3+usage0

#reading in command line

argnum = (len(sys.argv)-1)

if argnum == 0:
	print 'no filename was added, program stops'
	print usage
	sys.exit(0)

fname = sys.argv[1]
path = workdir+'/'+fname
if os.path.isfile(path) == True:
	print('File '+path+' found')
else:
	print 'no file was found under '+path
	print usage
	sys.exit(0)

#print('flags detected:')
flag = ''
tflag = 0
cnt = 1
while cnt <= argnum:
	arg = sys.argv[cnt]
	if arg.startswith('@'):
		flag = arg.strip('@')
#		print flag
	elif flag == 'type':
		tflag = arg
	cnt = cnt + 1

#reading file

#print 'reading in info file'
i = open(fname,'rb')
mark = 0
Tnum = 0
TEMP = []
Lnum = 0
LAMBDA = []
STEPS = []

for line in i:
	line2 = line.strip().strip('\n')
	if line.startswith('#'):
		mark = 1	
	elif line.startswith('Number of temperatures:'):
		flag, Tnum = line2.split(":")
		Tnum = int(Tnum)
	#	print('number of temp: '+str(Tnum))		
	elif line.startswith('Number of lambda values:'):
		flag, Lnum = line2.split(":")
		Lnum = int(Lnum)
	elif line.startswith('T'):
		line3 = line2.split()
		tcnt = 1
		while tcnt <= Tnum:
			tempi = line3[tcnt]
			TEMP.append(tempi)
			tcnt = tcnt + 1
	elif line.startswith('lambda'):
		line3 = line2.split()
		lcnt = 1
		while lcnt <= Lnum:
			lamdi = line3[lcnt]
			LAMBDA.append(lamdi)
			lcnt = lcnt + 1
		flag = 'go'
	elif mark == 1 and flag == 'go':
		if tflag == 0:
			try:
				ID, part, run, L1, T1, E1, L2, T2, E2, prob, s = line2.split()
				ID = int(ID)
				part = int(part)
				run = int(run)
				s = int(s)
				step = [ID, part, run, L1, T1, E1, L2, T2, E2, prob, s]
				STEPS.append(step)
			#	print str(step)
			except ValueError:
				print 'could not read in line:\n'+line
		else:
			try:
				ID, run, T1, L1, E1, T2, L2, E2, prob, s = line2.split()
				print 'not yet implemented'
				print usage
				sys.exit(0)
				ID = int(ID)
				run = int(run)
				s = int(s)
				step = [ID, run, T1, L1, E1, T2, L2, E2, prob, s]
				STEPS.append(step)
				print str(step)
			except ValueError:

				print 'could not read in line:\n'+line
# sorting out data
runnum = run
#print ' number of last trial:'+str(runnum)

#if Lnum == len(LAMBDA):
#	print 'number of Lambda values found: '+str(Lnum)
#	print str(LAMBDA)
#else:
#	print 'number of Lambda values found: '+Lnum
#	print 'but only following values are found:\n'+str(LAMBDA)
#	print 'program stops'
#	sys.exit(0)


LISTS = []
REPLICA = []
ORDER = []
cnt = 0
for Lambda in LAMBDA:
#	info = '#replica exchange sceme for replica '+str(cnt)+'\n'
	List = [('0	'+Lambda+'\n')]
	LISTS.append(List)
#	print str(List)
	cnt = cnt + 1
	REPLICA.append(cnt)
	ORDER.append("0")

BLOCK = []
if tflag == 0:
	bcnt = STEPS[0][2]
elif tflag == 1:
	bcnt = STEPS[0][1]
scnt = 0
for step in STEPS:
	if tflag == 1:
		ID =  step[0]
		run = step[1]
		L1 = step[3]
		L2 = step[6]
		lock = step[9]
# grouping data to follow replica exchange
	else:
		run = step[2]
	if run == bcnt:
		BLOCK.append(step)
		scnt = scnt + 1
	if scnt == Lnum:	
	#	print 'Block '+str(bcnt)
		bcnt = bcnt + 1
		rcnt = 0
# making replicawise order
		while rcnt < Lnum:	
			rep = REPLICA[rcnt]
			for brick in BLOCK:
				ID = brick[0]
				part = brick[1]
				run = brick[2]
				L1 = brick[3]
				L2 = brick[6]
				lock = brick[10]
				if rep == ID:
					if lock == 1:
						ORDER[rcnt] = part
						newline = str(run)+'	'+str(L2)+'\n'
					#	print (brick),
					#	print part
						LISTS[rcnt].append(newline)
					else:
						ORDER[rcnt] = rep
						newline = str(run)+'    '+str(L1)+'\n'
					#	print (brick),
					#	print ID
						LISTS[rcnt].append(newline)
			rcnt = rcnt + 1
	#	print str(ORDER)
 		REPLICA = ORDER
		BLOCK = []
		scnt = 0	

#writing out data:
cnt = 0
o = open('rep_change.out','wb')
while cnt < Lnum:
	o.write('#replica exchange sceme for replica '+str(cnt+1)+'\n')
	for line in LISTS[cnt]:
		o.write(line)
	o.write('&\n')
	cnt = cnt + 1
