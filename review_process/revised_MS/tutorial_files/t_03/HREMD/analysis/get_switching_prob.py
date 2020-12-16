#!/usr/bin/python

import sys, os

if not len(sys.argv) == 3:
  print "#\n#usage: ./get_switching_prob.py replica.dat 24\n#"
  sys.exit()

replica_file = sys.argv[1]
nr_replicas = int(sys.argv[2])


f = open(replica_file,"r")

switches = []
pairs = []
for i in range(1,nr_replicas,1):
  switches.append([])
  pairs.append([i,i+1])

read = False
for l in f:
  if read:
    line = l.split()
    rep1, rep2, switch = int(line[0]), int(line[1]), float(line[10])
    if [rep1,rep2] in pairs:
      index = pairs.index([rep1,rep2])
      switches[index].append(switch)
  elif l.startswith("#"):
    read = True
  else:
    continue
f.close()

for i in range(len(switches)):
  prob = sum(switches[i])/len(switches[i])
  print "pair: ", pairs[i][0],"-", pairs[i][1], "\t", prob

  
