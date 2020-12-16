#!/usr/bin/python

import numpy as np
import sys, os
import argparse

parser = argparse.ArgumentParser(fromfile_prefix_chars='@')
parser.add_argument("topo", type=str,help="topology file (*top)")
parser.add_argument("df_file", type=str,help="file containing DF distance as a function of time")
parser.add_argument("traj", type=str,help="coordinate trajectory (*trc.gz)")
parser.add_argument("imd", type=str,help="input file of HREMD, containing the lambda values of the replicas")
parser.add_argument("disres", type=str,help="distance restraints file, containing the PERTDFRESSPEC block")
parser.add_argument("--gromos", type=str,help="directory where the gromos++ program frameout is located, if not specified it is assumed to be set in the path of the user",default="")
args = parser.parse_args()

#checking input:
if not args.topo.endswith("top"):
    print "Error: {0} should end with .top".format(args.topo)
    sys.exit()
if not (args.traj.endswith("trc") or args.traj.endswith("trc.gz")):
    print "Error: {0} should end with .trc or .trc.gz".format(args.traj)
    sys.exit()
if args.gromos:
  if not args.gromos.endswith("/"):
    args.gromos = args.gromos+'/'

#get value closest to reference value
def find_nearest(array,value):
    idx=(np.abs(array-value)).argmin()
    return idx+1, array[idx]

#return filename like FRAME_00001.cnf
def get_fileName(frNr):
  l = len(str(frNr))
  if l == 1:
    name = "FRAME_0000{0}.cnf".format(frNr)
  elif l == 2:
    name = "FRAME_000{0}.cnf".format(frNr)
  elif l == 3:
    name = "FRAME_00{0}.cnf".format(frNr)
  elif l == 4:
    name = "FRAME_0{0}.cnf".format(frNr)
  elif l == 5:
    name = "FRAME_{0}.cnf".format(frNr)
  elif l == 0:
    print "ERROR: no framenumber specified"
    sys.exit()
  return name

#############MAIN##############################################


#read DF distances from the file
df_dist = []
df_file = open(args.df_file,"r")
for line in df_file:
  if not line.startswith("#"):
    values = line.split() 
    df_dist.append(float(values[1]))

df_file.close()

#read lambda values from REPLICA block in input file of HREMD:
imd_file = open(args.imd,'r')
lambda_lines = False
lambda_values = []
for line in imd_file:
  if "RETS(1" in line:
    #all lambda values should be read, we can stop now
    break
  elif lambda_lines and not ("RETS(1" in line):
    #(continue) reading lambda values
    l = line.strip().split()
    for i in l:
      lambda_values.append(float(i))
  elif "RELAM(1" in line:
    #found the right spot
    lambda_lines = True

imd_file.close()
print "#lambda values:", lambda_values

#read distance restraint specification file
disres_file = open(args.disres,'r')
PERTDFRESSPEC_block = False
R0_line = False
A_R0, B_R0 = 0.0, 0.0
for line in disres_file:
  if "PERTDFRESSPEC" in line:
    #found right block
    PERDFRESSPEC_block = True
  elif "A_R0" in line and PERDFRESSPEC_block:
    #found right line
    R0_line = True
  elif R0_line:
    l = line.strip().split()
    A_R0 = float(l[1])
    B_R0 = float(l[3])
    break

disres_file.close()
print "#A_R0 = {0}, B_R0 = {1}".format(A_R0,B_R0)

#calculate restraining distances for each lambda
HREMD_dist = []
for L in lambda_values:
  x = (1-L)* A_R0 + L*B_R0
  HREMD_dist.append(x)

HREMD_dist= np.asarray(HREMD_dist)
print "#HREMD_dist:", HREMD_dist

#for each of the replicas, find the frame with the DF distance matching the restraining distance the closest
for i in HREMD_dist:
  frNr, dist = find_nearest(df_dist,i)
  print i, "frame: ", frNr, " actual distance: ", dist

  os.system("{0}frameout \
	@topo 		{1} \
	@pbc 		r cog \
	@outformat 	cnf \
	@include	ALL \
	@atomsfit	1:a \
	@traj		{2} \
	@spec		SPEC \
	@frames		{3}".format(args.gromos,args.topo,args.traj,frNr))



  #get filename
  name = get_fileName(frNr)
  #index of distance
  repNr = np.nonzero(HREMD_dist == i)[0][0] + 1
  #generate new name
  outfile = "start_{0}.cnf".format(repNr)
  #mv output of frameout (FRAME_*.cnf) to start_*.cnf)
  os.rename(name, outfile)
  print 'written frame {0} to {1}'.format(frNr,outfile)

