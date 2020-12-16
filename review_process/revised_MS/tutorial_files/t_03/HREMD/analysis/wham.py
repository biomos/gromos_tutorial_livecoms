#!/usr/bin/python
from string import *
from math import *
from array import *
import os

###################################################################################################
# read imd and dist.rest. specification file to determine restraining distances for each lambda:  #
###################################################################################################
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("imd", type=str,help="input file of HREMD, containing the lambda values of the replicas")
parser.add_argument("disres", type=str,help="distance restraints file, containing the PERTDFRESSPEC block")
args = parser.parse_args()


#read lambda values from REPLICA block in input file of HREMD:
imd_file = open(args.imd,'r')
temp_lines = False
lambda_lines = False
lambda_values = []
for line in imd_file:
  if "RETS(1" in line:
    #all lambda values should be read, we can stop now
    lambda_lines = False
  elif lambda_lines and not ("RETS(1" in line):
    #(continue) reading lambda values
    l = line.strip().split()
    for i in l:
      lambda_values.append(float(i))
  elif "RELAM(1" in line:
    #found the right spot
    lambda_lines = True
  elif temp_lines:
    l = line.strip().split()
    Temp = int(l[0])
    temp_lines = False
  elif "TEMP0" in line:
    #found the spot to read the temperature
    temp_lines = True
  

imd_file.close()

#read distance restraint specification file
disres_file = open(args.disres,'r')
PERTDFRESSPEC_block = False
R0_line = False
A_R0, B_R0 = 0.0, 0.0
Forceconstant = 0.0
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
    K_A = float(l[2])
    K_B = float(l[4])
    if K_A == K_B:
      Forceconstant = K_A
    else:
      print "ERROR: K_A is not equal K_B!"
      sys.exit()
    break

disres_file.close()

#calculate restraining distances for each lambda
rzeros = []
for L in lambda_values:
  x = (1-L)* A_R0 + L*B_R0
  rzeros.append(x)
####################################################################################################


#assuming that rlin = 1, change here if different value is used for RL in the DISTANCEFIELD block of HREMD.imd
rlin = 1.0
print "#rlin is set to 1.0, change in wham.py if different value was used for RL in the block DISTANCEFIELD of HREMD.imd"

Niterations = 5000	#number of iterations for WHAM

Nsims = len(rzeros)
Ni = 0.0
kBT = 0.00831441 * Temp
counts = []
Ubias = []
bins = []
prob = []
F = []
FEC = []


for i in range(Nsims):	#initial guess for F is zero
  F.append(0.0)
  counts.append([])
  Ubias.append([])
  for j in range(Nbins):	
    counts[i].append(0)	#set initials counts to zero
    Ubias[i].append(0)
    bins.append(0)
    prob.append(0)
    FEC.append(0)

for i in range(Nsims):
  f = open("tcf/tcf_" + str(i+1) +".out", "r")	#read in number of counts per bin per simulation
  Nbins = int(os.popen('wc -l tcf/tcf_' + str(i+1) +'.out').read().split()[0]) - 17
  for m in range(16):   #read unimportant lines
    f.readline()
  for j in range(Nbins):
    line = f.readline()
    s = split(line)
    counts[i][j] = float(s[1])
    if i == 0:
      bins[j] = float(s[0])     #read bins
      Ni = Ni + counts[i][j]    #calculate number of counts per simulation
  f.close()

for i in range(Nsims):	#calculate Ubias for all simulations and bins
  for j in range(Nbins):
    if bins[j] < (rzeros[i] - rlin):
      Ubias[i][j] = - Forceconstant * (bins[j] - rzeros[i] + 0.5 * rlin) * rlin
    elif bins[j] > (rzeros[i] + rlin):
      Ubias[i][j] = Forceconstant * (bins[j] - rzeros[i] - 0.5 * rlin) * rlin
    else:
      Ubias[i][j] = 0.5 * Forceconstant * (bins[j] - rzeros[i]) ** 2

for k in range(Niterations):	#repeat calculation of prob and F for Niterations
  for j in range(Nbins):
    denom = 0.0
    numerator = 0.0
    for i in range(Nsims):	#calculation of P(x)
      denom = denom + (Ni * exp((F[i] - Ubias[i][j]) / kBT))
      numerator = numerator + counts[i][j]
    if numerator == 0.0:
      prob[j] = 0
    else:
      prob[j] =  numerator / denom

  for i in range(Nsims):	#calculation of F[i]
    sumprob = 0
    for j in range(Nbins):
      sumprob = sumprob + prob[j] * exp(-Ubias[i][j] / kBT)
    F[i] = -kBT * log(sumprob)

outfile1 = "wham_FEC_{0}bins_{1}iter.dat".format(Nbins,Niterations)


g = open(outfile1, "w")

for j in range(Nbins):
  if prob[j] == 0:
    FEC[j] = 0
  else:
    FEC[j] = -kBT *log(prob[j])
    g.write(str(bins[j]) + "     " + str(FEC[j]) + "\n")

g.close()


