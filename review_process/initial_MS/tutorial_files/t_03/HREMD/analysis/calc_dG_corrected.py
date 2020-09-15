#!/usr/bin/python

from math import *
import sys
import argparse

parser = argparse.ArgumentParser(fromfile_prefix_chars='@')
parser.add_argument("wham", type=str,help="result from wham, with minimum set to 0")
parser.add_argument("start_b",type=float,help="start of bound region in nm, e.g. 0.0")
parser.add_argument("end_b",type=float,help="end of bound region in rnm, e.g. 1.0")
parser.add_argument("start_u",type=float,help="start of unbound region in nm, e.g. 3.0")
parser.add_argument("end_u",type=float,help="end of unbound region in nm, e.g. 5.0")
parser.add_argument("Vunb",type=float,help="sampled volume in the unbound region in ")
args = parser.parse_args()

############################################

#standard state volume
V0 = 1.6630

#define kBT
kBT = 0.00831441 * 298


#file with FEC for the unbound region
f=open(args.wham,'r')

dist_u, dist_b = [], []
Exp_u, Exp_b  = [], []
#calculate exp(-WR(Z)/kBT) for each unbound Z
for line in f:
  values = line.split()
  distX = float(values[0])
  #only use if in unbound region
  if (distX >= args.start_u) and (distX <= args.end_u):
    dist_u.append(distX)
    WRZ_u = float(values[1])
    Exp_u.append(exp(-WRZ_u/kBT))
  elif (distX >= args.start_b) and (distX <= args.end_b):
    dist_b.append(distX)
    WRZ_b = float(values[1])
    Exp_b.append(exp(-WRZ_b/kBT))

int_u = 0.0
int_b = 0.0
#integrate over ubound region
for i in range(len(dist_u)-1):
  int_u += 0.5 * (Exp_u[i+1] + Exp_u[i]) * (dist_u[i+1] - dist_u[i])
#integrate over bound region
for i in range(len(dist_b)-1):
  int_b += 0.5 * (Exp_b[i+1] + Exp_b[i]) * (dist_b[i+1] - dist_b[i])


dGbind_raw = -kBT * log(int_b/int_u)
dGstd = -kBT * log(args.Vunb/V0)
dGbind_std = dGbind_raw + dGstd
print 'dGbind_raw:\t{0:>8.1f}'.format(dGbind_raw)
print 'dGstd:\t\t{0:>8.1f}'.format(dGstd)
print 'dGbind_std:\t{0:>8.1f}'.format(dGbind_std)
