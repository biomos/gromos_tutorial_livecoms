#!/usr/bin/python3

import sys
import re
import math
import os

def truncate(n, decimals=0):
  multiplier = 10 ** decimals
  return int(n * multiplier) / multiplier

def read_topo(toposrc):
  #open topology
  topofile = open(toposrc, 'r')
  topolines = topofile.readlines()
  topofile.close()
  return topolines


def interpolate(TOPO_A, TOPO_B, LAM):
  ##########################################
  # loop through the lines of the topologies
  # if we're in the SOLUTEATOM BLOCK
  # check if the two charges are the same
  # if not, interpolate according to the gi-
  # ven LAM point
  ##########################################
  read=0 #determines if line should be read
  START="SOLUTEATOM"
  END="END"
  NEWTOPO = []
  for linenum in range(len(TOPO_A)):
    #write line from topo B to the NEWTOPO
    NEWTOPO.append(TOPO_A[linenum])
    if START in TOPO_A[linenum]:
      read=1
    if read is 1:
      if END in TOPO_A[linenum]:
        read=0
      elif not "#" in TOPO_A[linenum]:
        # check if we're in the line that actually contains charges
        if len(TOPO_A[linenum].split()) >= 8:
          # check if the two charges are the same
          if TOPO_A[linenum].split()[5] != TOPO_B[linenum].split()[5]:
            interpol_charges=float(TOPO_A[linenum].split()[5]) * (1-float(LAM)) + float(TOPO_B[linenum].split()[5]) * float(LAM)

            # set the strings to substitute and subsitute with such that the lenght of the strings is always the same - independent of the leading minus sign
            if ( float(TOPO_A[linenum].split()[5]) < 0 and interpol_charges < 0 ) or ( float(TOPO_A[linenum].split()[5]) >= 0 and interpol_charges >= 0 ):
              subsitute_=TOPO_A[linenum].split()[5]
              with_=format(interpol_charges,'.5f')
            elif float(TOPO_A[linenum].split()[5]) < 0 and interpol_charges >= 0:
              subsitute_=TOPO_A[linenum].split()[5]
              with_=format(interpol_charges,' .5f')
            elif float(TOPO_A[linenum].split()[5]) >= 0 and interpol_charges < 0:
              subsitute_=' ' + TOPO_A[linenum].split()[5]
              with_=format(interpol_charges,'.5f')
            # change the charge in the last line of the NEWTOPO
            NEWTOPO[-1] = re.sub(subsitute_,with_, NEWTOPO[-1])
              

  return NEWTOPO

          
          



def main():
  if len(sys.argv) != 4:
    print("# This script can interpolate charges between two given topologies")
    print("# it takes topologies at state A and B and a specific lambda-value")
    print("# it writes out the new topology to std.out")
    print("# note that all other (e.g. Vdw-)parameters are taken from TOPO_A")
    print()
    print("# usage: interpolate_topocharges.py TOPO_A TOPO_B LAM")
    sys.exit()
  
  TOPO_PATH_A = sys.argv[1]
  TOPO_PATH_B = sys.argv[2]
  LAM = sys.argv[3]
  
  TOPO_FILE_A = read_topo(TOPO_PATH_A)
  TOPO_FILE_B = read_topo(TOPO_PATH_B)

  TOPO_FILE_NEW = interpolate(TOPO_FILE_A, TOPO_FILE_B, LAM)

  for topoline in TOPO_FILE_NEW:
    print(topoline, end='')
  
if __name__ == '__main__':
  main()
