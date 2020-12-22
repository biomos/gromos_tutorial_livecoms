#!/usr/bin/python3 -tt
import sys
import re
import math

# define some global variables
# elementary charges per Coulomb [e / C]
charges = 6.242 * 10**18
#avogadro's, [mol^-1]
NA    = 6.022140857 * 10**(23) 
#vacuum permittivity, [e**2 s**-2 nm**-3 kg**-1]
eps0  = 8.854187817 * 10**-12 * charges**2 
#quadrup-mom trace of the water mod rel to its single vdW interact, [e m2]
gamW  = 0.0082 *10**(-18)       
#charge, [e]
#QG    = float(netcharge)                      
#cut-off distance [m]
#cut    = float(cut) * 10**(-9)         

###########################################################
# this script integrates over the rdf output file 'rdf.out'
# the returned number still needs to be multiplicated by
# the water number density in the box (water per cubic-nm)
###########################################################

##################################################
# simple function for reading in files by filename
def readfile(filename):
  try:
    f = open(filename, 'rU')
    file = f.read()
  except IOError:
    print("\t Hey there, I miss a file: " + filename + "\n")
    sys.exit()
  return file

########################################################
# integrates over the radial distribution function
# the output of this function times the total number of
# water molecules in the box is the number of water
# molecules within the cutoff sphere
def integrate(g):

  N = 0
  dr = float(g[0][0]) #the change of the radius dr (r+-dr) is just written in the first output

  for i in range(0, len(g)):
    r = float(g[i][0])
    g_r = float(g[i][1])
    n = g_r * ((r+dr)**3 - (r-dr)**3)
    N = N + n
  return N * 4./3 * math.pi

def volume(xyz):
  V = 1
  for l in xyz:
    V = V*float(l)
  return V

def dsm_correction(N,nsolv,V,q,epsBW,rcut):
  DG = 0
  DG = - NA / (6*eps0) * 2*(epsBW-1)/(2*epsBW+1) * gamW * q * N*nsolv/V / (rcut**3 * 4./3 * math.pi)
  return DG

def main():
  rdf_file = readfile("./rdf.out")
  info_file = readfile("./system.info")
  g = re.findall(r'(\d+\.\d+)\s+(\d+.*\d*)', rdf_file)
  xyz = re.findall(r'\(x,y,z\)=\((\d+\.\d+),(\d+\.\d+),(\d+\.\d+)\)', info_file)
  nsolv = re.findall(r'n=(\d+)', info_file)
  q = re.findall(r'q=([-+]?\d+)', info_file)
  rcut = re.findall(r'rcut=(\d+\.?\d*)', info_file)
  epsBW = re.findall(r'epsBW=(\d+\.\d*)', info_file)
  gamW = re.findall(r'quadrupole=(\d+\.\d*)', info_file)
  # clean up readouts, save as floats
  xyz = xyz[0] # save tupel with vector as vector
  nsolv = float(nsolv[0])
  q = float(q[0])
  rcut = float(rcut[0])
  epsBW = float(epsBW[0])
  gamW = float(gamW[0])
  
  # do basic calculations
  N = integrate(g)
  V = volume(xyz)

  #converstions from nm to m
  rcut = rcut * pow(10,-9) #nm
  gamW = gamW * pow(10,-18) #nm^2

  # calculate the actual correction
  DG = dsm_correction(N,nsolv,V,q,epsBW,rcut)
  print("DGdsm [kJ/mol]: " + str(round(DG/1000, 3)))


if __name__ == '__main__':
  main()
  
