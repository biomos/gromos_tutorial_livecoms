#!/usr/bin/python3 -tt
import sys

###########################################################
# this script integrates the numbers from ener output files
# it is very basic and everything is basically hardcoded,
# since it is only intended for GROMOS tutorial use
###########################################################

##################################################
# simple function for reading in files by filename
def readfile(filename):
  try:
    f = open(filename, 'rU')
    file = f.readlines()
  except IOError:
    print("\t Hey there, I miss the following file: " + filename + "\n")
    sys.exit()
  return file

########################################################
# takes two files (which are lists of various potentials
# and also takes two lambda points
# integrate() performes numerical integration over the
# lines separately
def integrate(fileA,fileB,lamA,lamB):
  pot = []
  for i in range(len(fileA)):
    if not "#" in fileA[i]:
      pot.append( ( float(fileB[i].split()[1]) + float(fileA[i].split()[1]) ) / 2 * (float(lamB)-float(lamA)) )
  return pot
      

def main():
  thisdict = {"NPBC": [],
              "PBC": []
  }
  LAMS=[]
  FILES=[]
  LAM = [0.0,0.2,0.4,0.6,0.8,1.0]
  for STATE in ['NPBC','PBC']:
    potentials = []
    for l in range(1,len(LAM)):
      fileA_name = 'ener_' + STATE + '_L_' + str(LAM[l-1]) + '.out'
      fileB_name = 'ener_' + STATE + '_L_' + str(LAM[l]) + '.out'
      fileA = readfile(fileA_name)
      fileB = readfile(fileB_name)
      if len(potentials) == 0: #first round
        potentials = integrate(fileA, fileB, LAM[l-1], LAM[l])
      else:
        potentials_new = integrate(fileA, fileB, LAM[l-1], LAM[l])
        for i in range(len(potentials)):
          potentials[i] = potentials[i] + potentials_new[i]
    thisdict[STATE] = potentials

  # calculate the final results
  results = []
  for i in range(len(thisdict["NPBC"])):
    results.append(thisdict["NPBC"][i] - thisdict["PBC"][i])
  for res in results:
    print("DGdir [kJ/mol]: " + str(round(res,3)))


if __name__ == '__main__':
  main()
  
