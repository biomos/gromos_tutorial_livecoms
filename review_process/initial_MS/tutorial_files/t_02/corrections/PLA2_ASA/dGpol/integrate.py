#!/usr/bin/python3 -tt

import sys
import re
from operator import itemgetter

def getpotentials(file, scheme, CHARGES, NPBC_SLV, NPBC_VAC, PBC_SLV, PBC_VAC, FFT_LS, FFT_BM):
  CHARGES_NEW = []
  NPBC_SLV_NEW = []
  NPBC_VAC_NEW = []
  PBC_SLV_NEW = []
  PBC_VAC_NEW = []
  FFT_LS_NEW = []
  FFT_BM_NEW = []
  
  for line in file[1:]:
    templine = line.split()
    if (len(templine) == 0): #for the case that the last line is empty
      break
    CHARGES_NEW.append(templine[3])
    NPBC_SLV_NEW.append(templine[4])
    NPBC_VAC_NEW.append(templine[5])
    PBC_SLV_NEW.append(templine[6])
    PBC_VAC_NEW.append(templine[7])
    if (scheme == 'BM'):
      FFT_LS_NEW.append(templine[8])
      FFT_BM_NEW.append(templine[9])

  CHARGES.append(CHARGES_NEW)
  NPBC_SLV.append(NPBC_SLV_NEW)
  NPBC_VAC.append(NPBC_VAC_NEW)
  PBC_SLV.append(PBC_SLV_NEW)
  PBC_VAC.append(PBC_VAC_NEW)
  FFT_LS.append(FFT_LS_NEW)
  FFT_BM.append(FFT_BM_NEW)
    
  
  return CHARGES, NPBC_SLV, NPBC_VAC, PBC_SLV, PBC_VAC, FFT_LS, FFT_BM

    

def DGcalculator(lambdas, charges, conditions):
  DG=0
  for lam in range(1,len(lambdas)):
    sum=0
    for atom in range(0, len(charges[lam])):
      summand = ((float(lambdas[lam]) - float(lambdas[lam-1])) * (float(conditions[lam-1][atom]) + float(conditions[lam][atom]))/2) #numerical integration
      sum += (float(charges[len(lambdas)-1][atom])-float(charges[0][atom])) * summand #according to equ.12 in Marias net charge changes paper
    DG += sum
  return DG
      

def main():
  i=1
  LAMS=[0.0,0.2,0.4,0.6,0.8,1.0]
  FILES=[]
  FILES.append("dGslv_pbsolv_L_1.0.out")
  FILES.append("dGslv_pbsolv_L_0.8.out")
  FILES.append("dGslv_pbsolv_L_0.6.out")
  FILES.append("dGslv_pbsolv_L_0.4.out")
  FILES.append("dGslv_pbsolv_L_0.2.out")
  FILES.append("dGslv_pbsolv_L_0.0.out")

  scheme="BM"
  #while (i < len(sys.argv)-1 ):
  #  LAMS.append(sys.argv[i])
  #  i+=1
  #  FILES.append(sys.argv[i])
  #  i+=1
  #scheme = sys.argv[-1]
  
  CHARGES=[]
  NPBC_SLV = []
  NPBC_VAC = []
  PBC_SLV = []
  PBC_VAC = []
  FFT_LS = []
  FFT_BM = []
  DG_NPBC_SLV = []
  DG_NPBC_VAC = []
  DG_PBC_SLV = []
  DG_PBC_VAC = []
  DG_FFT_LS = []
  DG_FFT_BM = []
  
  for filename in FILES:
    try:
      file = open(filename)
      state = file.readlines()
    except IOError:
      print("\t Hey there, I miss the following file: " + filename + "\n")
      sys.exit()
      
    CHARGES, NPBC_SLV, NPBC_VAC, PBC_SLV, PBC_VAC, FFT_LS, FFT_BM = getpotentials(state, scheme, CHARGES, NPBC_SLV, NPBC_VAC, PBC_SLV, PBC_VAC, FFT_LS, FFT_BM)
    
  DG_NPBC_SLV = DGcalculator(LAMS, CHARGES, NPBC_SLV)
  DG_NPBC_VAC = DGcalculator(LAMS, CHARGES,  NPBC_VAC)
  DG_PBC_SLV = DGcalculator(LAMS, CHARGES, PBC_SLV)
  DG_PBC_VAC = DGcalculator(LAMS, CHARGES, PBC_VAC)
  if scheme == 'BM':
    DG_FFT_LS = DGcalculator(LAMS, CHARGES, FFT_LS)
    DG_FFT_BM = DGcalculator(LAMS, CHARGES, FFT_BM)
    DG_FFT = DG_FFT_LS - DG_FFT_BM #actually, it should be BM-LS, but this is done only for PBC, and this enters into DG_POL as - DG_PBC, it gets -(BM-LS) = LS - BM 
  DG_NPBC = DG_NPBC_SLV - DG_NPBC_VAC
  DG_PBC = DG_PBC_SLV - DG_PBC_VAC
  DG_POL = DG_NPBC - DG_PBC
  print('\n********RESULTS NPBC********\n')
  print('DG_NPBC_SLV [kJ/mol]: ' + str(round(DG_NPBC_SLV, 3)) + '\n')
  print('DG_NPBC_VAC [kJ/mol]: ' + str(round(DG_NPBC_VAC, 3)) + '\n')
  print('DG_NPBC [kJ/mol]: ' + str(round(DG_NPBC, 3)) + '\n')
  print('\n********RESULTS PBC********\n')
  print('DG_PBC_SLV [kJ/mol]: ' + str(round(DG_PBC_SLV, 3)) + '\n')
  print('DG_PBC_VAC [kJ/mol]: ' + str(round(DG_PBC_VAC, 3)) + '\n')
  print('DG_PBC [kJ/mol]: ' + str(round(DG_PBC, 3)) + '\n')
  if scheme == 'BM':
    print('\n********RESULTS FFT PBC********\n')
    print('DG_FFT_LS [kJ/mol]: ' + str(round(DG_FFT_LS, 3)) + '\n')
    print('DG_FFT_BM [kJ/mol]: ' + str(round(DG_FFT_BM, 3)) + '\n')
    print('DG_FFT [kJ/mol]: ' + str(round(DG_FFT, 3)) + '\n')
    print('\n********RESULTS NPBC - PBC********\n')
  if scheme == 'LS':
    print('DGpol [kJ/mol]: ' + str(round(DG_POL, 3)) + '\n')
  if scheme == 'BM':
    print('DGpol [kJ/mol]: ' + str(round(DG_POL + DG_FFT, 3)) + '\n')
  print('**********************\n\n')


  
if __name__ == '__main__':
  main()
  
