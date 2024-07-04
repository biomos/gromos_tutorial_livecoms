import glob
import os
import math
import numpy as np
import pandas as pd

"""
Script the calculates:
    Number of contributing frames for each endstates based on energy: CONT_FRAMES
    Percentatge of frames contributing: PERCENTATGE
    Number of frames that the given endstate was the lowest energy one: CountMIN
    Percentatge of frames in which the endstate was the lowest energy one: PERCEN_MIN
    Number of frames that contribute to the free energy of that given endstate: DG

The variable "OFFSETS" needs to be adjusted to the offsets used on the run in order
The variable "FREE" needs to be adjusted to the free energy difference of each endstate to the reference endstate
The variable "REFERENCE" needs to b the file containing the time series of the reference hamiltonian
The script expects the energy timeseries of the diferent endstates to be named: e1.dat, e2.dat, ....
"""


# CONSTANT VARIABLES
temp = 300
boltzman = 0.00831441
BETA = 1.0/(temp * boltzman)
#OFFSETS = [0, 0] 

df = pd.read_table("df.out", names=["DF", "err"], index_col=0, delim_whitespace=True, header=0)
refstates = df.index.str.endswith("_R")
FREE = np.array(df[refstates]['DF'].values)

REFERENCE = "eds_vr.dat"

def read_energy_file(file_name):
    values_array = []
    with open(file_name) as in_file:
        #skip first line
        in_file.readline()
        for line in in_file:
            line = line.rstrip()
            value = float(line.split(" ")[-1]) 
            values_array.append(value)
    return np.array(values_array, dtype=np.float64)

def read_offset_file(file_name):
    read_on = 0
    with open(file_name) as in_file:
        #skip first line
        in_file.readline()
        for line in in_file:
            while read_on == 0:
                line = line.rstrip()
                value = float(line.split(" ")[-1]) 
                read_on = 1
    return value

def get_offsets():
    offset_files = glob.glob("e*[0-9]r.dat")
    offset_files = sorted(offset_files, key=lambda x: int(x.split(".")[0][1:-1]), reverse=False)
    offsets = [read_offset_file(x) for x in offset_files]
    return offsets


def main():
    #read files
    OFFSETS = get_offsets()
    endstates_files = glob.glob("e*[0-9].dat")
    endstates_files = sorted(endstates_files, key=lambda x: int(x.split(".")[0][1:]), reverse=False)
    endstates_e = [read_energy_file(x) for x in endstates_files]
    endstates_totals = np.array([], dtype=np.float64)
    enes = np.array([], dtype=np.float64)
    reference = read_energy_file(REFERENCE)
    #compute the energies for each endstate
    for i,hi in enumerate(endstates_e):
        #compute exponential term
        de = (hi - OFFSETS[i]) * BETA * -1.0
        #compute exp energy summation
        expde = np.exp(de)
        endstates_totals = np.append(endstates_totals, expde)
        enes = np.append(enes, hi)
    #normalize results
    #format array
    n_states = len(endstates_files)
    n_frames = int(len(endstates_totals)/n_states)
    endstates_totals = endstates_totals.reshape(n_states, n_frames)
    enes = enes.reshape(n_states, n_frames)
    contributions = {x:0.0 for x in range(len(endstates_files))}
    lowest_energy = {x:0.0 for x in range(len(endstates_files))}
    dG_diff = {x:len(np.where((enes[x] - reference) < (FREE[x] + (temp * boltzman)))[0]) for x in range(len(endstates_files))}
    # compute contributions per frame
    for i in range(n_frames):
        final_e = np.sum(endstates_totals[:,i])
        lowest_energy[np.where(endstates_totals[:,i] == np.max(endstates_totals[:,i]))[0][0]] += 1.0
        for j in range(n_states):
            contributions[j] += endstates_totals[j,i]/final_e
    #print results
    tot_con = 0.0
    tot_con_2 = 0.0
    tot_con_3 = 0.0
    tot_dG = 0.0
    for key in contributions:
        tot_con += contributions[key]
        tot_con_2 += lowest_energy[key]
        tot_dG += dG_diff[key]
    print(f"ENDSTATE\tCONT_FRAMES\tPERCENTATGE\tDG")
    for i in range(n_states):
        print(f"Endstate_{i+1}\t{round(contributions[i],2)}\t\t{round(contributions[i]*100/tot_con,2)}\t\t{dG_diff[i]}")

if __name__ == "__main__":
    main()
    
        

