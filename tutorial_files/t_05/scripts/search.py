import glob
import os
import math
import numpy as np

BLOCKS = 10
REFERENCE = "eds_vr.dat"
temp = 300
boltzman = 0.00831441
BETA = 1.0/(temp * boltzman)


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

def accelerate_hamiltonian(hi, emin, emax):
    ahi = []
    for i,e in enumerate(hi):
        if (e <= emin[i]):
            ahi.append(e)
        elif (e >= emax[i]):
            ahi.append(e - 0.5 * (emax[i] - emin[i]))
        else:
            edemix = e - emin[i]
            ekfac = 1.0 / (emax[i] - emin[i])
            ahi.append((e - (0.5 * ekfac * edemix * edemix)))
    return np.array(ahi, dtype=np.float64)

def compute_theoretical_offset(energies, emin, emax, offset, reference, chunk_size, blocks):
    theoric_offset = []
    for i in range(blocks):
        r_emin = emin[i*chunk_size:] + offset[i*chunk_size:] 
        r_emax = emax[i*chunk_size:] + offset[i*chunk_size:]
        hi = energies[i*chunk_size:]
        # accelerate hamiltonian
        ahi = accelerate_hamiltonian(hi, r_emin, r_emax)
        # calculate the exponential term
        de = (ahi - reference[i*chunk_size:]) * BETA * -1.0
        # calculate free energy difference respect to reference
        expde = np.exp(de)
        endstate_fren = np.log(np.sum(expde)/len(expde))
        endstate_fren *= (-1.0/BETA)
        theoric_offset.append(endstate_fren)
    return theoric_offset
     

def main():
    emax_array = read_energy_file("eds_emax.dat")
    emin_array = read_energy_file("eds_emin.dat")
    chunk_size = math.ceil(len(emin_array)/BLOCKS)
    chunk_size = chunk_size * (BLOCKS-1)
    mean_emax = np.mean(emax_array[chunk_size:])
    mean_emin = np.mean(emin_array[chunk_size:])
    std_emax = np.std(emax_array[chunk_size:])
    std_emin = np.std(emin_array[chunk_size:])
    print("emax_mean: %s, emax_std: %s, emin_mean: %s, emin_std: %s, " % (mean_emax,std_emax,mean_emin,std_emin))
    offsets_set = glob.glob("e*[0-9]r.dat")
    offsets_set = sorted(offsets_set, key=lambda x: x.split(".")[0][1:-1], reverse=False)
    endstates_files = glob.glob("e*[0-9].dat")
    endstates_files = sorted(endstates_files, key=lambda x: x.split(".")[0][1:], reverse=False)
    reference_array = read_energy_file(REFERENCE)
    for i,offset in enumerate(offsets_set):
        num = offset.split(".")[0][1:-1]
        offset_array = read_energy_file(offset)
        endstate_array = read_energy_file(endstates_files[i])
        # Compute theoretical offset based on blocks
        chunk_size = math.ceil(len(offset_array)/BLOCKS)
        theoric_offset = compute_theoretical_offset(endstate_array, emin_array, emax_array, offset_array, reference_array, chunk_size, BLOCKS)
        means = [np.mean(offset_array[i*chunk_size:]) for i in range(BLOCKS)]
        stds =  [np.std(offset_array[i*chunk_size:]) for i in range(BLOCKS)]
        print("Offset for %s is:" % num)
        print("block   mean   std    Theor offset")
        for i in range(len(means)):
            print("block %s : %s   %s    %s" % (i+1, round(means[i],2), round(stds[i],2), round(theoric_offset[i],2)))  
    
main()        

