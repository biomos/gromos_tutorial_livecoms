# This scripts prepares the files in the format that the Pyreweight module of Mio lab wants
import numpy as np
# INPUT FILE
xdata =  "phi.out" # Coordinate for the first dimension of the PMF
ydata =  "psi.out" #Coordinate for the second dimension of the PMF (optional)
zdata = "" # Coordinate for the third dimension of the PMF" # Total biasing energy of the GAMD
totene = "totgamd_dV.dat" # Biasing energy of GAMD
temp = 300
boltzman = 0.00831441
unit_conv = 0.239006 # for results in kcal

def read_gromos_output_file(in_file):
    time_series = []
    value_series = []
    with open(in_file, "r") as inn:
        for line in inn:
            line = line.rstrip()
            if not line.startswith("#"):
                fields = line.split()
                time = round(float(fields[0]),3)
                value = float(fields[1])
                time_series.append(time)
                value_series.append(value)
    return np.array(time_series), np.array(value_series)

def create_data_table(value_series, time_series, name):
    with open(name, "w") as out:
        for i in range(len(value_series)):
            out.write("%s %s %s\n" % (str(value_series[i]), 
                                      i,
                                      str(value_series[i]*unit_conv)))

def match_time_series(*args):
    # args are time + value series
    # get indexes of common elements in timeseries
    pairs = []
    for i in range(0,len(args),2):
        pairs.append([args[i], args[i+1]])

    for i  in range(len(pairs)):
        for j in range((i+1), len(pairs)):
            # change pair1
            common = np.nonzero(np.in1d(pairs[i][0], pairs[j][0]))[0]
            pairs[i][0] = pairs[i][0][common]
            pairs[i][1] = pairs[i][1][common]
            # change pair2
            common = np.nonzero(np.in1d(pairs[j][0], pairs[i][0]))[0]
            pairs[j][0] = pairs[j][0][common]
            pairs[j][1] = pairs[j][1][common]

    return pairs

def save_data2bin(data, name):
    # save bin data
    with open(name, "w") as inn:
        for e in data:
            inn.write("%s\n" % e) 
   
def save_data2bincomb(data, data2, name):
    # example funciton to save bin data of two arrays by columns
    with open(name, "w") as inn:
        for e,ee  in zip(data, data2):
            inn.write("%s %s\n" % (e,ee))  
        
def main(xdata, ydata, zdata, totene, temp, boltzman, unit_conv):

    if xdata:
        tx, dx = read_gromos_output_file(xdata)
    if ydata:
        ty, dy = read_gromos_output_file(ydata)
    if zdata:
        tz, dz = read_gromos_output_file(zdata)
    if totene:
        te, de = read_gromos_output_file(totene)

    if xdata and ydata and zdata and totene:
        data_pairs = match_time_series(tx,dx,ty,dy,tz,dz,te,de)
        save_data2bin(data_pairs[0][1], "binX.dat")
        save_data2bin(data_pairs[1][1], "binY.dat")
        save_data2bin(data_pairs[2][1], "binZ.dat")
        create_data_table(data_pairs[3][1], data_pairs[3][0], "weights.dat")

    if xdata and ydata and totene:
        data_pairs = match_time_series(tx,dx,ty,dy,te,de)
        save_data2bin(data_pairs[0][1], "binX.dat")
        save_data2bin(data_pairs[1][1], "binY.dat")
        save_data2bincomb(data_pairs[0][1], data_pairs[1][1], "binXY.dat")
        create_data_table(data_pairs[2][1], data_pairs[2][0], "weights.dat")

    if xdata and totene:
        data_pairs = match_time_series(tx,dx,te,de)
        save_data2bin(data_pairs[0][1], "binX.dat")
        create_data_table(data_pairs[1][1], data_pairs[1][0], "weights.dat")
    else:
        print("Xdata and totene should be provided")

if __name__ == "__main__":
    main(xdata, ydata, zdata, totene, temp, boltzman, unit_conv)

    

