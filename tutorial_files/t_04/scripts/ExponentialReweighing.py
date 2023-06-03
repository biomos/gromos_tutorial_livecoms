import numpy as np
import subprocess
from argparse import ArgumentParser


def cmdlineparse():
    parser = ArgumentParser(description="command line arguments")
    parser.add_argument("--cv1", dest="cv1", required=True, help="First CV time series")
    parser.add_argument("--cv2", dest="cv2", required=True, help="Second CV time series")
    parser.add_argument("--vr", dest="vr", required=True, help="Energy of the biased trajectory")
    parser.add_argument("--vy", dest="vy", required=True, help="Energy of the unbiased trajectory")
    parser.add_argument("--xdim", dest="xdim", required=True, nargs="+", help="Dimensions and space between bins for the first CV")
    parser.add_argument("--ydim", dest="ydim", required=True, nargs="+", help="Dimensions and space between bins for the second CV")
    args=parser.parse_args()
    return args

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

def create_bin_map(xmin, xmax, dimx, ymin, ymax, dimy):
    bin2value = {}
    binsx = np.array([i for i in np.arange(xmin, xmax, dimx)])
    binsy = np.array([i for i in np.arange(ymin, ymax, dimy)])
    count = 0
    # create mapping for the linearized bins
    for e in binsx:
        for ee in binsy:
           bin2value.setdefault(count, (e,ee))
           count += 1
    return binsx, binsy, bin2value

def bin_data(binsx, binsy, datax, datay):
    binedx = np.digitize(datax, binsx)
    binedy = np.digitize(datay, binsy)
    # substract 1 to bins to count from 0
    binedx -= 1
    binedy -= 1
    # slide first bin to be sapced in increments of len(binsY)
    binedx *= len(binsy)
    # add them together
    lbins = binedx + binedy
    return lbins

def create_reweight_input(numbins, vr, vy):
    with open("gromos_reweight_2d.arg", "w") as inn:
        inn.write("@temp 300\n")
        inn.write("@x linear_bins.dat\n")
        inn.write("@vr %s\n" % vr)
        inn.write("@vy %s\n" % vy)
        maxb = numbins - 0.5
        # to make rounded bins
        inn.write("@bounds -0.5 %s %s\n" % (maxb, numbins))

def create_bined_timeseries(bined_data, time):
    with open("linear_bins.dat","w") as inn:
        for i in range(len(bined_data)):
            inn.write("%s  %s\n" % (time[i], bined_data[i]))

def probins_to_energy(lbins):
    enebins = []
    # set max ene to 1/10th of the rarest event
    maxene = - np.log(np.min(lbins[np.where(lbins != 0)[0]])/10.0) * 300 * 8.314462618*10**-3
    for i,e in enumerate(lbins):
        if e == 0.0:
            prob = maxene
        else:
            prob = - np.log(e) * 300 * 8.314462618*10**-3
        enebins.append(prob)
    enebins = np.array(enebins)
    enebins -= min(enebins)
    return enebins

def delinearize_bins(ener_bins, bins2value, output_file):
    with open(output_file, "w") as out:
        for i, ene in enumerate(ener_bins):
            out.write("%s\t%s\t%s\n" % (bins2value[i][0], bins2value[i][1], ene))

def run_reweight():
    with open("gromos_reweight_output.out", "w") as out:
        subprocess.run(["reweight", "@f", "gromos_reweight_2d.arg"], stdout=out)
       
def main():
    # read cmd arguments
    args = cmdlineparse()
    # load data
    time, cv1dat = read_gromos_output_file(args.cv1)
    time, cv2dat = read_gromos_output_file(args.cv2)
    #vrdat  = read_gromos_output_file(args.vr)
    #vydat  = read_gromos_output_file(args.vy)
    # bin data
    binsx, binsy, bin2value = create_bin_map(float(args.xdim[0]), float(args.xdim[1]), float(args.xdim[2]),
                                             float(args.ydim[0]), float(args.ydim[1]), float(args.ydim[2]))
    linearbins = bin_data(binsx, binsy, cv1dat, cv2dat)
    # create reweight input file
    numbins = len(binsx) * len(binsy)
    create_reweight_input(numbins, args.vr, args.vy)
    create_bined_timeseries(linearbins, time)
    # run reweight
    run_reweight()
    _, reweighted_bins = read_gromos_output_file("gromos_reweight_output.out")
    # turn probabilities into energies
    ener_bins = probins_to_energy(reweighted_bins)
    # delinearize bins and save them to a file
    delinearize_bins(ener_bins, bin2value, "reweighted_PMF_2D.dat")
    # clean up data
    
if __name__ == "__main__":
    main()
