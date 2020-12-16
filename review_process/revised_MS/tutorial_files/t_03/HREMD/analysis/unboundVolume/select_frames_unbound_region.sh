#!/bin/bash

if [ "$#" -ne 4 ]; then
    echo "Illegal number of parameters: specify nr of replicas, nr of runs, and the range of the unbound region (in nm) e.g.:"
    echo "# ./select_frames_unbound_region.sh 24 25 3.0 5.0 "
    exit 2
fi


#for replica numbers in range of the unbound region
nr_replicas=$1
nr_runs=$2
start_unb=$3
end_unb=$4

for j in $(seq 1 $nr_replicas);

do
  mkdir -p rep${j}
  #filter df_dist file, such that data is written out for every 2 ps (as in the *trc.gz trajectory)
  sed -n '2~20p' ../df_dist/df_dist_${j}.dat > rep${j}/df_dist.dat

  cd rep${j}
  pwd
  rm -f list_frames_unbound_region.txt

  #select frames with df distance between start_unb and end_unb
  #writes out frame numbers (df distance at 10ps comes from frame number 6 in the *trc.gz files)
  frames=`awk -v start="$start_unb" -v end="$end_unb" ' $2 >= start && $2 <= end {ORS=" "; print ($1/2) +1;}' df_dist.dat`

  #write out framenumbers which correspond to df distances within the unbound region
  if [ -n "$frames" ]; then
    echo "run"$i $frames >> list_frames_unbound_region.txt

    #make a list of all trajectory files for a single replica
    traj_files=`for xx in $(seq 1 $nr_runs); do echo -n "../../../run${xx}/HREMD_${xx}_${j}.trc.gz "; done`

    #run frameout on the selected frames
    frameout @pbc r cog @topo ../../../../topo/PLA2_ASA_Ca_2Na.top @spec SPEC @frames $frames @single @traj $traj_files
    mv FRAME*.cnf unbound_region_frames.trj
    cd ..

  #clean-up; if no frames contribute, remove directory
  else
    cd ..
    rm -r rep${j}
  fi
done

