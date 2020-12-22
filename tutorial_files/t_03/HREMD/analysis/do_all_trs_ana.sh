#!/bin/bash

mkdir -p df_dist

if [ "$#" -ne 2 ]; then
    echo "Illegal number of parameters: specify nr_runs and nr_replicas"
    exit 2
fi

nr_runs=$1
nr_replicas=$2


#for each replica i;
for i in $(seq 1 $nr_replicas); do 

  #get a list of trs files of replica i for each of the runs
  trs_files=`for xx in $(seq 1 $nr_runs); do echo -n "../run${xx}/HREMD_${xx}_${i}.trs.gz "; done`

  trs_ana \
	@trs  		$trs_files \
	@prop 		df_dist \
	@library 	trs.lib \


  mv df_dist.dat df_dist/df_dist_${i}.dat

done
