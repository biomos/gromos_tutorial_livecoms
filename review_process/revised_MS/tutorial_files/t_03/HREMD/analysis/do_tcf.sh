#!/bin/bash

if [ "$#" -ne 1 ]; then
    echo "Illegal number of parameters: specify nr_replicas"
    exit 2
fi

nr_replicas=$1

mkdir -p tcf

for rep in $(seq 1 $nr_replicas); do
  tcf @files df_dist/df_dist_${rep}.dat @time 0 0.1 @distribution 2 @bounds 0 5 200 @normalize > tcf/tcf_${rep}.out
done

