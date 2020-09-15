#/bin/bash

if [ "$#" -ne 2 ]; then
    echo "Illegal number of parameters: specify nr_runs and nr_trials"
    exit 2
fi

nr_runs=$1
nr_trials=$2

head -6 ../run1/replica.dat | awk ' BEGIN {OFS="\t"} { print }' > replica_all.dat ; for i in $(seq 1 $nr_runs); do x=$(( $nr_trials*(i-1) )); tail -n +7 ../run${i}/replica.dat |awk  -v x="$x" ' BEGIN {OFS="\t"} { $3=$3+x; print } '  >> replica_all.dat ;done
