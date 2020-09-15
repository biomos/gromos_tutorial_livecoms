#!/bin/bash

for LAM in 0.0 0.2 0.4 0.6 0.8 1.0
do
  echo "Processing L = ${LAM}"
  ./interpolate_topocharges.py ../../../topo/ASA_Na.top  ../../../topo/DUM_Na.top ${LAM} > ASA_Na_L_${LAM}.top
done

