#!/bin/bash

for LAM in 0.0 0.2 0.4 0.6 0.8 1.0
do
  echo "Processing L = ${LAM}"
  ./interpolate_topocharges.py ../../../topo/PLA2_ASA_Ca_2Na.top ../../../topo/PLA2_DUM_Ca_2Na.top ${LAM} > PLA2_ASA_Ca_2Na_L_${LAM}.top
done

