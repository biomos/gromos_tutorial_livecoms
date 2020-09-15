#!/bin/bash

#perform wham and move the curve such that the minimum is at 0
./wham.py ../HREMD.imd ../disres.dat
awk 'NR==1 { min=$2 } 
     FNR==NR { if ($2<=min) min=$2 ; next}
     { $2=($2-min) ; print}' wham_FEC_200bins_5000iter.dat wham_FEC_200bins_5000iter.dat > wham_FEC_200bins_5000iter_min0.dat

