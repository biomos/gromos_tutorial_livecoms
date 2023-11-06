#!/bin/bash
cd ./ene_ana
ene_ana @f ene_ana_meoh.arg > ene_ana_meoh.out
cd ../rmsd
rmsd @f rmsd_meoh.arg > rmsd_meoh.out
cd ../rdf
rdf @f rdf_ob_ow.arg > rdf_ob_ow.out
cd ../tser
tser @f tser_meoh.arg > tser_meoh.out
cd ../hbond
hbond @f hbond_meoh.arg > hbond_meoh.out