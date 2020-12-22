#!/bin/bash

...your-path-to... /programs/bar \
@files \
../../md_TI/md_TI_PLA2_ASA/L_0.0/bar_* \
../../md_TI/md_TI_PLA2_ASA/L_0.1/bar_* \
../../md_TI/md_TI_PLA2_ASA/L_0.2/bar_* \
../../md_TI/md_TI_PLA2_ASA/L_0.3/bar_* \
../../md_TI/md_TI_PLA2_ASA/L_0.4/bar_* \
../../md_TI/md_TI_PLA2_ASA/L_0.5/bar_* \
../../md_TI/md_TI_PLA2_ASA/L_0.6/bar_* \
../../md_TI/md_TI_PLA2_ASA/L_0.7/bar_* \
../../md_TI/md_TI_PLA2_ASA/L_0.8/bar_* \
../../md_TI/md_TI_PLA2_ASA/L_0.9/bar_* \
../../md_TI/md_TI_PLA2_ASA/L_1.0/bar_* \
@temp 298.0	  \
@maxiterations 500 \
@convergence 1E-10 \
@bootstrap 100 \
> \
bar_PLA2_ASA.dat				   


