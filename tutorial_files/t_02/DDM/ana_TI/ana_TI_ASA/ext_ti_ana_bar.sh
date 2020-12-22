#!/bin/bash

# variable to define lambda directories 
LAMBDA=($(seq 0 0.1 1.0))

# variable to define lambda points (for bar)
LAMBDAV=$(seq -s " " 0 0.1 1.0)

# for loop over the lambda directories
for i in "${LAMBDA[@]}"
do
	# define job number for first and last job to be included in the analysis
	START=$(less ../../md_TI/md_TI.jobs | grep L_${i} | sed '2q;d' | awk '{printf $1}') 
	END=$(less ../../md_TI/md_TI.jobs | grep L_${i} | tail -1 | head -1 | awk '{printf $1}') 

	# change to respective lambda directory
	cd ../../md_TI/md_TI_ASA/L_${i} 

	# path to program
	# trajectory files (energy and free energy) from first to last job
	# define library (version number has to match md binary!)
	# bar is optional: pass the lambda points that should be considered
	...your-path-to...  /programs/ext_ti_ana \
	    @en_files	$(eval echo *{$START..$END}.tre.gz) \
	    @fr_files	$(eval echo *{$START..$END}.trg.gz) \
	    @library	../../../ana_TI/ene_ana.md++.lib \
            @bar_data     $(eval echo ${LAMBDAV}) \
	    @imd	  *_${END}.imd \
	    @lam_precision 4 
	mv predicted* ext_TI_L_${i}.dat
	mv bar_data* bar_L_${i}.dat
	cd ../../../ana_TI/ana_TI_ASA/
 	
done

