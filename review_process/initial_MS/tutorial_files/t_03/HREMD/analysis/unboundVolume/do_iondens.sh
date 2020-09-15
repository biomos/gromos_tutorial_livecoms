#!/bin/bash


traj_files=`eval ls rep*/unbound_region_frames.trj`

#remove any old results
rm -f grid.pdb

iondens\
	@topo 		../../../topo/PLA2_ASA_Ca_2Na.top \
	@pbc		r cog \
	@grspace 	0.1 \
	@ions 		'va(-1,2:4,13)' \
	@thresholds 	0.00005 0.00005 \
	@ref 		../../../eq/eq_PLA2_ASA_Ca_2Na_7.cnf \
	@atoms 		1:a \
	@traj 		$traj_files \
	> iondens_unbound_region.out 

Vunb=`wc -l grid.pdb |awk '{print $1}'`
echo "nr of lines in grid.pdb: $Vunb "
