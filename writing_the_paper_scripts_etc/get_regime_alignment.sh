#!/bin/bash

working_directory="/Users/fayweili/Desktop/mockdata/6_runsWchimeraCountsOutput"
cd $working_directory

rm -r All_* 
for locus in APP GAP IBR PGI
	do
		# copy all, e.g., APP_clustered_reconsensus.afa to the All_APP directory, as R2*_APP_clustered_reconsensus.afa
		mkdir All_${locus}
		for regime in R2a R2c R2e # modify here to add more regime
			do
				# copy
				cp ${regime}/${locus}_clustered_reconsensus.afa All_${locus}/${regime}_${locus}_clustered_reconsensus.afa
				# rename seq names to reflect regime 
				sed "s/>/>${regime}/" All_${locus}/${regime}_${locus}_clustered_reconsensus.afa > All_${locus}/${regime}_${locus}_clustered_reconsensus_renamed.afa
			done 
		
		# do the iterative profile-profile alignment; hardcoded here
		muscle -profile -in1 All_${locus}/R2a_${locus}_clustered_reconsensus_renamed.afa -in2 All_${locus}/R2c_${locus}_clustered_reconsensus_renamed.afa -out All_${locus}/R2ac_${locus}_clustered_reconsensus_renamed.afa
		muscle -profile -in1 All_${locus}/R2ac_${locus}_clustered_reconsensus_renamed.afa -in2 All_${locus}/R2e_${locus}_clustered_reconsensus_renamed.afa -out All_${locus}/R2ace_${locus}_clustered_reconsensus_renamed.afa
	done