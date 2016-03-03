#!/bin/bash

# This shell script takes the output from multiple PURC clustering runs (these runs need
# to be organized in a master "working_directory" with each analysis regime getting its
# own subdirectory, named with the regime name. E.g.,
# masterdirectory/reg1
# masterdirectory/reg2  etc.

working_directory="/Users/fayweili/Desktop/mockdata/6_runsWchimeraCountsOutput"
cd $working_directory

# rm -r All_* 
# rm command deletes all the All_* things (directories, etc) within the working directory. 

for locus in APP GAP IBR PGI
	do
		# copy all, e.g., APP_clustered_reconsensus.afa to the All_APP subdirectory, as R2*_APP_clustered_reconsensus.afa
		mkdir All_${locus} # E.g., "All_APP" etc.
		for regime in R2a R2c R2e # Modify this to add regimes or change regime names
			do
				# copy the locus alignment into the new directory
				# E.g., copy R2a/APP_clustered_reconsensus.afa to All_APP, and rename it as R2a_APP_clustered_reconsensus.afa
				cp ${regime}/${locus}_clustered_reconsensus.afa All_${locus}/${regime}_${locus}_clustered_reconsensus.afa
				# rename seq names to reflect regime 
				sed "s/>/>${regime}/" All_${locus}/${regime}_${locus}_clustered_reconsensus.afa > All_${locus}/${regime}_${locus}_clustered_reconsensus_renamed.afa
			done 
		
		# Now that we have all the alignments and the sequences names are changed,
		# do the iterative profile-profile alignment; hardcoded here. If you have more regimes, you'll need to add additional lines
		# (changing the two infile names and the outfile name each time)
		muscle -profile -in1 All_${locus}/R2a_${locus}_clustered_reconsensus_renamed.afa -in2 All_${locus}/R2c_${locus}_clustered_reconsensus_renamed.afa -out All_${locus}/R2ac_${locus}_clustered_reconsensus_renamed.afa
		muscle -profile -in1 All_${locus}/R2ac_${locus}_clustered_reconsensus_renamed.afa -in2 All_${locus}/R2e_${locus}_clustered_reconsensus_renamed.afa -out All_${locus}/R2ace_${locus}_clustered_reconsensus_renamed.afa
	done