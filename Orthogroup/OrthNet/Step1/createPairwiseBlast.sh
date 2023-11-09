###############
# diploid set #
###############
echo "creating script to then run pairwise mmseqs"

/path_to_OrthNet_scripts/create_pairwiseBLAST_commands.py \
	/path_to_files/230927_diploid.list \
	-d /path_to_files \
	-o /path_to_output/230927_pairs \
	-M -n "--max-seqs 10" \
	> /path_to_files/230927_diploid_pairwiseMMseqs2.sh
	
#####################
# higher ploidy set #
#####################

echo "creating script to then run pairwise mmseqs"

/path_to_OrthNet_scripts/create_pairwiseBLAST_commands.py \
	/path_to_files/230214_full.list \
	-d /path_to_files \
	-o /path_to_output/230214_pairs \
	-M -n "--max-seqs 10" \
	> /path_to_files/230214_full_pairwiseMMseqs2.sh
