###############
# diploid set #
###############

/path_to_OrthNet_scripts/create_OrthNet.py \
	230927_diploid \
	-t /path_to_list/ \
	-T .gtfParsed.TD.txt \
	-i /path_to_CLfinderOut/230927_diploid_out.1/230927_diploid.4OrthNet.input \
	-sd -o /path_to_output/ONfinder
	
#####################
# higher ploidy set #
#####################

/path_to_OrthNet_scripts/create_OrthNet.py \
	230214_full \
	-t /path_to_list/ \
	-T .gtfParsed.TD.txt \
	-i /path_to_CLfinderOut/230914_full_out.1/230214_full.4OrthNet.input \
	-sd -o /path_to_output/ONfinder
