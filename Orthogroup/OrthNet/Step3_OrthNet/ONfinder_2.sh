###############
# diploid set #
###############

/path_to_OrthNet_scripts/mcl_OrthNet.py \
	230927_diploid \
	-i /path_to_output/ONfinder/ \
	-o /path_to_output/mcl/ \
	-sc -w /path_to_OrthNet_scripts/weights4mcl_OrthNet.list \
	-I 1.5	
	
#####################
# higher ploidy set #
#####################

/path_to_OrthNet_scripts/mcl_OrthNet.py \
	230214_full \
	-i /path_to_output/ONfinder/ \
	-o /path_to_output/mcl/ \
	-sc -w /path_to_OrthNet_scripts/weights4mcl_OrthNet.list \
	-I 1.5
