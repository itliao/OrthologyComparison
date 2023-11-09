###############
# diploid set #
###############

# 230927
/path_to_OrthNet_scripts/update_OrthNet_after_mcl.py \
	230927_diploid \
	/path_to_files/230927_pairs/mcl/230927_diploid_TD1.5_rC1.2_rNC0.5_uC0.6_uNC0.25_I1.5_mclOutPC.txt \
	-i /path_to_files/230927_pairs/ONfinder/ \
	-t /path_to_list \
	-b /path_to_files/230927_pairs/BHPairs.1/ \
	-o1 /path_to_files/230927_pairs/BHPairs.2/ \
	-o2 /path_to_files/230927_pairs/ONfinder/230927_diploid.2 \
	-u -T .gtfParsed.TD.txt -W 20 -N 3 -G 20

#####################
# higher ploidy set #
#####################

# 230915
/path_to_OrthNet_scripts/update_OrthNet_after_mcl.py \
	230214_full \
	/path_to_files/230214_pairs/mcl/230214_full_TD1.5_rC1.2_rNC0.5_uC0.6_uNC0.25_I1.5_mclOutPC.txt \
	-i /path_to_files/230214_pairs/ONfinder/ \
	-t /path_to_list \
	-b /path_to_files/230214_pairs/BHPairs.1/ \
	-o1 /path_to_files/230214_pairs/BHPairs.2/ \
	-o2 /path_to_files/230214_pairs/ONfinder/230214_full.2 \
	-u -T .gtfParsed.TD.txt -W 20 -N 3 -G 20
