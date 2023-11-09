###############
# diploid set #
###############

/path_to_OrthNet_scripts/CL_finder_multi.py \
	230927_diploid.list \
	-t /path_to_list/ \
	-unrp -T .gtfParsed.TD.txt \
	-b /path_to_BestHitPairs/BHPairs.1/ \
	-W 20 -N 3 -G 20 -o /path_to_BestHitPairs/CLfinder/230927_diploid_out.1

#####################
# higher ploidy set #
#####################

/path_to_OrthNet_scripts/CL_finder_multi.py \
	230214_full.list \
	-t /path_to_list/ \
	-unrp -T .gtfParsed.TD.txt \
	-b /path_to_BestHitPairs/BHPairs.1/ \
	-W 20 -N 3 -G 20 -o /path_to_BestHitPairs/CLfinder/230914_full_out.1
