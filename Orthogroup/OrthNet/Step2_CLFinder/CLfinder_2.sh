###############
# diploid set #
###############

# 230927 - diploid set + Aar 
/path_to_OrthNet_script/CL_finder_multi.py \
	230927_diploid.list \
	-t /path_to_list/ \
	-nr -T .gtfParsed.TD.txt \
	-b /path_to_BestHitPairs/BHPairs/ \
	-W 20 -N 3 -G 20 -o /path_to_BestHitPairs/CLfinder/230927_diploid_out
	
#####################
# higher ploidy set #
#####################
/path_to_OrthNet_script/CL_finder_multi.py \
	230214_full.list \
	-t /path_to_list\
	-nr -T .gtfParsed.TD.txt \
	-b /path_to_BestHitPairs/BHPairs/ \
	-W 20 -N 3 -G 20 -o /path_to_BestHitPairs/CLfinder/230214_full_out
	
