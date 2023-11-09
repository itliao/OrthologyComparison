###############
# diploid set #
###############

/path_to_OrthNet_scripts/update_BestHitPairs.py \
	230927_diploid.list \
	/path_to_BestHitPairs/CLfinder/230927_diploid_out/230927_diploid.4OrthNet.input \
	-b /path_to_BestHitPairs/BHPairs/ \
	-o /path_to_BestHitPairs/BHPairs.1/
	
#####################
# higher ploidy set #
#####################

/path_to_OrthNet_scripts/update_BestHitPairs.py \
	230214_full.list \
	/path_to_files/230214_full.4OrthNet.input \
	-b /path_to_BestHitPairs/BHPairs/ \
	-o /path_to_BestHitPairs/BHPairs.1/
