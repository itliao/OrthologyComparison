# Categorizing ortholog relationships between two species

* script: 
	* CountingOrthologs_FINAL.ipynb
	
* inputs:
	* baseline - OrthoFinder-BLAST-MCL on every single species pair possible
	  * [link to zip file](/baseline_files/baseline_pairwise.zip)
	* orthogroups from diploid set and higher ploidy set
	  * [see Orthogroup_README.md](https://github.com/itliao/OrthologyComparison/blob/main/Orthogroup/Orthogroup_README.md#output-files---to-use-for-inputs-for-downstream-analyses)
	
* raw output tables: 
	* OF_b - dip_OF_blast_orthoCat_230927.csv, full_OF_blast_orthoCat_230213.csv
	* OF_d - dip_OF_diamond_orthoCat_230927.csv, full_OF_diamond_orthoCat_230213.csv
	* OF_m - dip_OF_mmseqs_orthoCat_230927.csv, full_OF_mmseqs_orthoCat_230213.csv
	* SP_d - dip_SP_diamond_orthoCat_230927.csv, full_SP_diamond_orthoCat_230119.csv
	* SP_m - dip_SP_mmseqs_orthoCat_230927.csv, full_SP_mmseqs_orthoCat_230119.csv
	* BR - dip_BR_orthoCat_230927.csv, full_BR_orthoCat_230927.csv
	* ON - dip_ON_orthoCat_230927.csv, full_ON_orthoCat_230915.csv

	* baselineOF_blast_230221.csv
	
* reformatted output tables 
	* combinedOrthologCount_230927.xlsx
	* combinedOrthologCount_reformat_230928.csv
		* for making plots
	
# pairwise gene comparison between algorithms - for species pairs

* script: 	
	* CompareOrthogroupComposition_SpeciesPairs_dip_FINAL.py
	* CompareOrthogroupComposition_SpeciesPairs_mix_FINAL.py
	
* output:
	* diploid_speciesPairsJI_230929.csv
	* mixed_speciesPairsJI_230929.csv
