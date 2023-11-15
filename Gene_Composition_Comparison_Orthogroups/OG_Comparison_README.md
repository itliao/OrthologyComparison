# comparing gene composition in orthogroups between different orthology inference algorithms 

## calculating Jaccard Index (JI), Rand Score (RS), and Adjusted Rand Score (ARS)
* script: 
	* CompareOrthogroupComposition_diploid_FINAL.py
	* CompareOrthogroupComposition_full_FINAL.py 

* input
	* orthogroups from diploid set and higher ploidy set
		* [see Orthogroup_README.md](https://github.com/itliao/OrthologyComparison/blob/main/Orthogroup/Orthogroup_README.md#output-files---to-use-for-inputs-for-downstream-analyses)
	* at_features.tsv
	
* output
	* [for diploid set](/Gene_Composition_Comparison_Orthogroups/diploid_output):
		* diploid_pairwiseMetrics_231003.csv - **found in DRYAD**
			* to use to plot and summarize similarity indices
		* diploid_geneClusterTable_230928.csv 
			* includes ON
			* all orthogroups with Arabidopsis gene concatenated together
		* conservativeOG_diploidSet_230929.csv 
			* does not include ON orthogroup inferences
			* gene is included in orthogroup if it is found in five of six orthogroups from  different orthology inference programs
		* inclusiveOG_diploidSet_230929.csv
			* does not include ON
			* gene is included in orthogroup if it is found in at least one of the orthogroups from  different orthology inference programs
		* geneOutput_230929.zip - **found in DRYAD**
			* individual orthogroup files, named by *Arabidopsis* gene geneOutput_230929.zip
			* does not include ON
		
	* [for higher ploidy set](/higher_ploidy_output):
		* full_pairwiseMetrics_231004.csv - **found in DRYAD**
			* to use to plot and summarize similarity indices
		* full_geneClusterTable_230928.csv 
			* includes ON
			* all orthogroups with Arabidopsis gene concatenated together
		* conservativeOG_fullSet_230928.csv 
			* does not include ON orthogroup inferences
			* gene is included in orthogroup if it is found in five of six orthogroups from  different orthology inference programs
		* inclusiveOG_fullSet_230928.csv 
			* does not include ON
			* gene is included in orthogroup if it is found in at least one of the orthogroups from  different orthology inference programs
		* geneOutput_230928.zip  - **found in DRYAD**
			* individual orthogroup files, named by *Arabidopsis* gene geneOutput_230928.zip
			* does not include ON
