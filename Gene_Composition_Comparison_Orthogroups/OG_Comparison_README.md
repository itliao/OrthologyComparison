# comparing gene composition in orthogroups between different orthology inference algorithms 

## calculating Jaccard Index (JI), Rand Score (RS), and Adjusted Rand Score (ARS)
* script: 
	* CompareOrthogroupComposition_diploid_240217.py
	* CompareOrthogroupComposition_full_240217.py 

* input
	* orthogroups from diploid set and higher ploidy set
		* [see Orthogroup_README.md](https://github.com/itliao/OrthologyComparison/blob/main/Orthogroup/Orthogroup_README.md#output-files---to-use-for-inputs-for-downstream-analyses)
	* at_features.tsv
	
* output
	* [for diploid set](/Gene_Composition_Comparison_Orthogroups/diploid_output):
		* diploid_pairwiseMetrics_240212.csv - **found in DRYAD**
  			* [https://doi.org/10.5061/dryad.8sf7m0cw8](https://doi.org/10.5061/dryad.8sf7m0cw8)
			* raw similarity score metrics, but includes duplicated orthogroup comparisons (see below for correction)
		* summary_diploid_240221.csv
        		* for tracking whether the correct number of duplicate comparisons were removed 		 	
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
			* [https://doi.org/10.5061/dryad.8sf7m0cw8](https://doi.org/10.5061/dryad.8sf7m0cw8)
   			* individual orthogroup files, named by *Arabidopsis* gene geneOutput_230929.zip
			* does not include ON
		
	* [for diploid+higher ploidy set](/Gene_Composition_Comparison_Orthogroups/higher_ploidy_output):
		* full_pairwiseMetrics_240209.csv - **found in DRYAD**
			* [https://doi.org/10.5061/dryad.8sf7m0cw8](https://doi.org/10.5061/dryad.8sf7m0cw8)
    			* raw similarity score metrics, but includes duplicated orthogroup comparisons (see below for correction)
		* summary_full_240221.csv
      			* for tracking whether the correct number of duplicate comparisons were removed 
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
			* [https://doi.org/10.5061/dryad.8sf7m0cw8](https://doi.org/10.5061/dryad.8sf7m0cw8)
    			* individual orthogroup files, named by *Arabidopsis* gene geneOutput_230928.zip
			* does not include ON
 
* script (for correcting overcounting/duplicate orthogroup comparisons and metrics: 
	* Correcting_Pairwise_Metrics_240220-21_cleaned.ipynb
 * output:
 	* diploid_newMetrics_noDup_240221.csv **found in DRYAD**
   		* [https://doi.org/10.5061/dryad.8sf7m0cw8](https://doi.org/10.5061/dryad.8sf7m0cw8)
      		* no duplicate orthogroup composition metrics; to use to plot and summarize similarity indices
        * full_newMetrics_noDup_240221.csv  - **found in DRYAD**		
		* [https://doi.org/10.5061/dryad.8sf7m0cw8](https://doi.org/10.5061/dryad.8sf7m0cw8)
    		* no duplicate orthogroup composition metrics; to use to plot and summarize similarity indices

