# summary plots of orthogroup composition 
## STEP 1: Convert OG into table of numbers to input for visualization

* script
	* ConvertOG_for_Visualization_FINAL.ipynb

* inputs		
	* orthogroups from diploid set and higher ploidy set
		* [see Orthogroup_README.md](https://github.com/itliao/OrthologyComparison/blob/main/Orthogroup/Orthogroup_README.md#output-files---to-use-for-inputs-for-downstream-analyses)
	
* outputs
	* algorithm + "_allCounts.tsv
	* algorithm + "_binary.tsv"
	* algorithm + "_tenMax.tsv"
		* where "algorithm" stands for orthofinder, sonicparanoid, broccoli, or orthnet and each variation on those algorithms
	* [diploid_set](/Summary_Statistics/diploidBrass.zip)
	* [higher ploidy set](/Summary_Statistics/fullBrass.zip)

## STEP 2: Make plots - histogram, stacked bar, heatmap, and Upset bar plot

* Figures/Tables:
	* Figure 2 (part): UpSet bar plots, number of OG with specific species composition
	* Figure S1: Distribution of the number of species per OG
	* Figure S5: Percentage of the number of genes per species in OG , Heatmap of the number of genes per species in OG 

* script:
	* orthologyPlots_counts_FINAL.R
	
* inputs - files from above
	* algorithm + "_allCounts.tsv
	* algorithm + "_binary.tsv"
	* algorithm + "_tenMax.tsv"
		* where "algorithm" stands for orthofinder, sonicparanoid, broccoli, or orthnet and each variation on those algorithms
	* [diploid_set](/Summary_Statistics/diploidBrass.zip)
	* [higher ploidy set](/Summary_Statistics/fullBrass.zip)
	
* plot output
	* [diploid_set](/Summary_Statistics/diploidBrass_Plots.zip)
	* [higher ploidy set](/Summary_Statistics/fullBrass_Plots.zip)
	* plot types 
		* Distribution of the number of species per OG
		* Distribution of the number of OG with specific species composition (upset)
		* Percentage of the number of genes per species in OG 
		* Heatmap of the number of genes per species in OG 
		
# Reformatting plots and figures, statistics on comparisons
* Figures:
	* Figure 2 (part): Stacked bars of the number of genes in an orthogroup across algorithms
	* Figure 3, S4: Distribution of the number of genes per species **grouped by species**
	* Figure 4, S6: Orthogroup gene composition comparisons 
	* Figure S2: Violin plot with average number of species in an orthogroup across algorithms
	* Figure S7: Distribution of pairwise comparisons between orthology inference algorithms
	* Figure S8: Heatmap of metrics comparison of orthology inference for species pairs against baseline method
	* Figure S9: Stacked bars of predicted orthology relationships between all species pairs across algorithms
* Tables:
	* Table 3: Comparison of the number of species per orthogroup detected across orthology interference algorithms
	* Table S1: Summary statistics of the number of species in an orthogroup across algorithms
	* Table S5: Summary statistics for metrics calculated for all-against-all comparisons of orthogroup compositions among all algorithms

* script
	* orthology_stats_FINAL.R
		
* inputs
	* files from above
		* algorithm + "_allCounts.tsv
		* algorithm + "_binary.tsv"
		* algorithm + "_tenMax.tsv"
			* where "algorithm" stands for orthofinder, sonicparanoid, broccoli, or orthnet and each variation on those algorithms
		* [diploid_set](/Summary_Statistics/diploidBrass.zip)
		* [higher ploidy set](/Summary_Statistics/fullBrass.zip)
	* [combinedOrthologCount_reformat_230928.csv](/Gene_Composition_Comparison_Species_Pairs/combinedOrthologCount_reformat_230928.csv)
	
	* diploid_pairwiseMetrics_231003.csv **found in DRYAD**
	* full_pairwiseMetrics_231004.csv **found in DRYAD**
	* [diploid_speciesPairsJI_230929.csv](/Gene_Composition_Comparison_Species_Pairs/diploid_speciesPairsJI_230929.csv)
	* [mixed_speciesPairsJI_230929.csv](/Gene_Composition_Comparison_Species_Pairs/mixed_speciesPairsJI_230929.csv)
	
* outputs:
	* [diploid_set](/Summary_Statistics/diploid_outputs)
		* diploid_numSpPerOG_230928.pdf (Fig. 2)
		* diploid_numSpPerOG_sumstats_230928.csv (Table S1)
		* diploid_numSpPerOG_violin_230928.pdf (Fig. S2)
		* diploid_numSpPerOG_kwTest_230928.txt (Table 3)
		* diploid_numGenesPerSpOG_230928.pdf (Fig. 3, S4)
		* diploid_numGenesPerSpOG_sumStats_230928.csv
			* distribution of number of genes per species merged together
		* diploid_numGenesPerSpOG_wilcox_230928.csv (Table S3)
		* diploid_numGenesPerSpOG_kwTest_230928.txt (Table S3)
		* orthologyCountsPlots_230928.pdf (Fig. S9) 
		* metricsSummary_DipMix_231004.xlsx (Table S5)
		* metricsSummaryDiploid_231004.pdf (Fig. 4, S6, S7)
		* heatmap_SpeciesPairs_Diploid_231003.pdf (Fig. S8)
		
		
	* [higher_ploidy_set](/Summary_Statistics/higher_ploidy_outputs)
		* full_numSpPerOG_230928.pdf (Fig. 2)
		* full_numSpPerOG_sumstats_230928.csv (Table S1)
		* full_numSpPerOG_violin_230928.pdf (Fig. S2)
		* full_numSpPerOG_kwTest_230928.txt (Table 3)
		* full_numGenesPerSpOG_230928.pdf (Fig. 3, S4)
		* full_numGenesPerSpOG_sumStats_230928.csv 
			* distribution of number of genes per species merged together
		* full_numGenesPerSpOG_wilcox_230928.csv (Table S3)
		* full_numGenesPerSpOG_kwTest_230928.txt (Table S3)
		* orthologyCountsPlots_230928.pdf (Fig. S9) 
		* metricsSummary_DipMix_231004.xlsx (Table S5)
		* metricsSummaryMixed_231004.pdf (Fig. 4, S6, S7)
		* heatmap_SpeciesPairs_Mixed_231003.pdf (Fig. S8)
		

