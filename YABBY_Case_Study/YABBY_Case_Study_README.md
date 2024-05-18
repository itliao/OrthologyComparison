# YABBY case study

## STEP 1: Extract protein sequences from list of genes from inclusive orthogroup file 
* script: 
	* ExtractProteinFastaFromSpecies_FINAL.ipynb

* input 
	* protein fasta files from all species
		* [see Orthogroup_README.md for protein fasta sources](/Orthogroup/Orthogroup_README.md)
	* [inclusiveOG_fullSet_230928.csv](/Gene_Composition_Comparison_Orthogroups/higher_ploidy_output/inclusiveOG_fullSet_230928.csv)
	* AT_YABBY_geneID.txt

* [output](/YABBY_Case_Study/sequences): 
	* AT1G08465.1_sequence.fa
	* AT1G23420.2_sequence_mod.fa 
		* manually deleted an extra copy of gene Sal07g28570L
	* AT1G69180.1_sequence.fa
	* AT2G26580.1_sequence.fa
	* AT2G45190.1_sequence.fa
	* AT4G00180.1_sequence.fa
	* ALL_sequence.fa 
		* combined version to run through MAFFT, added Picea YABBY sequence as an outgroup gene copy to root the tree

## STEP 2: Align sequences using MAFFT with default settings 
MAFFT: https://www.ebi.ac.uk/Tools/msa/mafft/

* input
	* files above

* [output (.aln files)](/YABBY_Case_Study/alignments)
	* ALL_YABBY_231017.aln
	* AT1G08465.1.aln
	* AT1G23420.2.aln
	* AT1G69180.1.aln
	* AT2G26580.1.aln
	* AT2G45190.1.aln
	* AT4G00180.1.aln
	
* Visualize MAFFT alignments with Aliview
Aliview: https://ormbunkar.se/aliview/
  * [screenshots of sequences (Fig. S10)](/YABBY_Case_Study/screenshot) 

## STEP 3: Build trees using RAxML

CIPRES: https://www.phylo.org/
* parameters changed
	* bootstraps: 1000
	* defined outgroups
		* All sequences: Picea_glaucaBT115385
		* AT1G08465.1: Aa31sc703G10
		* AT1G23420.2: Aa31LG7G7450
		* AT1G69180.1: Aa31LG4G11850
		* AT2G26580.1: Aa31LG5G9130
		* AT2G45190.1: Aa31LG1G26740
		* AT4G00180.1: none (no Aethionema sequence)

## STEP 4: Visualize and modify trees with FigTree
[RAxML outputs and trees](/YABBY_Case_Study/RAxML_out)
* open RAxML_bipartitions.brass file and make slight visual modifications
* resulting trees: **sequence_id**_tree.pdf

## STEP 5: Networks from OrthNet

* script to get the clusters:
	* Modify_sif_for_cytoscape_240222_cleaned.ipynb
		
* inputs:
	* OrthNet_YAB.txt
 	* 230927_diploid.clstrd.afterMCL.edges.sif
  	* 230214_full.clstrd.afterMCL.edges.sif

* outputs:
 	* diploid
  		* [sif files and pdfs](/YABBY_Case_Study/OrthNetCytoscape/diploid) 
  	* diploid+higher ploidy
  		* [sif files and pdfs](/YABBY_Case_Study/OrthNetCytoscape/higher_ploidy) 
   

