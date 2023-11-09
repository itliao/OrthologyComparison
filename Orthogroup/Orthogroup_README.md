# genome and annotation files 

The genome files used for each the of the species are:

- *Aethionema arabicum* (Aar)
  - Origin: Fernandez-Pozo et al., 2021
  - Version: 3.1
  - gff with UTR, CDS, and exons: Ae.arabicum_v3.1_annotations_utr.gff
  - protein FASTA, primary transcripts: Ae.arabicum_v3.1_proteins_noTE_23160.fasta 

- *Arabidopsis thaliana* (Ath)
  - Origin: Phytozome
  - Version: Araport11
  - gff with UTR, CDS, and exons: Athaliana_447_Araport11.gene_exons.gff3
  - protein FASTA, primary transcripts: Athaliana_447_Araport11.protein_primaryTranscriptOnly.fa

- *Brassica rapa* (Bra)
  - Origin: Phytozome
  - Version: v1.3
  - gff with UTR, CDS, and exons: BrapaFPsc_277_v1.3.gene_exons.gff3
  - protein FASTA, primary transcripts: BrapaFPsc_277_v1.3.protein_primaryTranscriptOnly.fa

- *Cardamine hirsuta* (Chi)
  - Origin: Gan et al., 2016
  - Version: v1.0
  - gff with UTR, CDS, and exons: carhr38.gff
  - protein FASTA, primary transcripts: carhr38.aa.fa  

- *Capsella rubella* (Cru)
  - Origin: Phytozome
  - Version: v1.1
  - gff with UTR, CDS, and exons: Crubella_474_v1.1.gene_exons.gff3
  - protein FASTA, primary transcripts: Crubella_474_v1.1.protein_primaryTranscriptOnly.fa

- *Camelina sativa* (Csa)
  - Origin: EnsemblPlants
  - Version: 55
  - gff with UTR, CDS, and exons: Camelina_sativa.Cs.55.chr.gff3
  - protein FASTA, primary transcripts: Camelina_sativa.Cs.pep.all.fa.gz 

- *Sinapis alba* (Sal)
  - Origin: Yang et al. 2022
  - Version: v1.0
  - gff with UTR, CDS, and exons: Sal.Chr.20210627.gff
  - protein FASTA, primary transcripts: Sal.pep

- *Thlaspi arvense* (Tar)
  - Origin: NCBI, directory GCA_911865555.2
  - Version: v2
  - gff with UTR, CDS, and exons: genomic.gff
  - protein FASTA, primary transcripts: protein.faa 

## file preprocessing (mostly for OrthNet) 

Only primary transcripts should be included. While most of the protein files from Phytozome should only have primary transcripts, resources from other genomic repositories may not be formatted this way. 

To make sure all the input files and the gene ID names are the same across all 4 orthology programs, most genome annotation gff3/gtf files and protein FASTA files were modified in the following scripts:

**STEP 1: Convert all gff to gtf**
- script: gff2gtf.sh using gffread v0.12.7 

**STEP 2: Extract primary transcripts from original gff3**
- generate list of primary transcripts
- script: ExtractPrimaryGeneID_230131.ipynb

**STEP 3: parse the gtf (for OrthNet)**
- script: parse_gtf.sh

**STEP 4: extract primary transcripts from FASTA**
- only ran if the protein fasta was not originally from phytozome or if it already contained only primary transcripts
- script: ExtractPrimaryTranscriptProteinFasta.py
- species: Cardamine hirsuta, Camelina sativa

**STEP 5: convert IDs so that they are exactly the same (primarily for OrthNet)**
- if gtf and FASTA files have completely different gene and protein IDs
- species: Thlaspi arvense
- script: Convert_Tarvense_v2_files_221214.ipynb
- inputs:
  - protein.faa
  - genomic.gff
  - Tar.gtfParsed.txt
- outputs:
  - Tarvense_v2_geneProtein_pair.tsv
  - Tar_v2_protein_geneID.fasta
  - Tar_v2_proteinID.gtfParsed.txt (used for OrthNet input)
- For all species, the names were slightly modified using sed:
	- sed_file_commands.txt

## Final files as inputs for all programs:

* *Aethionema arabicum*(Aar)
  * gtf with UTR, CDS, and exons: Aar.gtfParsed.mod.txt
  * protein FASTA, primary transcripts: Ae.arabicum_v3.1_proteins_noTE_23160.fasta

* *Arabidopsis thaliana* (Ath)
  * gtf with UTR, CDS, and exons: Ath.gtfParsed.mod.txt
  * protein FASTA, primary transcripts: Athaliana_447_Araport11.protein_primaryTranscriptOnly.fa

* *Brassica rapa* (Bra)
  * gtf with UTR, CDS, and exons: Bra.gtfParsed.mod.txt
  * protein FASTA, primary transcripts: BrapaFPsc_277_v1.3.protein_primaryTranscriptOnly.fa

* *Cardamine hirsuta* (Chi)
  * gtf with UTR, CDS, and exons: Chi.gtfParsed.mod.txt
  * protein FASTA, primary transcripts: Chirsuta.primary.pep.fa

* *Capsella rubella* (Cru)
  * gtf with UTR, CDS, and exons: Cru.gtfParsed.mod.txt
  * protein FASTA, primary transcripts: Crubella_474_v1.1.protein_primaryTranscriptOnly.fa

* *Camelina sativa* (Csa)
  * gtf with UTR, CDS, and exons: Csa.gtfParsed.mod.txt
  * protein FASTA, primary transcripts: Csativa.primary.pep.mod.fa

* *Sinapis alba* (Sal)
  * gtf with UTR, CDS, and exons: Sal.gtfParsed.mod.txt
  * protein FASTA, primary transcripts: Sal.pep

* *Thlaspi arvense* (Tar)
  * gtf with UTR, CDS, and exons: Tar.gtfParsed.mod.txt
  * protein FASTA, primary transcripts: Tarvense.primary.pep.fa (from Tar_v2_protein_geneID.fasta)

# scripts to run orthology inference algorithms
## Orthofinder 
- filename extensions: *.fa
- scripts: 
	- orthofinder_all.sh
- Note: create separate directories for each group of proteomes and for each type of analyses

## SonicParanoid 
- filename extensions: *.fa
- scripts:
	- sonicparanoid.sh
- Note: create separate directories for each group of proteomes and for each type of analyses

## Broccoli 
- filename extensions: *.fasta
- scripts:
	- broccoli.sh
- Note: create separate directories for each group of proteomes and for each type of analyses

## OrthNet 

**STEP 1: Prepare input files**
- Create a list of species: date_test.list
```
	echo 'Ath Chi Cru Tar Aar' | tr ' ' '\n' > 230927_diploid.list
	echo 'Ath Chi Cru Tar Bra Csa Sal Aar' | tr ' ' '\n' > 230214_fullset.list
````
- Input 1: obtain reduced gtf files
  - parse gtf with parse_gtf_2table.py
  - script:
  	- parse_gtf.sh
  - inputs:
  	- see above. File extensions: speciesName.gtf
  - outputs: 
  	- File extensions: XXX.gtfParsed.mod.txt
	
- Input 2: sequence search within species (paralogs)
  - script:
  	- mmseqs2.sh 
  - inputs: 
  	- protein filename extensions: XXX.pep.rep

- Input 3: sequence search between species
1) create a script to run all possible pairwise combinations of blast 
  - script: 
  	- createPairwiseBlast.sh
  - output:
  	- 230214_full_pairwiseMMseqs2.sh
  	- 230927_diploid_pairwiseMMseqs2.sh
  	
2) Modify output from createPairwiseBlast.sh with correct paths
  - modified script:
	  - mmseqs2_pairwise_search.sh 

3) run mmseqs2_pairwise_search.sh
	
4) combine output into a folder (run these lines in directly in the terminal)
- output directories:
	- /u/scratch/i/irenelia/orthology/OrthNet/230214_pairs/
	- /u/scratch/i/irenelia/orthology/OrthNet/230927_pairs/

  ```
	for f in out__*.txt; do f2=${f##*out__}; cut -f1,2 $f | uniq > BestHits__${f2%%.txt}.list; done
	```
	Then:
	```
	mkdir ./BHPairs; mv BestHits__*.list ./BHPairs
	```
**STEP 2: CLFinder**

1) For annotating tandem duplication information
  - script: 
  	- CLfinder_combine1.sh
  - input:
  	- XXX.gtfParsed.mod.txt
  	- XXX.PG
  - output:
  	- XXX.gtfParsed.PG.txt
  	- XXX.gtfParsed.TD.txt

2) run CLfinder
  - script:
  	- CLfinder_2.sh
  - input:
  	- files in BHPairs
  	- XXX.gtfParsed.TD.txt
  - output:
  	- files in 230214_full_out
  		- 230214_full.4OrthNet.input 
  	- files in 230927_full_out
  		- 230927_full.4OrthNet.input 

3) update best hit pairs
  - script:
  	- CLfinder_update3.sh
  - input:
  	- 230214_full.4OrthNet.input (230214_full_out)
  	- 230927_full.4OrthNet.input (230927_full_out)
  - output:
  	- new, updated files in BHPairs.1
	
4) run CLfinder again with updated best hit pairs
  - script:
  	- CLfinder_4.sh
  - input:
  	- updated files in BHPairs.1
  	- XXX.gtfParsed.TD.txt
  - output:
  	- files in 230214_full_out.1 (230914_full_out.1)
  		- 230214_full.4OrthNet.input
  	- files in 230927_diploid_out.1 (230927_diploid_out.1)
  		- 230927_diploid.4OrthNet.input

5) summary report for all pairwise CLfinder analyses
  - script:
  	- CLfinder_summary5.sh
  - input:
  	- files in 230214_full_out.1
  	- files in 230927_diploid_out.1
  - output:
  	- 230214_CLfinder_summary.txt
  	- 230927_CLfinder_summary.txt
	
**STEP 3: OrthNet**

1) Create initial hard clusters
  - script: 
  	- ONfinder_1.sh
  - input:
  	- 230214_full.4OrthNet.input 
  	- 230927_diploid.4OrthNet.input
  	- XXX.gtfParsed.TD.txt
  - output:
  	- 230214_full.4OrthNet.input.TDadded
  	- 230214_full.4OrthNet.input_TDid.list 
  	- files in 230214/ONfinder
  	
  	- 230927_full.4OrthNet.input.TDadded
  	- 230927_full.4OrthNet.input_TDid.list 
  	- files in 230927/ONfinder
	
2) Markov clustering (mcl) of hard clusters
  - NOTE: changed the inflation parameter from 1.2 (default) to 1.5 to match OrthoFinder and SonicParanoid
  - script:
  	- ONfinder_2.sh
  - input:
  	- files in ONfinder
  - output:
  	- files in mcl
  		- 230214_full_TD1.5_rC1.2_rNC0.5_uC0.6_uNC0.25_I1.5_mclOutPC.txt
  		- 230927_diploid_TD1.5_rC1.2_rNC0.5_uC0.6_uNC0.25_I1.5_mclOutPC.txt
		
3) Update best-hit pairs and OrthNets after mcl
  - script:
  	- ONfinder_3.sh
  - input:
  	- 230214_full_TD1.5_rC1.2_rNC0.5_uC0.6_uNC0.25_I1.5_mclOutPC.txt
  	- files in 230214/ONfinder
  	
  	- 230927_diploid_TD1.5_rC1.2_rNC0.5_uC0.6_uNC0.25_I1.5_mclOutPC.txt
  	- files in 230927/ONfinder
  - output:
  	- updated files in 230214/BHPairs.2
  	- updated files in 230214_full.2
  		- 230214_full.clstrd.afterMCL.nodes.mclOutput (for downstream analyses)
  		
  	- updated files in 230927/BHPairs.2
  	- updated files in 230927_diploid.2
  		- 230927_full.clstrd.afterMCL.nodes.mclOutput (for downstream analyses)

4) create file for reading into Cytoscape
  - script:
  	- ONfinder_4.sh
  - input:
  	- 230214_full.clstrd.afterMCL.edges 
  	- 230927_diploid.clstrd.afterMCL.edges 
  - output:
  	- 230214_full.clstrd.afterMCL.edges.sif
  	- 230927_diploid.clstrd.afterMCL.edges.sif
	
5) re-create the CLfinder summary report
  - script:
  	- ONfinder_5.sh
  - input:
  	- files from 230214_full.2
  	- files from 230927_diploid.2
  - output:
  	- 230215_CLfinder_summary.afterOrthNet.txt
  	- 230927_CLfinder_summary.afterOrthNet.txt

**formating files for comparisons** 
for OrthNet
  - script:
  	- OrthNet_to_Spreadsheet_230915.ipynb
  	- OrthNet_to_Spreadsheet_230927.ipynb
  - inputs:
  	- 230214_full.2_1.5/230214_full.clstrd.afterMCL.nodes.mclOutput
  	- 230927_diploid.2/230927_diploid.clstrd.afterMCL.nodes.mclOutput
  - output:
  	- OrthNet_groups_full_230915.txt
  	- OrthNet_groups_diploid_230927.txt
- Note: output format similar to that of orthofinder outputs (comma and space placement)

# output files - to use for inputs for downstream analyses 

- dip - refers to the "diploid set" of species
- full - refers to the "higher ploidy set" of species

- OF_b - OrthoFinder-BLAST
  - dip_OF_b-N0.tsv
  - full_OFb-N0.tsv
- OF_d - OrthoFinder-DIAMOND
  - dip_OF_d-N0.tsv
  - full_OFd-N0.tsv
- OF_m - OrthoFinder-MMseqs2
  - dip_OF_m-N0.tsv
  - full_OFm-N0.tsv
- SP_d - SonicParanoid-DIAMOND
  - dip_SP_d-flat.ortholog_groups.tsv
  - full_SPd-flat.ortholog_groups.tsv
- SP_m - SonicParanoid-MMseqs2
  - dip_SP_m-flat.ortholog_groups.tsv
  - full_SPm-flat.ortholog_groups.tsv
- BR - Broccoli
  - dip_BR-table_OGs_protein_names.txt
  - full_BR-table_OGs_protein_names.txt
- ON - OrthNet
  - OrthNet_groups_diploid_230915.txt
  - OrthNet_groups_full_230915.txt
