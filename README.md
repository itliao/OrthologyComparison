# OrthologyComparison
Customizable orthology inference algorithms generate similar predicted orthogroups among Brassicaceae species

Maintainer's contact info: 
Irene T. Liao - ireneliao@ucla.edu


# orthology inference algorithms

Orthofinder v2.5.2 
* note - the latest version, v2.5.4, has an issue for using MMseqs2, and thus, v2.5.2 was used instead
	https://github.com/davidemms/OrthoFinder
	
SonicParanoid v1.3.8 
	http://iwasakilab.k.u-tokyo.ac.jp/sonicparanoid/
	
Broccoli v1.2
	https://github.com/rderelle/Broccoli
	
OrthNet 
	https://github.com/ohdongha/OrthNet

# General workflow & directories for scripts and files

**Orthogroup directory** 
  - scripts, inputs, and outputs (orthogroups) from running orthology inference algorithms

**Summary_Statistics directory**
  - scripts and summary results describing orthogroup results

**Gene_Composition_Comparison_Orthogroups directory**
  - scripts and results for comparing the gene compositon across orthogroups from different orthology inference algorithms
  - files for conservative and inclusive orthogroups of Arabidopsis genes after comparison across algorithms

**Gene_Composition_Comparison_Species_Pairs directory** 
  - scripts and results for comparing the gene compositon across orthogroups from different orthology inference algorithms

**YABBY_Case_Study directory** 
  - scripts and results for comparing YABBY orthogroups across different orthology inference algorithms
  





