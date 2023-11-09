# parse_gtf_2table.py is a python script from the OrthNet workflow for converting the information from an annotation into a tabular format

echo "Aethionema"
/path_to_OrthNet_scripts/parse_gtf_2table.py \
    -pr /path_to_file/Aarabicum.gtf \
    -e /path_to_file/Aar_primaryGeneID.txt \
    /path_to_file/Aar.filterGTF.txt > \
    /path_to_file/Aar.filterGTF.log

/path_to_OrthNet_scripts/parse_gtf_2table.py \
    -pr /path_to_file/Aar.filterGTF.txt \
    /path_to_file/Aar.gtfParsed.txt > \
    /path_to_file/Aar.gtfParsed.log
        
echo "Arabidopsis"
/path_to_OrthNet_scripts/parse_gtf_2table.py \
    -pr /path_to_file/Athaliana.gtf \
    -e /path_to_file/Ath_primaryGeneID.txt \
    /path_to_file/Ath.filterGTF.txt > \
    /path_to_file/Ath.filterGTF.log

/path_to_OrthNet_scripts/parse_gtf_2table.py \
    -pr /path_to_file/Ath.filterGTF.txt \
    /path_to_file/Ath.gtfParsed.txt > \
    /path_to_file/Ath.gtfParsed.log
        
echo "Brassica"
/path_to_OrthNet_scripts/parse_gtf_2table.py \
	-pr /path_to_file/Brapa.gtf \
	-e /path_to_file/Bra_primaryGeneID.txt \
	/path_to_file/Bra.filterGTF.txt > \
	/path_to_file/Bra.filterGTF.log

/path_to_OrthNet_scripts/parse_gtf_2table.py \
	-pr /path_to_file/Bra.filterGTF.txt \
	/path_to_file/Bra.gtfParsed.txt > \
	/path_to_file/Bra.gtfParsed.log

echo "Camelina"
/path_to_OrthNet_scripts/parse_gtf_2table.py \
	-pr /path_to_file/Csativa.gtf \
	-e /path_to_file/Csa_primaryGeneID.txt \
	/path_to_file/Csa.filterGTF.txt > \
	/path_to_file/Csa.filterGTF.log

/path_to_OrthNet_scripts/parse_gtf_2table.py \
	-pr /path_to_file/Csa.filterGTF.txt  \
	/path_to_file/Csa.gtfParsed.txt > \
	/path_to_file/Csa.gtfParsed.log
	
echo "Capsella" 
/path_to_OrthNet_scripts/parse_gtf_2table.py \
    -pr /path_to_file/Crubella.gtf \
    -e /path_to_file/Cru_primaryGeneID.txt \
    /path_to_file/Cru.filterGTF.txt > \
    /path_to_file/Cru.filterGTF.log

/path_to_OrthNet_scripts/parse_gtf_2table.py \
    -pr /path_to_file/Cru.filterGTF.txt \
    /path_to_file/Cru.gtfParsed.txt > \
    /path_to_file/Cru.gtfParsed.log

echo "Cardamine"	
/path_to_OrthNet_scripts/parse_gtf_2table.py \
	-pr /path_to_file/Chirsuta.gtf \
	-e /path_to_file/Chi_primaryGeneID.txt \
	/path_to_file/Chi.filterGTF.txt > \
	/path_to_file/Chi.filterGTF.log

/path_to_OrthNet_scripts/parse_gtf_2table.py \
	-pr /path_to_file/Chi.filterGTF.txt \
	/path_to_file/Chi.gtfParsed.txt > \
	/path_to_file/Chi.gtfParsed.log

echo "Sinapis"	
/path_to_OrthNet_scripts/parse_gtf_2table.py \
	-pr /path_to_file/Salba.gtf \
	-e /path_to_file/Sal_primaryGeneID.txt \
	/path_to_file/Sal.filterGTF.txt > \
	/path_to_file/Sal.filterGTF.log
	
/path_to_OrthNet_scripts/parse_gtf_2table.py \
	-pr /path_to_file/Sal.filterGTF.txt \
	/path_to_file/Sal.gtfParsed.txt > \
	/path_to_file/Sal.gtfParsed.log

echo "Thlaspi"
/path_to_OrthNet_scripts/parse_gtf_2table.py \
	-pr /path_to_file/Tarvense.gtf \
	-e /path_to_file/Tar_primaryGeneID.txt \
	/path_to_file/Tar.filterGTF.txt > \
	/path_to_file/Tar.filterGTF.log
	
/path_to_OrthNet_scripts/parse_gtf_2table.py \
    -pr /path_to_file/Tarvense.gtf  \
    /path_to_file/Tar.gtfParsed.txt > \
    /path_to_file/Tar.gtfParsed.log
