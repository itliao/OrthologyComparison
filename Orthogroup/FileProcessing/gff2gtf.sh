echo "Aethionema"
gffread \
	/path_to_gff/Ae.arabicum_v3.1_annotations_utr.gff  \
	-T -o /path_to_output/Aarabicum.gtf 

echo "Arabidopsis"
gffread \
	/path_to_gff/Athaliana_447_Araport11.gene_exons.gff3 \
	-T -o /path_to_output/Athaliana.gtf
	
echo "Brassica"
gffread \
	/path_to_gff/BrapaFPsc_277_v1.3.gene_exons.gff3 \
	-T -o /path_to_output/Brapa.gtf

echo "Camelina"
gffread \
	/path_to_gff/Camelina_sativa.Cs.55.chr.gff3  \
	-T -o /path_to_output/Csativa.gtf 

echo "Capsella"	
gffread \
	/path_to_gff/Crubella_474_v1.1.gene_exons.gff3 \
	-T -o /path_to_output/Crubella.gtf  
 
echo "Cardamine"	
gffread \
	/path_to_gff/carhr38.gff \
	-T -o /path_to_output/Chirsuta.gtf  

 echo "Sinapis"	
gffread \
	/path_to_gff/Sal.Chr.20210627.gff \
	-T -o /path_to_output/Salba.gtf 
