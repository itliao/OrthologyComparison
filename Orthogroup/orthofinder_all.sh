echo "Running orthofinder - brassicaceae - diploid, blast"
/path_to_orthofinder/orthofinder.py \
	-f /path_to_proteomes/dip_b/ \
	-S blast -M msa 

echo "Running orthofinder - brassicaceae - diploid, diamond"
/path_to_orthofinder/orthofinder.py \
	-f /path_to_proteomes/dip_d/ \
	-S diamond -M msa 
		
echo "Running orthofinder - brassicaceae - diploid, mmseqs"
/path_to_orthofinder/orthofinder.py \
	-f /path_to_proteomes/dip_m/ \
	-S mmseqs -M msa 
	
echo "Running orthofinder - brassicaceae - full, blast"
/path_to_orthofinder/orthofinder.py \
	-f /path_to_proteomes/full_b/ \
	-S blast -M msa 

echo "Running orthofinder - brassicaceae - full, diamond"
/path_to_orthofinder/orthofinder.py \
	-f /path_to_proteomes/full_d/ \
	-S diamond -M msa 
		
echo "Running orthofinder - brassicaceae - full, mmseqs"
/path_to_orthofinder/orthofinder.py \
	-f /path_to_proteomes/full_m/ \
	-S mmseqs -M msa 
