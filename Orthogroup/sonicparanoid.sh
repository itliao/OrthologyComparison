echo "sonicparanoid, brassicaceae - full, diamond"
sonicparanoid \
	-i /path_to_proteomes/full_SPd \
	-o /path_to_proteomes/sonicparanoid/full_SPd_out \
	--aln-tool diamond \
	-p 230209_d
	
echo "sonicparanoid, brassicaceae - full, mmseqs"
sonicparanoid \
	-i /path_to_proteomes/full_SPm \
	-o /path_to_proteomes/sonicparanoid/full_SPm_out \
	--aln-tool mmseqs \
	-p 230209_m

echo "sonicparanoid, brassicaceae - diploid, diamond"
sonicparanoid \
	-i /path_to_proteomes/dipSP_d \
	-o /path_to_proteomes/sonicparanoid/dipD_out \
	--aln-tool diamond \
	-p 230922_d
	
echo "sonicparanoid, brassicaceae - diploid, mmseqs"
sonicparanoid \
	-i /path_to_proteomes/dipSP_m \
	-o /path_to_proteomes/sonicparanoid/dipM_out \
	--aln-tool mmseqs \
	-p 230922_m
