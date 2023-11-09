###################
# for diploid set #
###################

mkdir /path_to_tmp/tmp # temporary folder for MMseqs2 runs

while read g

do 

mmseqs createdb /path_to_files/${g}.pep.rep \
	/path_to_files/${g}_DB
	
mmseqs createindex /path_to_files/${g}_DB \
	/path_to_tmp/tmp
	
mmseqs cluster /path_to_files/${g}_DB \
	/path_to_files/${g}_c \
	/path_to_tmp/tmp --max-seqs 50000 -c 0.5
	
mmseqs createtsv /path_to_files/${g}_DB \
	/path_to_files/${g}_DB \
	/path_to_files/${g}_c \
	/path_to_files/${g}_c.tsv

/path_to_OrthNet_scripts/parse_mmseqs_clusters.py \
	-H /path_to_files/${g}_c.tsv \
	/path_to_files/${g}.PG

done < /path_to_files/230927_diploid.list

#########################
# for higher ploidy set #
#########################

echo "running MMseqs2"

mkdir /path_to_tmp/tmp # temporary folder for MMseqs2 runs

while read g

do 

mmseqs createdb /path_to_files/${g}.pep.rep \
	/path_to_files/${g}_DB
	
mmseqs createindex /path_to_files/${g}_DB \
	/path_to_tmp/tmp
	
mmseqs cluster /path_to_files/${g}_DB \
	/path_to_files/${g}_c \
	/path_to_tmp/tmp --max-seqs 50000 -c 0.5
	
mmseqs createtsv /path_to_files/${g}_DB \
	/path_to_files/${g}_DB \
	/path_to_files/${g}_c \
	/path_to_files/${g}_c.tsv

/path_to_OrthNet_scripts/parse_mmseqs_clusters.py \
	-H /path_to_files/${g}_c.tsv \
	/path_to_files/${g}.PG

done < /path_to_files/230214_fullset.list

