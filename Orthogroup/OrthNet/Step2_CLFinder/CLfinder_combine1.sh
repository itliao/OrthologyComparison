###############
# diploid set #
###############

echo "CLfinder Step 1: combine input 1 and input 2 and tandem duplication info"

while read g

do 

/path_to_OrthNet_scripts/join_files_by_NthCol.py \
    /path_to_files/${g}.gtfParsed.mod.txt 1 1 \
    /path_to_files/${g}.PG \
    /path_to_files/${g}.gtfParsed.PG.txt

/path_to_OrthNet_scripts/TD_finder.py \
    /path_to_files/${g}.gtfParsed.PG.txt \
    $g 4 /path_to_files/${g}.gtfParsed.TD.txt

done < /path_to_files/230927_diploid.list

#####################
# higher ploidy set #
#####################

echo "CLfinder Step 1: combine input 1 and input 2 and tandem duplication info"

while read g

do 

/path_to_OrthNet_scripts/join_files_by_NthCol.py \
    /path_to_files/${g}.gtfParsed.mod.txt 1 1 \
    /path_to_files/${g}.PG \
    /path_to_files/${g}.gtfParsed.PG.txt

/path_to_OrthNet_scripts/TD_finder.py \
    /path_to_files/${g}.gtfParsed.PG.txt \
    $g 4 /path_to_files/${g}.gtfParsed.TD.txt

done < /path_to_files/230214_full.list
