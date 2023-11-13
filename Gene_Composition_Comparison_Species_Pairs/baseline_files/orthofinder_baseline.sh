#############################################
#scripts for pairwise species - diploid	    #
# testing BLAST + mcl output				#
# use Orthofinder v2.5.2 (download source)  #
#############################################

# 230118
echo "Running orthofinder - AthChi, blast"
/path_to_orthofinder/orthofinder.py \
	-f /path_to_proteome/AthChi/ \
	-S blast -og

echo "Running orthofinder - AthCru, blast"
/path_to_orthofinder/orthofinder.py \
	-f /path_to_proteome/AthCru/ \
	-S blast -og
	
echo "Running orthofinder - AthTar, blast"
/path_to_orthofinder/orthofinder.py \
	-f /path_to_proteome/AthTar/ \
	-S blast -og

echo "Running orthofinder - CruChi, blast"
/path_to_orthofinder/orthofinder.py \
	-f /path_to_proteome/CruChi/ \
	-S blast -og
	
echo "Running orthofinder - CruTar, blast"
/path_to_orthofinder/orthofinder.py \
	-f /path_to_proteome/CruTar/ \
	-S blast -og
	
echo "Running orthofinder - ChiTar, blast"
/path_to_orthofinder/orthofinder.py \
	-f /path_to_proteome/ChiTar/ \
	-S blast -og

########################################################
# scripts for pairwise species - remainder of full set #
# testing BLAST + mcl output						   #
# use Orthofinder v2.5.2 (download source)             #
########################################################

# 230209
echo "Running orthofinder - AthSal, blast"
/path_to_orthofinder/orthofinder.py \
	-f /path_to_proteome/AthSal/ \
	-S blast -og

echo "Running orthofinder - AthBra, blast"
/path_to_orthofinder/orthofinder.py \
	-f /path_to_proteome/AthBra/ \
	-S blast -og
	
echo "Running orthofinder - AthCsa, blast"
/path_to_orthofinder/orthofinder.py \
	-f /path_to_proteome/AthCsa/ \
	-S blast -og
	
echo "Running orthofinder - AthAar, blast"
/path_to_orthofinder/orthofinder.py \
	-f /path_to_proteome/AthAar/ \
	-S blast -og

echo "Running orthofinder - CruSal, blast"
/path_to_orthofinder/orthofinder.py \
	-f /path_to_proteome/CruSal/ \
	-S blast -og
	
echo "Running orthofinder - CruBra, blast"
/path_to_orthofinder/orthofinder.py \
	-f /path_to_proteome/CruBra/ \
	-S blast -og

echo "Running orthofinder - CruCsa, blast"
/path_to_orthofinder/orthofinder.py \
	-f /path_to_proteome/CruCsa/ \
	-S blast -og
	
echo "Running orthofinder - CruAar, blast"
/path_to_orthofinder/orthofinder.py \
	-f /path_to_proteome/CruAar/ \
	-S blast -og
	
echo "Running orthofinder - ChiSal, blast"
/path_to_orthofinder/orthofinder.py \
	-f /path_to_proteome/ChiSal/ \
	-S blast -og

echo "Running orthofinder - ChiBra, blast"
/path_to_orthofinder/orthofinder.py \
	-f /path_to_proteome/ChiBra/ \
	-S blast -og

echo "Running orthofinder - ChiCsa, blast"
/path_to_orthofinder/orthofinder.py \
	-f /path_to_proteome/ChiCsa/ \
	-S blast -og

echo "Running orthofinder - ChiAar, blast"
/path_to_orthofinder/orthofinder.py \
	-f /path_to_proteome/ChiAar/ \
	-S blast -og
	
echo "Running orthofinder - TarSal, blast"
/path_to_orthofinder/orthofinder.py \
	-f /path_to_proteome/TarSal/ \
	-S blast -og
	
echo "Running orthofinder - TarBra, blast"
/path_to_orthofinder/orthofinder.py \
	-f /path_to_proteome/TarBra/ \
	-S blast -og

echo "Running orthofinder - TarCsa, blast"
/path_to_orthofinder/orthofinder.py \
	-f /path_to_proteome/TarCsa/ \
	-S blast -og

echo "Running orthofinder - TarAar, blast"
/path_to_orthofinder/orthofinder.py \
	-f /path_to_proteome/TarAar/ \
	-S blast -og
	
echo "Running orthofinder - SalBra, blast"
/path_to_orthofinder/orthofinder.py \
	-f /path_to_proteome/SalBra/ \
	-S blast -og
	
echo "Running orthofinder - SalCsa, blast"
/path_to_orthofinder/orthofinder.py \
	-f /path_to_proteome/SalCsa/ \
	-S blast -og
	
echo "Running orthofinder - SalAar, blast"
/path_to_orthofinder/orthofinder.py \
	-f /path_to_proteome/SalAar/ \
	-S blast -og
	
echo "Running orthofinder - BraCsa, blast"
/path_to_orthofinder/orthofinder.py \
	-f /path_to_proteome/BraCsa/ \
	-S blast -og
	
echo "Running orthofinder - BraAar, blast"
/path_to_orthofinder/orthofinder.py \
	-f /path_to_proteome/BraAar/ \
	-S blast -og

echo "Running orthofinder - CsaAar, blast"
/path_to_orthofinder/orthofinder.py \
	-f /path_to_proteome/CsaAar/ \
	-S blast -og