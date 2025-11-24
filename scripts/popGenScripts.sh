module load vcftools

## author: Forrest Li
## scripts to calculate population genetics statistics such as genome-wide Fst and SweeD selection scans

## get Fst values using vcftools between different language families

# vcftools --vcf /group/jrigrp11/seeds_vcf/zeaGBS20161020ZeaAlbertoTaxaAGPv4_forBeagle.imputed.vcf.vcf \
# --maf 0.05 \
# --weir-fst-pop /group/jrigrp10/maize-linguistics/results/NL/aztecan_NL_accessions.txt \
# --weir-fst-pop /group/jrigrp10/maize-linguistics/results/NL/mayan_NL_accessions.txt \
# --weir-fst-pop /group/jrigrp10/maize-linguistics/results/NL/otomanguean_NL_accessions.txt \
# --out /group/jrigrp10/maize-linguistics/results/NL/NL.weir

# vcftools --vcf /group/jrigrp11/seeds_vcf/zeaGBS20161020ZeaAlbertoTaxaAGPv4_forBeagle.imputed.vcf.vcf \
# --maf 0.05 \
# --weir-fst-pop /group/jrigrp10/maize-linguistics/results/Haynie/aztecan_Haynie_accessions.txt \
# --weir-fst-pop /group/jrigrp10/maize-linguistics/results/Haynie/mayan_Haynie_accessions.txt \
# --weir-fst-pop /group/jrigrp10/maize-linguistics/results/Haynie/otomanguean_Haynie_accessions.txt \
# --out /group/jrigrp10/maize-linguistics/results/Haynie/Haynie.weir

## generate vcf files for individual language sub-groupings to perform selection tests such as 

##################################################
## for Yucatan Mayan

# grep -e 'Yucatec Maya' -e 'Itzá' -e 'Mopán Maya' /group/jrigrp10/maize-linguistics/data/mayan/0_1_mayan.txt | cut -f 1  > /group/jrigrp10/maize-linguistics/data/mayan/Yucatec_Maya-only.txt

# # vcftools --vcf /group/jrigrp11/seeds_vcf/zeaGBS20161020ZeaAlbertoTaxaAGPv4_forBeagle.imputed.vcf.vcf \
# # --maf 0.01 \
# # --keep /group/jrigrp10/maize-linguistics/data/mayan/Yucatec_Maya-only.txt \
# # --chr 1 \
# # --recode \
# # --recode-INFO-all \
# # --out /group/jrigrp10/maize-linguistics/data/mayan/Yucatec_Maya-only-chr01

# # ## run SweeD to perform selection test via SFS

# # SweeD -name Yucatec-01 -input /group/jrigrp10/maize-linguistics/data/mayan/Yucatec_Maya-only-chr01.recode.vcf -grid 40000 -folded

# # Yucatec chr 1-9

# for i in $(seq 2 9);
# do 
# 	## run SweeD to perform selection test via SFS
# 	echo "doing ${i}"

# 	vcftools --vcf /group/jrigrp11/seeds_vcf/zeaGBS20161020ZeaAlbertoTaxaAGPv4_forBeagle.imputed.vcf.vcf \
# 	--maf 0.01 \
# 	--keep /group/jrigrp10/maize-linguistics/data/mayan/Yucatec_Maya-only.txt \
# 	--chr ${i} \
# 	--recode \
# 	--recode-INFO-all \
# 	--out /group/jrigrp10/maize-linguistics/data/mayan/Yucatec_Maya-only-chr0${i}

# 	## run SweeD to perform selection test via SFS

# 	SweeD -name Yucatec-0${i} -input /group/jrigrp10/maize-linguistics/data/mayan/Yucatec_Maya-only-chr0${i}.recode.vcf -grid 40000 -folded
# done

# ## Yucatec chr10

# vcftools --vcf /group/jrigrp11/seeds_vcf/zeaGBS20161020ZeaAlbertoTaxaAGPv4_forBeagle.imputed.vcf.vcf \
# --maf 0.01 \
# --keep /group/jrigrp10/maize-linguistics/data/mayan/Yucatec_Maya-only.txt \
# --chr 10 \
# --recode \
# --recode-INFO-all \
# --out /group/jrigrp10/maize-linguistics/data/mayan/Yucatec_Maya-only-chr10

# ## run SweeD to perform selection test via SFS

# SweeD -name Yucatec-10 -input /group/jrigrp10/maize-linguistics/data/mayan/Yucatec_Maya-only-chr10.recode.vcf -grid 40000 -folded


##################################################
## for Otomi

grep -e 'Querétaro Otomi' -e 'Tilapa Otomi' -e 'Mezquital Otomi' -e 'Ixtenco Otomi' /group/jrigrp10/maize-linguistics/data/otomanguean/0_1_otomanguean.txt | cut -f 1  > /group/jrigrp10/maize-linguistics/data/otomanguean/Otomi-only.txt


# Otomi chr4

vcftools --vcf /group/jrigrp11/seeds_vcf/zeaGBS20161020ZeaAlbertoTaxaAGPv4_forBeagle.imputed.vcf.vcf \
--maf 0.01 \
--keep /group/jrigrp10/maize-linguistics/data/otomanguean/Otomi-only.txt \
--chr 4 \
--recode \
--recode-INFO-all \
--out /group/jrigrp10/maize-linguistics/data/otomanguean/Otomi-only-chr04

SweeD -name Otomi-04 -input /group/jrigrp10/maize-linguistics/data/otomanguean/Otomi-only-chr04.recode.vcf -grid 40000 -folded

# Otomi chr5

vcftools --vcf /group/jrigrp11/seeds_vcf/zeaGBS20161020ZeaAlbertoTaxaAGPv4_forBeagle.imputed.vcf.vcf \
--maf 0.01 \
--keep /group/jrigrp10/maize-linguistics/data/otomanguean/Otomi-only.txt \
--chr 5 \
--recode \
--recode-INFO-all \
--out /group/jrigrp10/maize-linguistics/data/otomanguean/Otomi-only-chr05

SweeD -name Otomi-05 -input /group/jrigrp10/maize-linguistics/data/otomanguean/Otomi-only-chr05.recode.vcf -grid 40000 -folded

# Otomi chr9

vcftools --vcf /group/jrigrp11/seeds_vcf/zeaGBS20161020ZeaAlbertoTaxaAGPv4_forBeagle.imputed.vcf.vcf \
--maf 0.01 \
--keep /group/jrigrp10/maize-linguistics/data/otomanguean/Otomi-only.txt \
--chr 9 
--recode \
--recode-INFO-all \
--out /group/jrigrp10/maize-linguistics/data/otomanguean/Otomi-only-chr09

SweeD -name Otomi-09 -input /group/jrigrp10/maize-linguistics/data/otomanguean/Otomi-only-chr09.recode.vcf -grid 40000 -folded

##################################################
## for Cora/Huichol

grep -e 'Cora' -e 'Huichol' /group/jrigrp10/maize-linguistics/data/aztecan/0_1_aztecan.txt | cut -f 1  > /group/jrigrp10/maize-linguistics/data/aztecan/Cora-Huichol-only.txt

# Cora chr4

vcftools --vcf /group/jrigrp11/seeds_vcf/zeaGBS20161020ZeaAlbertoTaxaAGPv4_forBeagle.imputed.vcf.vcf \
--maf 0.01 \
--keep /group/jrigrp10/maize-linguistics/data/aztecan/Cora-Huichol-only.txt \
--chr 4 \
--recode \
--recode-INFO-all \
--out /group/jrigrp10/maize-linguistics/data/aztecan/Cora-Huichol-only-chr04

SweeD -name Cora-Huichol-04 -input /group/jrigrp10/maize-linguistics/data/aztecan/Cora-Huichol-only-chr04.recode.vcf -grid 40000 -folded

##################################################
## for Chol

grep -e 'Chol' -e 'Tabasco Chontal' -e 'Chortí' -e 'Cholti' /group/jrigrp10/maize-linguistics/data/mayan/0_1_mayan.txt | cut -f 1  > /group/jrigrp10/maize-linguistics/data/mayan/Chol-only.txt

# Chol chr5

vcftools --vcf /group/jrigrp11/seeds_vcf/zeaGBS20161020ZeaAlbertoTaxaAGPv4_forBeagle.imputed.vcf.vcf \
--maf 0.01 \
--keep /group/jrigrp10/maize-linguistics/data/mayan/Chol-only.txt \
--chr 5 \
--recode \
--recode-INFO-all \
--out /group/jrigrp10/maize-linguistics/data/mayan/Chol-only-chr05

SweeD -name Chol-05 -input /group/jrigrp10/maize-linguistics/data/mayan/Chol-only-chr05.recode.vcf -grid 40000 -folded

# Chol chr7

vcftools --vcf /group/jrigrp11/seeds_vcf/zeaGBS20161020ZeaAlbertoTaxaAGPv4_forBeagle.imputed.vcf.vcf \
--maf 0.01 \
--keep /group/jrigrp10/maize-linguistics/data/mayan/Chol-only.txt \
--chr 7 \
--recode \
--recode-INFO-all \
--out /group/jrigrp10/maize-linguistics/data/mayan/Chol-only-chr07

SweeD -name Chol-07 -input /group/jrigrp10/maize-linguistics/data/mayan/Chol-only-chr07.recode.vcf -grid 40000 -folded

#################################################
## for Pipil

grep -e 'Eastern Nahuatl' -e 'Northern Puebla Nahuatl' -e 'Pipil' /group/jrigrp10/maize-linguistics/data/aztecan/0_1_aztecan.txt | cut -f 1  > /group/jrigrp10/maize-linguistics/data/aztecan/Pipil-only.txt

## Pipl chr9

vcftools --vcf /group/jrigrp11/seeds_vcf/zeaGBS20161020ZeaAlbertoTaxaAGPv4_forBeagle.imputed.vcf.vcf \
--maf 0.01 \
--keep /group/jrigrp10/maize-linguistics/data/aztecan/Pipil-only.txt \
--chr 9 \
--recode \
--recode-INFO-all \
--out /group/jrigrp10/maize-linguistics/data/aztecan/Pipil-only-chr09

SweeD -name Pipil-09 -input /group/jrigrp10/maize-linguistics/data/aztecan/Pipil-only-chr09.recode.vcf -grid 40000 -folded

##################################################
## for Mangue
grep -e 'Mangue' -e 'Chiapenec' -e 'Mephaa'  /group/jrigrp10/maize-linguistics/data/otomanguean/0_1_otomanguean.txt | cut -f 1  > /group/jrigrp10/maize-linguistics/data/otomanguean/Mangue-only.txt

# Mangue chr1

vcftools --vcf /group/jrigrp11/seeds_vcf/zeaGBS20161020ZeaAlbertoTaxaAGPv4_forBeagle.imputed.vcf.vcf \
--maf 0.01 \
--keep /group/jrigrp10/maize-linguistics/data/otomanguean/Mangue-only.txt \
--chr 1 \
--recode \
--recode-INFO-all \
--out /group/jrigrp10/maize-linguistics/data/otomanguean/Mangue-only-chr01 

SweeD -name Mangue-01 -input /group/jrigrp10/maize-linguistics/data/otomanguean/Mangue-only-chr01.recode.vcf -grid 40000 -folded

# Mangue chr3

vcftools --vcf /group/jrigrp11/seeds_vcf/zeaGBS20161020ZeaAlbertoTaxaAGPv4_forBeagle.imputed.vcf.vcf \
--maf 0.01 \
--keep /group/jrigrp10/maize-linguistics/data/otomanguean/Mangue-only.txt \
--chr 3 \
--recode \
--recode-INFO-all \
--out /group/jrigrp10/maize-linguistics/data/otomanguean/Mangue-only-chr03

SweeD -name Mangue-03 -input /group/jrigrp10/maize-linguistics/data/otomanguean/Mangue-only-chr03.recode.vcf -grid 40000 -folded


##################################################
## for Matlazinca
grep -e 'San Francisco Matlatzinca' -e 'Atzingo Matlatzinca' -e 'Michoacán Mazahua' -e 'Central Mazahua'  /group/jrigrp10/maize-linguistics/data/otomanguean/0_1_otomanguean.txt | cut -f 1  > /group/jrigrp10/maize-linguistics/data/otomanguean/Matlazinca-only.txt

# Matlazinca chr8
## MATLATZINCA WAS INITIALLY MISSPELLED BUT DID NOT RERUN THESE FILES SO KEPT AS MATLAZINCA [sic]

vcftools --vcf /group/jrigrp11/seeds_vcf/zeaGBS20161020ZeaAlbertoTaxaAGPv4_forBeagle.imputed.vcf.vcf \
--maf 0.01 \
--keep /group/jrigrp10/maize-linguistics/data/otomanguean/Matlazinca-only.txt \
--chr 8 \
--recode \
--recode-INFO-all \
--out /group/jrigrp10/maize-linguistics/data/otomanguean/Matlazinca-only-chr08

SweeD -name Matlazinca-08 -input /group/jrigrp10/maize-linguistics/data/otomanguean/Matlazinca-only-chr08.recode.vcf -grid 40000 -folded

