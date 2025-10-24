module load vcftools

## author: Forrest Li
## scripts to calculate population genetics statistics such as genome-wide Fst and SweeD selection scans

## get Fst values using vcftools between different language families

vcftools --vcf /group/jrigrp11/seeds_vcf/zeaGBS20161020ZeaAlbertoTaxaAGPv4_forBeagle.imputed.vcf.vcf \
--maf 0.05 \
--weir-fst-pop /group/jrigrp10/maize-linguistics/results/NL/aztecan_NL_accessions.txt \
--weir-fst-pop /group/jrigrp10/maize-linguistics/results/NL/mayan_NL_accessions.txt \
--weir-fst-pop /group/jrigrp10/maize-linguistics/results/NL/otomanguean_NL_accessions.txt \
--out /group/jrigrp10/maize-linguistics/results/NL/NL.weir

vcftools --vcf /group/jrigrp11/seeds_vcf/zeaGBS20161020ZeaAlbertoTaxaAGPv4_forBeagle.imputed.vcf.vcf \
--maf 0.05 \
--weir-fst-pop /group/jrigrp10/maize-linguistics/results/Haynie/aztecan_Haynie_accessions.txt \
--weir-fst-pop /group/jrigrp10/maize-linguistics/results/Haynie/mayan_Haynie_accessions.txt \
--weir-fst-pop /group/jrigrp10/maize-linguistics/results/Haynie/otomanguean_Haynie_accessions.txt \
--out /group/jrigrp10/maize-linguistics/results/Haynie/Haynie.weir

## generate vcf files for individual language sub-groupings to perform selection tests such as 

##################################################
## for Yucatan Mayan

grep -e 'Yucatec Maya' -e 'Itzá' -e 'Mopán Maya' /group/jrigrp10/maize-linguistics/data/mayan/0_1_mayan.txt | cut -f 1  > /group/jrigrp10/maize-linguistics/data/mayan/Yucatec_Maya-only.txt

vcftools --vcf /group/jrigrp11/seeds_vcf/zeaGBS20161020ZeaAlbertoTaxaAGPv4_forBeagle.imputed.vcf.vcf \
--maf 0.01 \
--keep /group/jrigrp10/maize-linguistics/data/mayan/Yucatec_Maya-only.txt \
--chr 1 \
--recode \
--recode-INFO-all \
--out /group/jrigrp10/maize-linguistics/data/mayan/Yucatec_Maya-only-chr01

## run SweeD to perform selection test via SFS

SweeD -name Yucatec-01 -input /group/jrigrp10/maize-linguistics/data/mayan/Yucatec_Maya-only-chr01.recode.vcf -grid 40000 -folded

##################################################
## for Otomi

grep -e 'Querétaro Otomi' -e 'Tilapa Otomi' -e 'Mezquital Otomi' -e 'Ixtenco Otomi' /group/jrigrp10/maize-linguistics/data/otomanguean/0_1_otomanguean.txt | cut -f 1  > /group/jrigrp10/maize-linguistics/data/otomanguean/Otomi-only.txt

vcftools --vcf /group/jrigrp11/seeds_vcf/zeaGBS20161020ZeaAlbertoTaxaAGPv4_forBeagle.imputed.vcf.vcf \
--maf 0.01 \
--keep /group/jrigrp10/maize-linguistics/data/otomanguean/Otomi-only.txt \
--chr 9 
--recode \
--recode-INFO-all \
--out /group/jrigrp10/maize-linguistics/data/otomanguean/Otomi-only-chr09

## run SweeD to perform selection test via SFS

SweeD -name Otomi-09 -input /group/jrigrp10/maize-linguistics/data/otomanguean/Otomi-only-chr09.recode.vcf -grid 40000 -folded

vcftools --vcf /group/jrigrp11/seeds_vcf/zeaGBS20161020ZeaAlbertoTaxaAGPv4_forBeagle.imputed.vcf.vcf \
--maf 0.01 \
--keep /group/jrigrp10/maize-linguistics/data/otomanguean/Otomi-only.txt \
--chr 4 \
--recode \
--recode-INFO-all \
--out /group/jrigrp10/maize-linguistics/data/otomanguean/Otomi-only-chr04

##################################################
## for Cora/Huichol

grep -e 'Cora' -e 'Huichol' /group/jrigrp10/maize-linguistics/data/aztecan/0_1_aztecan.txt | cut -f 1  > /group/jrigrp10/maize-linguistics/data/aztecan/Cora-Huichol-only.txt

vcftools --vcf /group/jrigrp11/seeds_vcf/zeaGBS20161020ZeaAlbertoTaxaAGPv4_forBeagle.imputed.vcf.vcf \
--maf 0.01 \
--keep /group/jrigrp10/maize-linguistics/data/aztecan/Cora-Huichol-only.txt \
--chr 4 \
--recode \
--recode-INFO-all \
--out /group/jrigrp10/maize-linguistics/data/aztecan/Cora-Huichol-only-chr09

## run SweeD to perform selection test via SFS

SweeD -name Cora-Huichol-04 -input /group/jrigrp10/maize-linguistics/data/aztecan/Cora-Huichol-only-chr09.recode.vcf -grid 40000 -folded

