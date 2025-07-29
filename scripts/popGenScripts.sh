module load vcftools

## get Fst values using vcftools between different language families

vcftools --vcf /group/jrigrp11/seeds_vcf/zeaGBS20161020ZeaAlbertoTaxaAGPv4_forBeagle.imputed.vcf.vcf \
--maf 0.05 \
--weir-fst-pop /group/jrigrp10/maize-linguistics/results/NL/aztecan_NL_accessions.txt \
--weir-fst-pop /group/jrigrp10/maize-linguistics/results/NL/mayan_NL_accessions.txt \
--weir-fst-pop /group/jrigrp10/maize-linguistics/results/NL/otomanguean_NL_accessions.txt

## generate vcf files for individual language sub-groupings to perform selection tests such as for Yucatan Mayan

grep -e 'Yucatec Maya' -e 'Itzá' -e 'Mopán Maya' /group/jrigrp10/maize-linguistics/data/mayan/0_1_mayan.txt | cut -f 1  > /group/jrigrp10/maize-linguistics/data/mayan/Yucatec_Maya-only.txt

vcftools --vcf /group/jrigrp11/seeds_vcf/zeaGBS20161020ZeaAlbertoTaxaAGPv4_forBeagle.imputed.vcf.vcf \
--maf 0.01 \
--keep /group/jrigrp10/maize-linguistics/data/mayan/Yucatec_Maya-only.txt \
--recode \
--recode-INFO-all \
--out /group/jrigrp10/maize-linguistics/data/Yucatec_Maya-only

## run SweeD to perform selection test via SFS

SweeD -name Test-Yucatec -input /group/jrigrp10/maize-linguistics/data/Yucatec_Maya-only.recode.vcf -grid 40000 -folded