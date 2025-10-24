#!/bin/bash -l

module load vcftools


language='mayan'
sub_family='yucatec'

language_vcf=/group/jrigrp10/maize-linguistics/data/${language}/${language}.v4.imputed.recode.vcf

site_list_directory=/group/jrigrp10/maize-linguistics/results/${language}_gemma_output/

vcftools --vcf ${language_vcf} --snps ${site_list_directory}/${sub_family}_highlight.txt \
--out ${site_list_directory}/${sub_family}_genotype --recode 