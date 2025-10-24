#!/bin/bash -l

#SBATCH -D /group/jrigrp10/maize-linguistics/scripts
#SBATCH -o /home/fli21/slurm-log/prepare-FEEMS=%j.txt
#SBATCH -e /home/fli21/slurm-log/prepare-FEEMS=%j.txt
#SBATCH -J feems
#SBATCH -t 4:00:00
#SBATCH --mem 10GB
#SBATCH -n 8
#SBATCH --account=jrigrp
#SBATCH --partition=high2
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=frrli@ucdavis.edu

dataset='NL'
data_dir='/group/jrigrp10/maize-linguistics/data/'
results_dir='/group/jrigrp10/maize-linguistics/results/'


## author: Forrest Li
## this script prepares data formats for running FEEMS, such as creating coordinate files for samples and plink format for genotypes.
## this script includes preparation for all SeeDs coordinates, for NL and Haynie separately, and NL and Haynie together.

## for ALL OF SEEDS

# create master coordinate file
python -c "import pandas as pd 
clim = pd.read_table('/group/jrigrp10/maize-linguistics/data/SeeDs_passport_fulll_2017.04.25_16.08.12.txt')
SeeD_accessions = pd.read_table('/group/jrigrp10/maize-linguistics/data/feems/all_Mesoamerican_SeeDs_accessions.txt', names=['Unique.ID', 'Duplicate'])
id_mapping = pd.read_csv('/group/jrigrp10/maize-linguistics/data/selected_genotypeIDs.csv')
seed_linguistics_samples_coordinates = id_mapping[id_mapping['V1'].isin(SeeD_accessions['Unique.ID'])].merge(clim, 
                                                             left_on = 'Sample', 
                                                             right_on = 'Sample ID of DNA from single plants used in GWAS')[['V1', 'Sample', 
                                                                                                                             'locations_longitude', 'locations_latitude', 'locations_elevation']]
seed_linguistics_samples_coordinates.to_csv('/group/jrigrp10/maize-linguistics/data//feems/seed_linguistics_coordinates.txt', sep = '\t', quoting = False, index = False, header = True)"


module load conda

python -c "import pandas as pd
clim = pd.read_table('${data_dir}/feems/seed_linguistics_coordinates.txt')
samples = pd.read_table('${data_dir}/feems/all_Mesoamerican_SeeDs_accessions.txt', names=['Unique.ID', 'Duplicate'])
coordinates = samples.merge(clim, left_on = 'Unique.ID', right_on = 'V1')[['locations_longitude', 'locations_latitude']]
coordinates.to_csv('${data_dir}/feems/all_samples_coordinates_SeeD.txt', sep = '\t', quoting = False, index = False, header = False)"

module load plink
module load bcftools

plink --vcf /group/jrigrp10/smambakk/language/zeaGBS20161020ZeaAlbertoTaxaAGPv4_forBeagle.imputed.vcf.vcf \
--make-bed \
--keep ${data_dir}/feems/all_Mesoamerican_SeeDs_accessions.txt \
--maf 0.01 \
--chr 1-10 \
--out ${data_dir}/feems/all_linguistics_samples_SeeD \
--allow-extra-chr \
--indep-pairwise 10kb 1 0.3



## for one dataset only
# cat ${results_dir}/${dataset}/mayan_${dataset}_accessions.txt ${results_dir}/${dataset}/aztecan_${dataset}_accessions.txt ${results_dir}/${dataset}/otomanguean_${dataset}_accessions.txt | 
# uniq | 
# sort |
# awk -v FS='\t' -v OFS='\t' 'BEGIN {OFS = FS} {print $0 FS $1}' |
# uniq > ${data_dir}/all_samples_list_${dataset}.txt

# ## create master coordinate file
# # python -c "import pandas as pd 
# # clim = pd.read_table('/group/jrigrp10/maize-linguistics/data/SeeDs_passport_fulll_2017.04.25_16.08.12.txt')
# # NL_accessions = pd.read_table('/group/jrigrp10/maize-linguistics/data/all_samples_list_NL.txt', names=['Unique.ID', 'Duplicate'])
# # Haynie_accessions = pd.read_table('/group/jrigrp10/maize-linguistics/data/all_samples_list_Haynie.txt', names=['Unique.ID', 'Duplicate'])
# # id_mapping = pd.read_csv('/group/jrigrp10/maize-linguistics/data/selected_genotypeIDs.csv')
# # all_linguistics_samples_coordinates = id_mapping[id_mapping['V1'].isin(NL_accessions['Unique.ID']) | 
# # id_mapping['V1'].isin(Haynie_accessions['Unique.ID'])].merge(clim, 
# #                                                              left_on = 'Sample', 
# #                                                              right_on = 'Sample ID of DNA from single plants used in GWAS')[['V1', 'Sample', 
# #                                                                                                                              'locations_longitude', 'locations_latitude']]
# # all_linguistics_samples_coordinates.to_csv('/group/jrigrp10/maize-linguistics/data//feems/all_linguistics_coordinates.txt', sep = '\t', quoting = False, index = False, header = True)"


# module load conda

# python -c "import pandas as pd
# clim = pd.read_table('${data_dir}/all_linguistics_coordinates.txt')
# samples = pd.read_table('${data_dir}/all_samples_list_${dataset}.txt', names=['Unique.ID', 'Duplicate'])
# coordinates = samples.merge(clim, left_on = 'Unique.ID', right_on = 'V1')[['locations_longitude', 'locations_latitude']]
# coordinates.to_csv('${data_dir}/feems/all_samples_coordinates_${dataset}.txt', sep = '\t', quoting = False, index = False, header = False)"

# module load plink
# module load bcftools

# plink --vcf /group/jrigrp10/smambakk/language/zeaGBS20161020ZeaAlbertoTaxaAGPv4_forBeagle.imputed.vcf.vcf \
# --make-bed \
# --keep ${data_dir}/all_samples_list_${dataset}.txt \
# --maf 0.05 \
# --chr 1-10 \
# --out ${data_dir}/feems/all_linguistics_samples_${dataset} \
# --allow-extra-chr \
# --indep-pairwise 10kb 1 0.3


## for both datasets
# cat ${results_dir}/NL/mayan_NL_accessions.txt \
# ${results_dir}/NL/aztecan_NL_accessions.txt \
# ${results_dir}/NL/otomanguean_NL_accessions.txt \
# ${results_dir}/Haynie/mayan_Haynie_accessions.txt \
# ${results_dir}/Haynie/aztecan_Haynie_accessions.txt \
# ${results_dir}/Haynie/otomanguean_Haynie_accessions.txt | 
# uniq | 
# sort |
# awk -v FS='\t' -v OFS='\t' 'BEGIN {OFS = FS} {print $0 FS $1}' |
# uniq > ${data_dir}/all_samples_list_HaynieNL.txt

# # create master coordinate file
# python -c "import pandas as pd 
# clim = pd.read_table('/group/jrigrp10/maize-linguistics/data/SeeDs_passport_fulll_2017.04.25_16.08.12.txt')
# NL_accessions = pd.read_table('/group/jrigrp10/maize-linguistics/data/all_samples_list_NL.txt', names=['Unique.ID', 'Duplicate'])
# Haynie_accessions = pd.read_table('/group/jrigrp10/maize-linguistics/data/all_samples_list_Haynie.txt', names=['Unique.ID', 'Duplicate'])
# id_mapping = pd.read_csv('/group/jrigrp10/maize-linguistics/data/selected_genotypeIDs.csv')
# all_linguistics_samples_coordinates = id_mapping[id_mapping['V1'].isin(NL_accessions['Unique.ID']) | 
# id_mapping['V1'].isin(Haynie_accessions['Unique.ID'])].merge(clim, 
#                                                              left_on = 'Sample', 
#                                                              right_on = 'Sample ID of DNA from single plants used in GWAS')[['V1', 'Sample', 
#                                                                                                                              'locations_longitude', 'locations_latitude', 'locations_elevation']]
# all_linguistics_samples_coordinates.to_csv('/group/jrigrp10/maize-linguistics/data//feems/all_linguistics_coordinates.txt', sep = '\t', quoting = False, index = False, header = True)"


# module load conda

# python -c "import pandas as pd
# clim = pd.read_table('${data_dir}/all_linguistics_coordinates.txt')
# samples = pd.read_table('${data_dir}/all_samples_list_HaynieNL.txt', names=['Unique.ID', 'Duplicate'])
# coordinates = samples.merge(clim, left_on = 'Unique.ID', right_on = 'V1')[['locations_longitude', 'locations_latitude']]
# coordinates.to_csv('${data_dir}/feems/all_samples_coordinates_HaynieNL.txt', sep = '\t', quoting = False, index = False, header = False)"

# module load plink
# module load bcftools

# plink --vcf /group/jrigrp10/smambakk/language/zeaGBS20161020ZeaAlbertoTaxaAGPv4_forBeagle.imputed.vcf.vcf \
# --make-bed \
# --keep ${data_dir}/all_samples_list_HaynieNL.txt \
# --maf 0.01 \
# --chr 1-10 \
# --out ${data_dir}/feems/all_linguistics_samples_HaynieNL \
# --allow-extra-chr \
# --indep-pairwise 10kb 1 0.3

