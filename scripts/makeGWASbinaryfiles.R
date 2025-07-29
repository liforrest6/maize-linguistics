library(readxl)
library(dplyr)
library(data.table)

language_families = c('aztecan', 'mayan', 'otomanguean')

## read sheets from the excel grouping the languages together
gwas_groupings_aztecan = read_excel('../data/GWAS_language_groupings.xlsx', sheet = 'Aztecan')
gwas_groupings_mayan = read_excel('../data/GWAS_language_groupings.xlsx', sheet = 'Mayan')
gwas_groupings_otomanguean = read_excel('../data/GWAS_language_groupings.xlsx', sheet = 'Otomanguean')

## double check the number of accessions per language grouping
gwas_groupings_aztecan %>% group_by(GWAS_grouping) %>% summarize(total_count = sum(Sample_size))
gwas_groupings_mayan %>% group_by(GWAS_grouping) %>% summarize(total_count = sum(Sample_size))
gwas_groupings_otomanguean %>% group_by(GWAS_grouping) %>% summarize(total_count = sum(Sample_size))

## read binary 0_1 files for presence of given language from haynie
binary_file_aztecan = vroom('../data/aztecan/0_1_aztecan.txt')
binary_file_mayan = vroom('../data/mayan/0_1_mayan.txt')
binary_file_otomanguean = vroom('../data/otomanguean/0_1_otomanguean.txt')

### create new columns based on groupings, use rowSums for all languages under a grouping since there is no overlap in Haynie ###########
compileGWAScolumns = function(grouping_name, binary_file, gwas_groupings_file) {
  colnames(binary_file) = c('Unique.ID', 'Polygon.Language', gwas_groupings_file$Language_name)
  individual_names = gwas_groupings_file[gwas_groupings_file$GWAS_grouping == grouping_name, ] %>% pull(Language_name)
  new_column = rowSums(binary_file[individual_names])
  new_column
}

aztecan_GWAS_binary = cbind(binary_file_aztecan %>% dplyr::select(Unique.ID),
                            as.data.frame(sapply(unique(gwas_groupings_aztecan$GWAS_grouping), 
                                          compileGWAScolumns, 
                                          binary_file = binary_file_aztecan, 
                                          gwas_groupings_file = gwas_groupings_aztecan)
))
mayan_GWAS_binary = cbind(binary_file_mayan %>% dplyr::select(Unique.ID),
                          as.data.frame(sapply(unique(gwas_groupings_mayan$GWAS_grouping), 
                                           compileGWAScolumns, 
                                           binary_file = binary_file_mayan, 
                                           gwas_groupings_file = gwas_groupings_mayan)
))
otomanguean_GWAS_binary = cbind(binary_file_otomanguean %>% dplyr::select(Unique.ID),
                                as.data.frame(sapply(unique(gwas_groupings_otomanguean$GWAS_grouping), 
                                           compileGWAScolumns, 
                                           binary_file = binary_file_otomanguean,
                                           gwas_groupings_file = gwas_groupings_otomanguean)
))

### write tables ###########
write.table(mayan_GWAS_binary, '/group/jrigrp10/maize-linguistics/data/mayan/0_1_GWAS_mayan.txt',sep = '\t', row.names = F)
write.table(aztecan_GWAS_binary, '/group/jrigrp10/maize-linguistics/data/aztecan/0_1_GWAS_aztecan.txt',sep = '\t', row.names = F)
write.table(otomanguean_GWAS_binary, '/group/jrigrp10/maize-linguistics/data/otomanguean/0_1_GWAS_otomanguean.txt',sep = '\t', row.names = F)


