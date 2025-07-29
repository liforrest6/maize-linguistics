###OTOMANGUEAN/ORIGINAL###
## adapted to calculate average branch lengths
## author: Forrest Li

#loading packages
library(rgdal)
library(ape)
library(ggtree)
library(treeio)
library(ggplot2)
library(rgeos)
# library(plyr)
library(dplyr)
library(sf)
library(vroom)
library(data.table)
library(raster)
library(RRphylo)
library(readxl)
library(tibble)
library(tidyr)
library(gdata)
library(geosphere)

source('languageFunctions.R')

#get maize seed data from SeeDs passport
maize_seed_data <- read.delim("../data/SeeDs_passport_fulll_2017.04.25_16.08.12.txt")
bankUpdated = read.csv('../data/selected_genotypeIDs.csv')
colnames(bankUpdated) = c('Unique.ID', 'Sample')
#reading genetic data file for mapping genotype to accession number
genetic_data <- read.delim("../data/Target_sample_ID_mapping_DOI_and_information.txt", skip = 2)
genetic_data <- rename(genetic_data, c('bank_number' = 'BankAccessionNumber..description'))

filtered_maize_coord = filterMaizeCoordinatesbyCountry(maize_seed_data)

# make a copy of the dataframe for left joining/merging later
filtered_maize_coord_df = filtered_maize_coord

filtered_maize_coord = projectMaizeCoordinates(filtered_maize_coord)

# read indigenous languages shapes from Native Lands
haynie_lang<-readOGR(dsn=path.expand("../data/Languages_NAm"))
haynie_lang = processHayniePolygons(haynie_lang)

# find overlap between maize coordinates and Native Lands polygons via spatial join
haynie_joined_coord_language = st_join(filtered_maize_coord, haynie_lang)

haynie_joined_coord_language = haynie_joined_coord_language[!is.na(haynie_joined_coord_language$Name),]
haynie_joined_coord_language = rename(haynie_joined_coord_language, c('Bank ID' = 'names',
                                                        'Polygon Language Name' = 'Name',
                                                        'GWAS ID' = "Sample.ID.of.DNA.from.single.plants.used.in.GWAS",
                                                        'Language Family' = 'family_nam'))

######  Import all of the mapping sheets ################################
aztecan_haynie_mapping <- read_excel("../data/master_paired_languages_2.xlsx", sheet = 'Aztecan-Haynie')
mayan_haynie_mapping <- read_excel("../data/master_paired_languages_2.xlsx", sheet = 'Mayan-Haynie')
otomanguean_haynie_mapping <- read_excel("../data/master_paired_languages_2.xlsx", sheet = 'Otomanguean-Haynie')

mrca_mapping = read_excel("../data/master_paired_languages_2.xlsx", sheet = 'Haynie-Tip Mapping')


aztecan_tree = read.newick('../data/trees/Uto-Aztecan.tre')
mayan_tree = read.newick('../data/trees/Mayan.tre')
otomanguean_tree = read.newick('../data/trees/Otomanguean.tre')

######  Native Lands Otomanguean ################################
otomanguean_accessions = merge(haynie_joined_coord_language, 
                               otomanguean_haynie_mapping, 
                               by.x = 'Polygon Language Name',
                               by.y = 'Polygon Name')

##### special case:  some polygons in otomanguean don't have corresponding tree tip languages - going to remove these
otomanguean_accessions = otomanguean_accessions %>% filter(!`Tree Name` == 'No comparable')

otomangueanMaster = otomanguean_accessions
otomangueanMaster$LangName = gsub('OM\\..*\\.', '', otomangueanMaster$`Tree Name`, fixed = F)
otomangueanMaster$LangName = gsub('OM.', '', otomangueanMaster$`LangName`, fixed = T)

otomangueanMasterGWAS = merge(merge(otomangueanMaster, maize_seed_data[c(5)], by.x = 'Bank ID', by.y = 'bank_number'), 
                              bankUpdated, by.x = 'GWAS ID', by.y = 'Sample')

## take all single matches to tree as well as multiple matches to tree
otomanguean_tip_languages = gsub('OM.', '',
                                 gsub('OM\\..*\\.', '',
                                      mrca_mapping %>%
                                     filter(`Language Family` == 'Otomanguean') %>%
                                     pull(`Tip Name`) %>%
                                     unique(),
                                 fixed = F),
                           fixed = T)
otomanguean_mrca_languages = gsub('OM.', '',
                                  gsub('OM\\..*\\.', '',
                                       mrca_mapping %>%
                                         filter(`Language Family` == 'Otomanguean') %>%
                                         pull(`MRCA Name`) %>%
                                         unique(),
                                       fixed = F),
                                  fixed = T)

otomanguean_levenshtein_matrix_list = unlist(c(unique(otomangueanMaster$`LangName`), 
                                            otomanguean_tip_languages))
otomanguean_levenshtein_distances = obtainLevenshteinDistances(otomanguean_levenshtein_matrix_list)
otomanguean_levenshtein_matrix = generateLanguageMatrix(otomanguean_levenshtein_distances, otomanguean_mrca_languages, mrca_mapping)


## calculate the distance between each accession using the language lookup matrix
otomanguean_levenshtein_langlist = createAccessionDistanceList(otomangueanMasterGWAS, otomanguean_levenshtein_matrix)
otomanguean_accession_levenshtein_matrix = generateAccessionMatrix(otomanguean_levenshtein_langlist)


#### for average branch length ########################################################################
otomanguean_tip_languages_extended = mrca_mapping %>%
  filter(`Language Family` == 'Otomanguean') %>%
  pull(`Tip Name`) %>%
  unique()

otomanguean_branch_distances = calculateBranchDistances(c(unique(otomangueanMaster$`Tree Name`),
                                                          otomanguean_tip_languages_extended), otomanguean_tree)
colnames(otomanguean_branch_distances) = gsub('OM.', '',
                                              gsub('OM\\..*\\.', '',
                                                   colnames(otomanguean_branch_distances),
                                                   fixed = F),
                                              fixed = T)
rownames(otomanguean_branch_distances) = colnames(otomanguean_branch_distances)
otomanguean_branch_distance_matrix = generateLanguageMatrix(otomanguean_branch_distances, otomanguean_mrca_languages, mrca_mapping)
otomanguean_branch_langlist = createAccessionDistanceList(otomangueanMasterGWAS, otomanguean_branch_distance_matrix)
otomanguean_accession_branch_lang_matrix = generateAccessionMatrix(otomanguean_branch_langlist)

########################################################################

otomanguean_bankUpdated_latlon = extract(otomangueanMasterGWAS, geometry, 
                                     into = c('Longitude', 'Latitude'), '\\((.*),(.*)\\)', 
                                     conv = T, remove = F)[,c('GWAS ID', 'Unique.ID', 
                                                              'Longitude', 'Latitude', 'Elevation',
                                                              'Tree Name', 'LangName', 'Polygon Language Name')]
otomanguean_bankUpdated_latlon = otomanguean_bankUpdated_latlon[!duplicated(otomanguean_bankUpdated_latlon$Unique.ID),]
## calculate haversine and elevation distance
otomanguean_haversine_distance_df = calculateHaversineDistance(otomanguean_bankUpdated_latlon)
otomanguean_elevation_distance_df = calculateElevationDistance(otomanguean_bankUpdated_latlon)
otomanguean_elevationSigned_distance_df = calculateElevationDistanceSigned(otomanguean_bankUpdated_latlon)

######  Haynie Aztecan ################################
aztecan_accessions = merge(haynie_joined_coord_language, 
                           aztecan_haynie_mapping, 
                           by.x = 'Polygon Language Name',
                           by.y = 'Polygon Name')

aztecanMaster = aztecan_accessions
aztecanMaster$LangName = gsub('UA\\..*\\.', '', aztecanMaster$`Tree Name`, fixed = F)
aztecanMaster$LangName = gsub('UA.', '', aztecanMaster$`LangName`, fixed = T)

aztecanMasterGWAS = merge(merge(aztecanMaster, maize_seed_data[c(5)], by.x = 'Bank ID', by.y = 'bank_number'), 
                          bankUpdated, by.x = 'GWAS ID', by.y = 'Sample')

aztecan_tip_languages = gsub('UA.', '',
                                 gsub('UA\\..*\\.', '',
                                      mrca_mapping %>%
                                        filter(`Language Family` == 'Aztecan') %>%
                                        pull(`Tip Name`) %>%
                                        unique(),
                                      fixed = F),
                                 fixed = T)
aztecan_mrca_languages = gsub('UA.', '',
                                  gsub('UA\\..*\\.', '',
                                       mrca_mapping %>%
                                         filter(`Language Family` == 'Aztecan') %>%
                                         pull(`MRCA Name`) %>%
                                         unique(),
                                       fixed = F),
                                  fixed = T)
aztecan_levenshtein_matrix_list = unlist(c(unique(aztecanMaster$`LangName`), 
                                        aztecan_tip_languages))
aztecan_levenshtein_distances = obtainLevenshteinDistances(aztecan_levenshtein_matrix_list)
aztecan_levenshtein_matrix = generateLanguageMatrix(aztecan_levenshtein_distances, aztecan_mrca_languages, mrca_mapping)

## calculate the distance between each accession using the language lookup matrix
aztecan_levenshtein_langlist = createAccessionDistanceList(aztecanMasterGWAS, aztecan_levenshtein_matrix)
aztecan_accession_levenshtein_matrix = generateAccessionMatrix(aztecan_levenshtein_langlist)

#### for average branch length ########################################################################
aztecan_tip_languages_extended = mrca_mapping %>%
  filter(`Language Family` == 'Aztecan') %>%
  pull(`Tip Name`) %>%
  unique()

aztecan_branch_distances = calculateBranchDistances(c(unique(aztecanMaster$`Tree Name`),
                                                      aztecan_tip_languages_extended), aztecan_tree)
colnames(aztecan_branch_distances) = gsub('UA.', '',
                                          gsub('UA\\..*\\.', '',
                                               colnames(aztecan_branch_distances),
                                               fixed = F),
                                          fixed = T)
rownames(aztecan_branch_distances) = colnames(aztecan_branch_distances)
aztecan_branch_distance_matrix = generateLanguageMatrix(aztecan_branch_distances, aztecan_mrca_languages, mrca_mapping)
aztecan_branch_langlist = createAccessionDistanceList(aztecanMasterGWAS, aztecan_branch_distance_matrix)
aztecan_accession_branch_lang_matrix = generateAccessionMatrix(aztecan_branch_langlist)

########################################################################

aztecan_bankUpdated_latlon = extract(aztecanMasterGWAS, geometry, 
                                   into = c('Longitude', 'Latitude'), '\\((.*),(.*)\\)', 
                                   conv = T, remove = F)[,c('GWAS ID', 'Unique.ID', 
                                                            'Longitude', 'Latitude', 'Elevation',
                                                            'Tree Name', 'LangName', 'Polygon Language Name')]
aztecan_bankUpdated_latlon = aztecan_bankUpdated_latlon[!duplicated(aztecan_bankUpdated_latlon$Unique.ID),]
## calculate haversine and elevation distance 
aztecan_haversine_distance_df = calculateHaversineDistance(aztecan_bankUpdated_latlon)
aztecan_elevation_distance_df = calculateElevationDistance(aztecan_bankUpdated_latlon)
aztecan_elevationSigned_distance_df = calculateElevationDistanceSigned(aztecan_bankUpdated_latlon)


######  Native Lands Mayan ################################
mayan_accessions = merge(haynie_joined_coord_language, 
                         mayan_haynie_mapping, 
                           by.x = 'Polygon Language Name',
                           by.y = 'Polygon Name')

mayanMaster = mayan_accessions
# mayanMaster$LangName = gsub('May.MAYAN\\..*\\.', '', mayanMaster$`Tree Name`, fixed = F)
mayanMaster$LangName = gsub('May.MAYAN.', '', mayanMaster$`Tree Name`, fixed = F)

mayanMasterGWAS = merge(merge(mayanMaster, maize_seed_data[c(5)], by.x = 'Bank ID', by.y = 'bank_number'), bankUpdated, by.x = 'GWAS ID', by.y = 'Sample')

mayan_tip_languages = gsub('May.MAYAN.', '',
                           mrca_mapping %>%
                                    filter(`Language Family` == 'Mayan') %>%
                                    pull(`Tip Name`) %>%
                                    unique(),
                                  fixed = F)
mayan_mrca_languages = gsub('May.MAYAN.', '',
                            mrca_mapping %>%
                                     filter(`Language Family` == 'Mayan') %>%
                                     pull(`MRCA Name`) %>%
                                     unique(),
                                   fixed = F)
mayan_levenshtein_matrix_list = unlist(c(unique(mayanMaster$`LangName`), 
                                      mayan_tip_languages))
mayan_levenshtein_distances = obtainLevenshteinDistances(mayan_levenshtein_matrix_list)
mayan_levenshtein_matrix = generateLanguageMatrix(mayan_levenshtein_distances, mayan_mrca_languages, mrca_mapping)

## calculate the distance between each accession using the language lookup matrix
mayan_levenshtein_langlist = createAccessionDistanceList(mayanMasterGWAS, mayan_levenshtein_matrix)

## make it long-form dataframe: most of the code below was to make this faster since it took a long time before (much faster now)
mayan_accession_levenshtein_matrix = generateAccessionMatrix(mayan_levenshtein_langlist)

#### for average branch length ########################################################################
mayan_tip_languages_extended = mrca_mapping %>%
  filter(`Language Family` == 'Mayan') %>%
  pull(`Tip Name`) %>%
  unique()

mayan_branch_distances = calculateBranchDistances(c(unique(mayanMaster$`Tree Name`),
                                                    mayan_tip_languages_extended), mayan_tree)
colnames(mayan_branch_distances) = gsub('May.MAYAN.', '',
                                        colnames(mayan_branch_distances),
                                        fixed = F)
rownames(mayan_branch_distances) = colnames(mayan_branch_distances)
mayan_branch_distance_matrix = generateLanguageMatrix(mayan_branch_distances, mayan_mrca_languages, mrca_mapping)
mayan_branch_langlist = createAccessionDistanceList(mayanMasterGWAS, mayan_branch_distance_matrix)
mayan_accession_branch_lang_matrix = generateAccessionMatrix(mayan_branch_langlist)

########################################################################


mayan_bankUpdated_latlon = extract(mayanMasterGWAS, geometry, 
                                   into = c('Longitude', 'Latitude'), '\\((.*),(.*)\\)', 
                                   conv = T, remove = F)[,c('GWAS ID', 'Unique.ID', 
                                                            'Longitude', 'Latitude', 'Elevation',
                                                            'Tree Name', 'LangName', 'Polygon Language Name')]
mayan_bankUpdated_latlon = mayan_bankUpdated_latlon[!duplicated(mayan_bankUpdated_latlon$Unique.ID),]

mayan_haversine_distance_df = calculateHaversineDistance(mayan_bankUpdated_latlon)
mayan_elevation_distance_df = calculateElevationDistance(mayan_bankUpdated_latlon)
mayan_elevationSigned_distance_df = calculateElevationDistanceSigned(mayan_bankUpdated_latlon)


#### write functions ########################################

## aztecan
fwrite(as.data.frame(aztecan_accession_levenshtein_matrix), '../results/Haynie/aztecanHaynie_levenshtein_language_distances.csv', row.names = T, quote = F)
fwrite(as.data.frame(aztecan_accession_branch_lang_matrix), '../results/Haynie/aztecanHaynie_branch_language_distances.csv', row.names = T, quote = F)
write.csv(extract(aztecanMasterGWAS, geometry, into = c('Longitude', 'Latitude'), '\\((.*),(.*)\\)', conv = T), '../results/Haynie/aztecanHaynieMaster.csv', row.names = F, quote = F)
fwrite(aztecan_haversine_distance_df, file= "../results/Haynie/aztecanHaynie_haversine_distances.csv", row.names = TRUE, quote = F)
fwrite(aztecan_elevation_distance_df, '../results/Haynie/aztecanHaynie_elevation_distances.csv', row.names = T, quote = F)

## otomanguean
fwrite(as.data.frame(otomanguean_accession_levenshtein_matrix), '../results/Haynie/otomangueanHaynie_levenshtein_language_distances.csv', row.names = T, quote = F)
fwrite(as.data.frame(otomanguean_accession_branch_lang_matrix), '../results/Haynie/otomangueanHaynie_branch_language_distances.csv', row.names = T, quote = F)
write.csv(extract(otomangueanMasterGWAS, geometry, into = c('Longitude', 'Latitude'), '\\((.*),(.*)\\)', conv = T), '../results/Haynie/otomangueanHaynieMaster.csv', row.names = F, quote = F)
fwrite(otomanguean_haversine_distance_df, file= "../results/Haynie/otomangueanHaynie_haversine_distances.csv", row.names = TRUE, quote = F)
fwrite(otomanguean_elevation_distance_df, '../results/Haynie/otomangueanHaynie_elevation_distances.csv', row.names = T, quote = F)

## mayan
fwrite(as.data.frame(mayan_accession_levenshtein_matrix), '../results/Haynie/mayanHaynie_levenshtein_language_distances.csv', row.names = T, quote = F)
fwrite(as.data.frame(mayan_accession_branch_lang_matrix), '../results/Haynie/mayanHaynie_branch_language_distances.csv', row.names = T, quote = F)
write.csv(extract(mayanMasterGWAS, geometry, into = c('Longitude', 'Latitude'), '\\((.*),(.*)\\)', conv = T), '../results/Haynie/mayanHaynieMaster.csv', row.names = F, quote = F)
fwrite(mayan_haversine_distance_df, file= "../results/Haynie/mayanHaynie_haversine_distances.csv", row.names = TRUE, quote = F)
fwrite(mayan_elevation_distance_df, '../results/Haynie/mayanHaynie_elevation_distances.csv', row.names = T, quote = F)

fwrite(aztecan_elevationSigned_distance_df, '../results/Haynie/aztecanHaynie_elevationSigned_distances.csv', row.names = T, quote = F)
fwrite(otomanguean_elevationSigned_distance_df, '../results/Haynie/otomangueanHaynie_elevationSigned_distances.csv', row.names = T, quote = F)
fwrite(mayan_elevationSigned_distance_df, '../results/Haynie/mayanHaynie_elevationSigned_distances.csv', row.names = T, quote = F)

