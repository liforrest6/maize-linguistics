###HAYNIE UTO-AZTECAN###

library(rgdal)
library(ape)
library(ggtree)
library(treeio)
library(ggplot2)
library(raster)
library(rgeos)
library(plyr)
library(dplyr)
library(data.table)
library(tidyr)

source('languageFunctions.R')
#working directory = R Work
#load tree first

maize_seed_data <- read.delim("../data/SeeDs_passport_fulll_2017.04.25_16.08.12.txt")
#reading genetic data file for mapping genotype to accession number
genetic_data <- read.delim("../data/Target_sample_ID_mapping_DOI_and_information.txt", skip = 2)
genetic_data <- rename(genetic_data, c('bank_number' = 'BankAccessionNumber..description'))

filtered_maize_coord = filterMaizeCoordinatesbyCountry(maize_seed_data)

# make a copy of the dataframe for left joining/merging later
filtered_maize_coord_df = filtered_maize_coord

filtered_maize_coord = projectMaizeCoordinates(filtered_maize_coord)

## collect Haynie language polygons
haynie_lang<-readOGR(dsn=path.expand("../data/Languages_NAm"))

haynie_lang = processHayniePolygons(haynie_lang)

joined_coord_language = st_join(filtered_maize_coord, haynie_lang)

# proj4string(filtered_maize_coord) = proj4string(lang2)
# extracted_langs = raster::extract(lang2, filtered_maize_coord, df = T)
# # do a join of the copied dataframe with the names with the extracted languages
# coord_with_language_data <- merge(filtered_maize_coord_df, extracted_langs, by.x = 0, by.y = 'point.ID', all.x = T)
# View(coord_with_language_data)
# 
# #reading genetic data file
# genetic_data <- read.delim("Target_sample_ID_mapping_DOI_and_information.txt", skip = 2)
# View(genetic_data)
# genetic_data <- rename(genetic_data, c('bank_number' = 'BankAccessionNumber..description'))

#adding in genetic data
combined_coord_genetic_language <- merge(joined_coord_language, genetic_data, by.x= 'names', by.y = "bank_number", all.x = T)

#only include languages in the Uto-Aztecan family
combined_coord_genetic_language <- subset(combined_coord_genetic_language, family_name == 'Uto-Aztecan')

#remove languages that don't have matches
library(readxl)
SGpaired_language_vectors <- read_excel("../data/HaynieAztecanMatches.xlsx")

# $-Operator
vector_column_d <- SGpaired_language_vectors[ , "Fully Paired Polygons"]

coord_with_language_data_matched_languages <- subset(combined_coord_genetic_language, (Name %in% vector_column_d$"Fully Paired Polygons"))

#remove if N/A in language
coord_with_language_data_matched_languages_no_NA <- coord_with_language_data_matched_languages %>% drop_na("Name")

#cleaning up labels
coord_with_language_data_matched_languages_no_NA <- rename(coord_with_language_data_matched_languages_no_NA, c('Bank ID' = 'names'))
coord_with_language_data_matched_languages_no_NA <- rename(coord_with_language_data_matched_languages_no_NA, c('Polygon Language Name' = 'Name'))
coord_with_language_data_matched_languages_no_NA <- rename(coord_with_language_data_matched_languages_no_NA, c('GWAS ID' = "Sample.ID.of.DNA.from.single.plants.used.in.GWAS"))
# sort(unique(coord_with_language_data_matched_languages_no_NA$'Polygon Language Name'))
# unique(coord_with_language_data_matched_languages_no_NA$'Bank ID')
# View(coord_with_language_data_matched_languages_no_NA)

#filter out accessions missing SampleIDs (there were no duplicated Sample IDs this time around)
coord_with_language_data_matched_languages_no_NA <- coord_with_language_data_matched_languages_no_NA[complete.cases(coord_with_language_data_matched_languages_no_NA$'Sample.ID'), ]


#adding in tip/root to coord_with_language_data_matched_languages_no_NA dataframe
coord_with_language_data_matched_languages_no_NA <- coord_with_language_data_matched_languages_no_NA %>%
  mutate('Tree Tip/Root Name' = case_when(
    endsWith(`Polygon Language Name`, "Akimel O'odham") ~ "UA.TEPIMAN.TOHONO_OODHAM",
    endsWith(`Polygon Language Name`, "Central Nahuatl") ~ "UA.AZTECAN.CLASSICAL_NAHUATL",
    endsWith(`Polygon Language Name`, "Coatepec Nahuatl") ~ "UA.AZTECAN.NAHUATL_COATEPEC_COSTALES",
    endsWith(`Polygon Language Name`, "Cora") ~ "UA.CORACHOL.EL_NAYAR_CORA",
    endsWith(`Polygon Language Name`, "Eastern Durango Nahuatl") ~ "UA.AZTECAN.NAHUATL_SAN_PEDRO_JICORA",
    endsWith(`Polygon Language Name`, "Eastern Nahuatl") ~ "UA.AZTECAN.TABASCO_NAHUATL_CUPILCO",
    endsWith(`Polygon Language Name`, "Huichol") ~ "UA.CORACHOL.HUICHOL",
    endsWith(`Polygon Language Name`, "Mayo") ~ "UA.CAHITA.MAYO",
    endsWith(`Polygon Language Name`, "Michoacán Nahuatl") ~ "UA.AZTECAN.NAHUATL_POMARO_AQUILA",
    endsWith(`Polygon Language Name`, "Northern Puebla Nahuatl") ~ "UA.AZTECAN.HIGHLAND_PUEBLA_NAHUATL",
    endsWith(`Polygon Language Name`, "Northern Tepehuan") ~ "UA.TEPIMAN.NORTHERN_TEPEHUAN",
    endsWith(`Polygon Language Name`, "Opata") ~ "UA.CAHITA.OPATA",
    endsWith(`Polygon Language Name`, "Pima Bajo") ~ "UA.TEPIMAN.PIMA_BAJO",
    endsWith(`Polygon Language Name`, "Pipil") ~ "UA.AZTECAN.PIPIL",
    endsWith(`Polygon Language Name`, "Pochutec") ~ "UA.AZTECAN.POCHUTLA_NAHUATL",
    endsWith(`Polygon Language Name`, "Southern Tepehuan") ~ "UA.TEPIMAN.SOUTHERN_TEPEHUAN",
    endsWith(`Polygon Language Name`, "Tahue") ~ "UA.TUBAR.TUBAR",
    endsWith(`Polygon Language Name`, "Tarahumara") ~ "UA.TARAHUMARAN.CENTRAL_TARAHUMARA",
    endsWith(`Polygon Language Name`, "Tepecano") ~ "UA.TEPIMAN.TEPECANO",
    endsWith(`Polygon Language Name`, "Western Durango Nahuatl") ~ "UA.AZTECAN.NAHUATL_SAN_PEDRO_JICORA",
  ))



#adding numbers to haynie df
# coord_with_language_data_matched_languages_no_NA <- coord_with_language_data_matched_languages_no_NA %>%
#   mutate('Tree Tip/Root Associated Number' = case_when(
#     endsWith(`Polygon Language Name`, "Akimel O'odham") ~ 87,
#     endsWith(`Polygon Language Name`, "Central Nahuatl") ~ 23,
#     endsWith(`Polygon Language Name`, "Coatepec Nahuatl") ~ 5,
#     endsWith(`Polygon Language Name`, "Cora") ~ 73,
#     endsWith(`Polygon Language Name`, "Eastern Durango Nahuatl") ~ 11,
#     endsWith(`Polygon Language Name`, "Eastern Nahuatl") ~ 44,
#     endsWith(`Polygon Language Name`, "Huichol") ~ 72,
#     endsWith(`Polygon Language Name`, "Mayo") ~ 79,
#     endsWith(`Polygon Language Name`, "Michoacán Nahuatl") ~ 9,
#     endsWith(`Polygon Language Name`, "Northern Puebla Nahuatl") ~ 34,
#     endsWith(`Polygon Language Name`, "Northern Tepehuan") ~ 89,
#     endsWith(`Polygon Language Name`, "Opata") ~ 84,
#     endsWith(`Polygon Language Name`, "Pima Bajo") ~ 86,
#     endsWith(`Polygon Language Name`, "Pipil") ~ 41,
#     endsWith(`Polygon Language Name`, "Pochutec") ~ 20,
#     endsWith(`Polygon Language Name`, "Southern Tepehuan") ~ 90,
#     endsWith(`Polygon Language Name`, "Tahue") ~ 78,
#     endsWith(`Polygon Language Name`, "Tarahumara") ~ 77,
#     endsWith(`Polygon Language Name`, "Tepecano") ~ 91,
#     endsWith(`Polygon Language Name`, "Western Durango Nahuatl") ~ 11,
#   ))
# View(coord_with_language_data_matched_languages_no_NA)

aztecanHaynieMaster = coord_with_language_data_matched_languages_no_NA
aztecanHaynie_language_matrix_list = unlist(unique(aztecanHaynieMaster$`Tree Tip/Root Name`))
#converting all tree tip/root associated numbers to integers
#only needed to write csv #haynie_aztecan_all_accessions_df <- apply(haynie_aztecan_all_accessions_df,2,as.character)
# haynie_aztecan_all_accessions_df$"Tree Tip/Root Associated Number" <- lapply(haynie_aztecan_all_accessions_df$"Tree Tip/Root Associated Number", as.numeric)

#write.xlsx(simongreenhill_all_accessions_df)
#library("writexl")
#write_xlsx(simongreenhill_all_accessions_df,"simongreenhill_all_accessions_df.xlsx")
#write.csv(haynie_aztecan_all_accessions_df, "haynie_aztecan_all_accessions_df.csv")

#distance matrices:

aztecan_tree = read.newick('../data/trees/Uto-Aztecan.tre')
# mrca(mayan_tree, full = FALSE)

## calculate distances between all languages in language_matrix_list using the Mayan tree
aztecanHaynie_distance_list = calculateLanguageDistances(aztecanHaynie_language_matrix_list, aztecan_tree)

## ugly code but basically first take all of the non-MRCA languages and put it into its own matrix
aztecanHaynie_lang_distances = aztecanHaynie_distance_list[colnames(aztecanHaynie_distance_list), ]

aztecanHaynie_language_matrix = aztecanHaynie_lang_distances

## read updated GWAS IDs 
bankUpdated = read.csv('../data/selected_genotypeIDs.csv')
## limit to what we have SeeDs GWAS IDs for - unique 300 GWAS IDs, 521 entries for overlapping languages
aztecanHaynieMasterGWAS = merge(merge(aztecanHaynieMaster, maize_seed_data[c(5)], by.x = 'Bank ID', by.y = 'bank_number'), bankUpdated, by.x = 'GWAS ID', by.y = 'Sample')

aztecanHaynie_langlist = c()

## calculate the distance between each accession using the language lookup matrix
for(i in 1:length(aztecanHaynieMasterGWAS$`V1`)) {
  for (j in 1:i) {
    x = aztecanHaynie_language_matrix[aztecanHaynieMasterGWAS[[i, 'Tree Tip/Root Name']],
                                      aztecanHaynieMasterGWAS[[j, 'Tree Tip/Root Name']]]
    aztecanHaynie_langlist[[length(aztecanHaynie_langlist) + 1]] = data.frame(A = aztecanHaynieMasterGWAS[[i, 'V1']], B = aztecanHaynieMasterGWAS[[j, 'V1']], distance = x)
  }
}

library(tidyr) ## DO NOT LOAD UNTIL AFTER YOU USE RASTER::EXTRACT
library(gdata)
aztecanHaynie_accession_lang_matrix = generateAccessionMatrix(aztecanHaynie_langlist)

fwrite(as.data.frame(aztecanHaynie_accession_lang_matrix), '../results/aztecanHaynie_average_branch_language_distances.csv', row.names = T, quote = F)

aztecanHaynieMaster = extract(aztecanHaynieMaster, geometry, into = c('Longitude', 'Latitude'), '\\((.*),(.*)\\)', conv = T, remove = F)
aztecanHaynie_bankUpdated_latlon = merge(bankUpdated, aztecanHaynieMaster[,c('GWAS ID', 'Longitude', 'Latitude', 'Elevation')],by.x = 'Sample', by.y = 'GWAS ID')
aztecanHaynie_bankUpdated_latlon = aztecanHaynie_bankUpdated_latlon[!duplicated(aztecanHaynie_bankUpdated_latlon),]
aztecanHaynie_bankUpdated_latlon = aztecanHaynie_bankUpdated_latlon[which(aztecanHaynie_bankUpdated_latlon$V1 == rownames(aztecanHaynie_accession_lang_matrix)), ]

## calculate haversine distance - will still need elevation
library(geosphere)
aztecanHaynie_haversine_distance_df = calculateHaversineDistance(aztecanHaynie_bankUpdated_latlon)

fwrite(aztecanHaynie_haversine_distance_df, file= "../results/aztecanHaynie_haversine_distances.csv", row.names = TRUE, quote = F)

## elevation
aztecanHaynie_elevation_distance_df = calculateElevationDistance(aztecanHaynie_bankUpdated_latlon)

fwrite(aztecanHaynie_elevation_distance_df, '../results/aztecanHaynie_elevation_distances.csv', row.names = T, quote = F)


## for GWAS work

fwrite(aztecanHaynieMasterGWAS, '../results/aztecanHaynie_GWAS_allAccessions.csv', row.names = F, quote = F)
# write.csv(aztecanHaynieMasterGWAS, '../results/aztecanHaynie_GWAS_allAccessions_fix.csv', row.names = F, quote = F, eol = '\n')

# library(tibble)
# haynie_aztecan_distances_matrix <- haynie_aztecan_distances_matrix %>% rownames_to_column("Languages")
# haynie_aztecan_distances_matrix$"Languages" <- gsub("UA.AZTECAN.NAHUATL_POMARO_AQUILA", "9", haynie_aztecan_distances_matrix$"Languages")
# haynie_aztecan_distances_matrix$"Languages" <- gsub("UA.AZTECAN.CLASSICAL_NAHUATL", "23", haynie_aztecan_distances_matrix$"Languages")
# haynie_aztecan_distances_matrix$"Languages" <- gsub("UA.AZTECAN.POCHUTLA_NAHUATL", "20", haynie_aztecan_distances_matrix$"Languages")
# haynie_aztecan_distances_matrix$"Languages" <- gsub("UA.AZTECAN.NAHUATL_COATEPEC_COSTALES", "5", haynie_aztecan_distances_matrix$"Languages")
# haynie_aztecan_distances_matrix$"Languages" <- gsub("UA.AZTECAN.TABASCO_NAHUATL_CUPILCO", "44", haynie_aztecan_distances_matrix$"Languages")
# haynie_aztecan_distances_matrix$"Languages" <- gsub("UA.AZTECAN.NAHUATL_SAN_PEDRO_JICORA", "11", haynie_aztecan_distances_matrix$"Languages")
# haynie_aztecan_distances_matrix$"Languages" <- gsub("UA.CORACHOL.HUICHOL", "72", haynie_aztecan_distances_matrix$"Languages")
# haynie_aztecan_distances_matrix$"Languages" <- gsub("UA.CAHITA.OPATA", "84", haynie_aztecan_distances_matrix$"Languages")
# haynie_aztecan_distances_matrix$"Languages" <- gsub("UA.TEPIMAN.TOHONO_OODHAM", "87", haynie_aztecan_distances_matrix$"Languages")
# haynie_aztecan_distances_matrix$"Languages" <- gsub("UA.TEPIMAN.SOUTHERN_TEPEHUAN", "90", haynie_aztecan_distances_matrix$"Languages")
# haynie_aztecan_distances_matrix$"Languages" <- gsub("UA.AZTECAN.HIGHLAND_PUEBLA_NAHUATL", "34", haynie_aztecan_distances_matrix$"Languages")
# haynie_aztecan_distances_matrix$"Languages" <- gsub("UA.TARAHUMARAN.CENTRAL_TARAHUMARA", "77", haynie_aztecan_distances_matrix$"Languages")
# haynie_aztecan_distances_matrix$"Languages" <- gsub("UA.AZTECAN.PIPIL", "41", haynie_aztecan_distances_matrix$"Languages")
# haynie_aztecan_distances_matrix$"Languages" <- gsub("UA.TEPIMAN.PIMA_BAJO", "86", haynie_aztecan_distances_matrix$"Languages")
# haynie_aztecan_distances_matrix$"Languages" <- gsub("UA.TEPIMAN.TEPECANO", "91", haynie_aztecan_distances_matrix$"Languages")
# haynie_aztecan_distances_matrix$"Languages" <- gsub("UA.CORACHOL.EL_NAYAR_CORA", "73", haynie_aztecan_distances_matrix$"Languages")
# haynie_aztecan_distances_matrix$"Languages" <- gsub("UA.TUBAR.TUBAR", "78", haynie_aztecan_distances_matrix$"Languages")
# haynie_aztecan_distances_matrix$"Languages" <- gsub("UA.CAHITA.MAYO", "79", haynie_aztecan_distances_matrix$"Languages")
# haynie_aztecan_distances_matrix$"Languages" <- gsub("UA.TEPIMAN.NORTHERN_TEPEHUAN", "89", haynie_aztecan_distances_matrix$"Languages")
# haynie_aztecan_distances_matrix$"Languages" <- gsub("UA.AZTECAN.NAHUATL_SAN_PEDRO_JICORA", "11", haynie_aztecan_distances_matrix$"Languages")
# 
# haynie_aztecan_distances_matrix <- haynie_aztecan_distances_matrix[haynie_aztecan_distances_matrix$"Languages" %in% c("9","23","20","5", "44","11","72","84","87","90","34","77","41","86","91","73", "78", "79", "89"), ]
# haynie_aztecan_distances_matrix <- haynie_aztecan_distances_matrix %>% remove_rownames %>% column_to_rownames(var="Languages")
# View(haynie_aztecan_distances_matrix)
# dim(haynie_aztecan_distances_matrix)
#this code isn't needed? #rownames(simongreenhill_distances_matrix) = c("113","116","149","176", "179","183","185","186","187","190","195","196","34","41","54","72","84","85","87","90")

#match accessions to language distances
haynie_aztecan_language_distance_df = matrix(nrow= length(haynie_aztecan_all_accessions_df$"Sample.ID"), ncol= length(haynie_aztecan_all_accessions_df$"Sample.ID"))

for(i in 1:length(haynie_aztecan_all_accessions_df$"Sample.ID")) {
  for (j in 1:i) {
    ith = as.character(haynie_aztecan_all_accessions_df[[i, 'Tree Tip/Root Associated Number']])
    jth = as.character(haynie_aztecan_all_accessions_df[[j, 'Tree Tip/Root Associated Number']])
    x = haynie_aztecan_distances_matrix[ith, jth]
    print(c(ith, jth))
    print(x)
    haynie_aztecan_language_distance_df[i,j] = x[[1]]
    haynie_aztecan_language_distance_df[j,i] = x[[1]]
  }
}

#adding column and row names (the same) back to Seed IDs, saving as dataframe
colnames(haynie_aztecan_language_distance_df ) <- as.character(haynie_aztecan_all_accessions_df$"Sample.ID")
rownames(haynie_aztecan_language_distance_df) <- as.character(haynie_aztecan_all_accessions_df$"Sample.ID")
View(haynie_aztecan_language_distance_df)
haynie_aztecan_language_distance_df <- data.frame(haynie_aztecan_language_distance_df)

#saving Seed ID distance dataframe to CSV
fwrite(haynie_aztecan_language_distance_df, file= "haynie_aztecan_language_distance_df.csv", row.names = TRUE)

#elevation distance matrix
#get differences in elevation using AWS Terrain Tiles

library(elevatr)
##https://cran.r-project.org/web/packages/elevatr/vignettes/introduction_to_elevatr.html
haynie_aztecan_elevs = get_elev_point(haynie_aztecan_all_accessions_df[, c('Longitude', 'Latitude')], prj = 'EPSG:4326', src = 'aws')
haynie_aztecan_elevs_df = cbind(haynie_aztecan_all_accessions_df$'Sample.ID', haynie_aztecan_elevs@coords, haynie_aztecan_elevs@data)
colnames(haynie_aztecan_elevs_df) = c('Sample.ID', 'Longitude', 'Latitude', 'elevation', 'elev_units')

dim(haynie_aztecan_elevs_df)

haynie_aztecan_elev_diffs = matrix(ncol = 712, nrow = 712)

for (i in 1:712){
  for (j in 1:712){
    haynie_aztecan_elev_diffs[i, j] = abs(haynie_aztecan_elevs_df[i, 'elevation'] - haynie_aztecan_elevs_df[j, 'elevation'])
  }
}
View(haynie_aztecan_elev_diffs)

haynie_aztecan_elev_diffs_df = as.data.frame(haynie_aztecan_elev_diffs)

colnames(haynie_aztecan_elev_diffs_df) <- unlist(haynie_aztecan_elevs_df['Sample.ID'])
rownames(haynie_aztecan_elev_diffs_df) <- unlist(haynie_aztecan_elevs_df['Sample.ID'])
View(haynie_aztecan_elev_diffs_df)

fwrite(haynie_aztecan_elev_diffs_df, 'haynie_aztecan_elevation_distances.csv', row.names = T, quote = F)

#haversine distance matrix
library(geosphere)
library(readxl)

haynie_aztecan_haversine_distance_df = matrix(nrow= length(haynie_aztecan_all_accessions_df$"Sample.ID"), ncol= length(haynie_aztecan_all_accessions_df$"Sample.ID"))
for(i in 1:length(haynie_aztecan_all_accessions_df$"Sample.ID")) {
  for (j in 1:i) {
    x=
      distm(haynie_aztecan_all_accessions_df[i, c('Longitude', 'Latitude')], haynie_aztecan_all_accessions_df[j, c('Longitude', 'Latitude')], fun = distHaversine)
    haynie_aztecan_haversine_distance_df[i,j] = x
    haynie_aztecan_haversine_distance_df[j,i] = x
  }
}

haynie_aztecan_haversine_distance_df <- data.frame(haynie_aztecan_haversine_distance_df)
colnames(haynie_aztecan_haversine_distance_df) <- as.character(haynie_aztecan_all_accessions_df$"Sample.ID")
rownames(haynie_aztecan_haversine_distance_df) <- as.character(haynie_aztecan_all_accessions_df$"Sample.ID")
View(haynie_aztecan_haversine_distance_df)

fwrite(haynie_aztecan_haversine_distance_df, file= "haynie_aztecan_haversine_distance_df.csv", row.names = TRUE)