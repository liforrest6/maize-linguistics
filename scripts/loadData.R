## Load data for map plotting and analyses
## author: Forrest Li

source('apikey.R')
source('languageFunctions.R')


# Generate tiles for maps ---------------------------------------------------------------------

# General Mexico terrain kernels
mexico_terrain_kernels = get_stadiamap(c(left = -110, bottom = 10, right = -80, top = 27),
                                       zoom = 6, maptype = 'stamen_toner_lite', color = 'bw', messaging = F)

mexico_terrain_kernels_transformed = ggmap_bbox(mexico_terrain_kernels)

mesoamerica_terrain_kernels = get_stadiamap(c(left = -114, bottom = 8, right = -82, top = 33),
                                            zoom = 6, maptype = 'stamen_toner_lite', color = 'bw', messaging = F)
mesoamerica_terrain_kernels_transformed = ggmap_bbox(mesoamerica_terrain_kernels)

## otomanguean nativelands c(left = -105, bottom = 8, right = -82, top = 24)
otomanguean_terrain_kernels = get_stadiamap(c(left = -105, bottom = 8, right = -82, top = 24),
                                            zoom = 6, maptype = 'stamen_toner_lite', color = 'bw', messaging = F)
otomanguean_terrain_kernels_transformed = ggmap_bbox(otomanguean_terrain_kernels)
## aztecan nativelands c(left = -114, bottom = 10, right = -87, top = 33)
aztecan_terrain_kernels = get_stadiamap(c(left = -114, bottom = 10, right = -87, top = 33),
                                        zoom = 6, maptype = 'stamen_toner_lite', color = 'bw', messaging = F)
aztecan_terrain_kernels_transformed = ggmap_bbox(aztecan_terrain_kernels)

## mayan nativelands c(left = -101, bottom = 12, right = -85, top = 25)
mayan_terrain_kernels = get_stadiamap(c(left = -101, bottom = 12, right = -85, top = 25),
                                      zoom = 6, maptype = 'stamen_toner_lite', color = 'bw', messaging = F)
mayan_terrain_kernels_transformed = ggmap_bbox(mayan_terrain_kernels)

americas_terrain_kernels = get_stadiamap(c(left = -120, bottom = -40, right = -35, top = 33),
                                         zoom = 5, maptype = 'stamen_toner_lite', color = 'bw', messaging = F)
americas_terrain_kernels_transformed = ggmap_bbox(americas_terrain_kernels)


# Load polygons -------------------------------------------------------------------------------

nativelands_lang<-readOGR(dsn=path.expand("../data/indigenousLanguages_shp"))
nativelands_lang = processNativeLandsPolygons(nativelands_lang)

haynie_lang<-readOGR(dsn=path.expand("../data/Languages_NAm"))
haynie_lang = processHayniePolygons(haynie_lang)

# Load mapping files --------------------------------------------------------------------------

aztecan_haynie_mapping <- read_excel("../data/master_paired_languages_2.xlsx", sheet = 'Aztecan-Haynie')
mayan_haynie_mapping <- read_excel("../data/master_paired_languages_2.xlsx", sheet = 'Mayan-Haynie')
otomanguean_haynie_mapping <- read_excel("../data/master_paired_languages_2.xlsx", sheet = 'Otomanguean-Haynie')

aztecan_nativelands_mapping <- read_excel("../data/master_paired_languages_2.xlsx", sheet = 'Aztecan')
mayan_nativelands_mapping <- read_excel("../data/master_paired_languages_2.xlsx", sheet = 'Mayan')
otomanguean_nativelands_mapping <- read_excel("../data/master_paired_languages_2.xlsx", sheet = 'Otomanguean')

native_mapping = read_excel("../data/master_paired_languages_2.xlsx", sheet = 'MRCA-Tip Mapping')[,c(2, 1, 3, 4, 5)]
haynie_mapping = read_excel("../data/master_paired_languages_2.xlsx", sheet = 'Haynie-Tip Mapping')[,c(2, 1, 3, 4, 5)]


# Load trees ----------------------------------------------------------------------------------

aztecan_tree = read.newick('../data/trees/Uto-Aztecan.tre')
mayan_tree = read.newick('../data/trees/Mayan.tre')
otomanguean_tree = read.newick('../data/trees/Otomanguean.tre')

aztecan_haynie_tree_filtered = keep.tip(aztecan_tree, 
                                        c(aztecan_haynie_mapping$`Tree Name`[-which(!aztecan_haynie_mapping$`Tree Name` %in% aztecan_tree$tip.label)],
                                          haynie_mapping %>% filter(`Language Family` == 'Aztecan') %>% pull(`Tip Name`)))
mayan_haynie_tree_filtered = keep.tip(mayan_tree, 
                                      c(mayan_haynie_mapping$`Tree Name`[-which(!mayan_haynie_mapping$`Tree Name` %in% mayan_tree$tip.label)],
                                        haynie_mapping %>% filter(`Language Family` == 'Mayan') %>% pull(`Tip Name`)))

otomanguean_haynie_tree_filtered = keep.tip(otomanguean_tree, 
                                            c(otomanguean_haynie_mapping$`Tree Name`[-which(!otomanguean_haynie_mapping$`Tree Name` %in% otomanguean_tree$tip.label)],
                                              haynie_mapping %>% filter(`Language Family` == 'Otomanguean') %>% pull(`Tip Name`)))
aztecan_native_tree_filtered = keep.tip(aztecan_tree, 
                                        c(aztecan_nativelands_mapping$`Tree Name`[-which(!aztecan_nativelands_mapping$`Tree Name` %in% aztecan_tree$tip.label)],
                                          native_mapping %>% filter(`Language Family` == 'Aztecan') %>% pull(`Tip Name`)))
mayan_native_tree_filtered = keep.tip(mayan_tree, 
                                      c(mayan_nativelands_mapping$`Tree Name`[-which(!mayan_nativelands_mapping$`Tree Name` %in% mayan_tree$tip.label)],
                                        native_mapping %>% filter(`Language Family` == 'Mayan') %>% pull(`Tip Name`)))

otomanguean_native_tree_filtered = keep.tip(otomanguean_tree, 
                                            c(otomanguean_nativelands_mapping$`Tree Name`[-which(!otomanguean_nativelands_mapping$`Tree Name` %in% otomanguean_tree$tip.label)],
                                              native_mapping %>% filter(`Language Family` == 'Otomanguean') %>% pull(`Tip Name`)))


# Load accession lists ------------------------------------------------------------------------
{
  otomanguean_NL_accessions = projectMaizeCoordinates(vroom('../results/NL/otomangueanMaster.csv'))
  otomanguean_haynie_accessions = projectMaizeCoordinates(vroom('../results/Haynie/otomangueanHaynieMaster.csv'))
  aztecan_NL_accessions = projectMaizeCoordinates(vroom('../results/NL/aztecanMaster.txt'))
  aztecan_haynie_accessions = projectMaizeCoordinates(vroom('../results/Haynie/aztecanHaynieMaster.csv'))
  mayan_NL_accessions = projectMaizeCoordinates(vroom('../results/NL/mayanMaster.csv'))
  mayan_haynie_accessions = projectMaizeCoordinates(vroom('../results/Haynie/mayanHaynieMaster.csv'))
}

all_haynie_accessions = rbind(otomanguean_haynie_accessions, aztecan_haynie_accessions, mayan_haynie_accessions)
all_NL_accessions = rbind(otomanguean_NL_accessions, aztecan_NL_accessions, mayan_NL_accessions)
