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
library(tibble)
library(gdata)
library(geosphere)
library(ggmap)
library(readxl)

source('languageFunctions.R')
source('apikey.R')

#get maize seed data from SeeDs passport
maize_seed_data <- read.delim("../data/SeeDs_passport_fulll_2017.04.25_16.08.12.txt")
bankUpdated = read.csv('../data/selected_genotypeIDs.csv')
colnames(bankUpdated) = c('Unique.ID', 'Sample')
#reading genetic data file for mapping genotype to accession number
genetic_data <- read.delim("../data/Target_sample_ID_mapping_DOI_and_information.txt", skip = 2)
genetic_data <- rename(genetic_data, c('bank_number' = 'BankAccessionNumber..description'))

filtered_maize_coord = filterMaizeCoordinatesbyCountry(maize_seed_data)


mayan_glottolog = read.csv('../data/david_centroids/Sheet 2-Mayan languages.csv')
aztecan_glottolog = read.csv('../data/david_centroids/Sheet 1-Uto-Aztecan Yuto-Nawan.csv')
otomanguean_glottolog = read.csv('../data/otomanguean_centroids.csv')

haynie_lang<-readOGR(dsn=path.expand("../data/Languages_NAm"))
haynie_lang = processHayniePolygons(haynie_lang)

# find overlap between maize coordinates and Native Lands polygons via spatial join
haynie_joined_coord_language = st_join(filtered_maize_coord, haynie_lang)

nativelands_lang<-readOGR(dsn=path.expand("../data/indigenousLanguages_shp"))
nativelands_lang = processNativeLandsPolygons(nativelands_lang)

# find overlap between maize coordinates and Native Lands polygons via spatial join
joined_coord_language = st_join(filtered_maize_coord, nativelands_lang)

ggmap_bbox <- function(map) {
  if (!inherits(map, "ggmap")) stop("map must be a ggmap object")
  # Extract the bounding box (in lat/lon) from the ggmap to a numeric vector, 
  # and set the names to what sf::st_bbox expects:
  map_bbox <- setNames(unlist(attr(map, "bb")), 
                       c("ymin", "xmin", "ymax", "xmax"))
  
  # Coonvert the bbox to an sf polygon, transform it to 3857, 
  # and convert back to a bbox (convoluted, but it works)
  bbox_3857 <- st_bbox(st_transform(st_as_sfc(st_bbox(map_bbox, crs = 4326)), 3857))
  
  # Overwrite the bbox of the ggmap object with the transformed coordinates 
  attr(map, "bb")$ll.lat <- bbox_3857["ymin"]
  attr(map, "bb")$ll.lon <- bbox_3857["xmin"]
  attr(map, "bb")$ur.lat <- bbox_3857["ymax"]
  attr(map, "bb")$ur.lon <- bbox_3857["xmax"]
  map
}

projectCoordinates = function(coord) {
  coordinates(coord) = ~ Longitude + Latitude
  proj4string(coord) = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
  coord <- st_as_sf(coord,coords = 1:2)
  coord <- coord %>% 
    st_transform(crs = 4326)
  return(coord)
}

mayan_glottolog = projectCoordinates(mayan_glottolog)
aztecan_glottolog = projectCoordinates(aztecan_glottolog)
otomanguean_glottolog = projectCoordinates(otomanguean_glottolog)

mexico_terrain_kernels = get_stadiamap(c(left = -110, bottom = 10, right = -80, top = 27),
                                       zoom = 6, maptype = 'stamen_toner_lite', color = 'bw', messaging = F)

mexico_terrain_kernels_transformed = ggmap_bbox(mexico_terrain_kernels)


ggmap(mexico_terrain_kernels_transformed) +
  coord_sf(crs = st_crs(3857)) +
  geom_sf(aes(fill = Name),
          data = st_transform(haynie_lang, 
                              3857),
          inherit.aes = F,
          alpha = 0.4) +
  geom_sf(data = st_transform(mayan_glottolog, 3857), 
          inherit.aes = F,
          size = .8) +
  theme(legend.position = "bottom")

mayan_haynie = read_excel('../data/master_paired_languages_2.xlsx', sheet = 'Mayan-Haynie')
aztecan_haynie = read_excel('../data/master_paired_languages_2.xlsx', sheet = 'Aztecan-Haynie')
otomanguean_haynie = read_excel('../data/master_paired_languages_2.xlsx', sheet = 'Otomanguean-Haynie')


mayan_haynie_overlap = st_join(haynie_lang %>% filter(Name %in% mayan_haynie$`Polygon Name`), mayan_glottolog, left = F) %>% pull(Name.x) %>% unique()
aztecan_haynie_overlap = st_join(haynie_lang %>% filter(Name %in% aztecan_haynie$`Polygon Name`), aztecan_glottolog, left = F) %>% pull(Name.x) %>% unique()
otomanguean_haynie_overlap = st_join(haynie_lang %>% filter(Name %in% otomanguean_haynie$`Polygon Name`), otomanguean_glottolog, left = F) %>% pull(Name.x) %>% unique()

setdiff(mayan_haynie$`Polygon Name`, mayan_haynie_overlap) %>% length()
setdiff(aztecan_haynie$`Polygon Name`, aztecan_haynie_overlap) %>% length()
setdiff(otomanguean_haynie$`Polygon Name`, otomanguean_haynie_overlap) %>% length()

mayan_NL = read_excel('../data/master_paired_languages_2.xlsx', sheet = 'Mayan')
aztecan_NL = read_excel('../data/master_paired_languages_2.xlsx', sheet = 'Aztecan')
otomanguean_NL = read_excel('../data/master_paired_languages_2.xlsx', sheet = 'Otomanguean')

mayan_overlap = st_join(nativelands_lang %>% filter(Name %in% mayan_NL$`Polygon Name`), mayan_glottolog, left = F) %>% pull(Name.x) %>% unique()
aztecan_overlap = st_join(nativelands_lang %>% filter(Name %in% aztecan_NL$`Polygon Name`), aztecan_glottolog, left = F) %>% pull(Name.x) %>% unique()
otomanguean_overlap = st_join(nativelands_lang %>% filter(Name %in% otomanguean_NL$`Polygon Name`), otomanguean_glottolog, left = F) %>% pull(Name.x) %>% unique()

setdiff(mayan_NL$`Polygon Name`, mayan_overlap) %>% length() 
setdiff(aztecan_NL$`Polygon Name`, aztecan_overlap) %>% length()
setdiff(otomanguean_NL$`Polygon Name`, otomanguean_overlap) %>% length()


ggmap(mexico_terrain_kernels_transformed) +
  coord_sf(crs = st_crs(3857)) +
  geom_sf(aes(fill = Name),
          data = st_transform(nativelands_lang %>% filter(Name %in% setdiff(aztecan_NL$`Polygon Name`, aztecan_overlap)), 
                              3857),
          inherit.aes = F,
          alpha = 0.4) +
  geom_sf(data = st_transform(aztecan_glottolog ,
                                                         3857), 
          inherit.aes = F,
          size = .8) +
  theme(legend.position = "bottom")
