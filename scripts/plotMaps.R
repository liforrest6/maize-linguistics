library(sp)
library(ggmap)
library(sf)
library(ggtree)
library(readxl)
library(ape)
library(treeio)
library(rgdal)

source('apikey.R')



mexico_terrain_kernels = get_stadiamap(c(left = -110, bottom = 10, right = -80, top = 27),
                                zoom = 6, maptype = 'stamen_toner_lite', color = 'bw', messaging = F)

mexico_terrain_kernels_transformed = ggmap_bbox(mexico_terrain_kernels)



### make ggtree plots for haynie polygons relative to the Jager dataset ####

aztecan_haynie_mapping <- read_excel("../data/master_paired_languages_2.xlsx", sheet = 'Aztecan-Haynie')
mayan_haynie_mapping <- read_excel("../data/master_paired_languages_2.xlsx", sheet = 'Mayan-Haynie')
otomanguean_haynie_mapping <- read_excel("../data/master_paired_languages_2.xlsx", sheet = 'Otomanguean-Haynie')

aztecan_nativelands_mapping <- read_excel("../data/master_paired_languages_2.xlsx", sheet = 'Aztecan')
mayan_nativelands_mapping <- read_excel("../data/master_paired_languages_2.xlsx", sheet = 'Mayan')
otomanguean_nativelands_mapping <- read_excel("../data/master_paired_languages_2.xlsx", sheet = 'Otomanguean')

native_mapping = read_excel("../data/master_paired_languages_2.xlsx", sheet = 'MRCA-Tip Mapping')[,c(2, 1, 3, 4, 5)]
haynie_mapping = read_excel("../data/master_paired_languages_2.xlsx", sheet = 'Haynie-Tip Mapping')[,c(2, 1, 3, 4, 5)]

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


ggtree(mayan_haynie_tree_filtered) %<+% haynie_mapping +
  geom_tiplab(aes(fill = factor(MRCA_abbreviated)),
              color = 'black',
              geom = 'label', 
              label.padding = unit(0.15, 'lines'),
              label.size = 0,
              size = 2) + 
  xlim(0,1) + 
  labs(fill = 'sub_grouping',
       title = 'Mayan Haynie')

ggtree(aztecan_haynie_tree_filtered) %<+% haynie_mapping +
  geom_tiplab(aes(fill = factor(MRCA_abbreviated)),
              color = 'black',
              geom = 'label', 
              label.padding = unit(0.15, 'lines'),
              label.size = 0,
              size = 2) + 
  xlim(0,1) + 
  labs(fill = 'sub_grouping',
       title = 'Aztecan Haynie')

ggtree(otomanguean_haynie_tree_filtered) %<+% haynie_mapping +
  geom_tiplab(aes(fill = factor(MRCA_abbreviated)),
              color = 'black',
              geom = 'label', 
              label.padding = unit(0.15, 'lines'),
              label.size = 0,
              size = 2) + 
  xlim(0,.3)+ 
  labs(fill = 'sub_grouping',
       title = 'Otomanguean Haynie')


aztecan_native_tree_filtered = keep.tip(aztecan_tree, 
                                 c(aztecan_nativelands_mapping$`Tree Name`[-which(!aztecan_nativelands_mapping$`Tree Name` %in% aztecan_tree$tip.label)],
                                   native_mapping %>% filter(`Language Family` == 'Aztecan') %>% pull(`Tip Name`)))
mayan_native_tree_filtered = keep.tip(mayan_tree, 
                               c(mayan_nativelands_mapping$`Tree Name`[-which(!mayan_nativelands_mapping$`Tree Name` %in% mayan_tree$tip.label)],
                                 native_mapping %>% filter(`Language Family` == 'Mayan') %>% pull(`Tip Name`)))

otomanguean_native_tree_filtered = keep.tip(otomanguean_tree, 
                                     c(otomanguean_nativelands_mapping$`Tree Name`[-which(!otomanguean_nativelands_mapping$`Tree Name` %in% otomanguean_tree$tip.label)],
                                       native_mapping %>% filter(`Language Family` == 'Otomanguean') %>% pull(`Tip Name`)))

ggtree(mayan_native_tree_filtered) %<+% native_mapping +
  geom_tiplab(aes(fill = factor(MRCA_abbreviated)),
              color = 'black',
              geom = 'label', 
              label.padding = unit(0.15, 'lines'),
              label.size = 0,
              size = 2) + 
  xlim(0,1) + 
  labs(fill = 'sub_grouping',
       title = 'Mayan NL')

ggtree(aztecan_native_tree_filtered) %<+% native_mapping +
  geom_tiplab(aes(fill = factor(MRCA_abbreviated)),
              color = 'black',
              geom = 'label', 
              label.padding = unit(0.15, 'lines'),
              label.size = 0,
              size = 2) + 
  xlim(0,1)+ 
  labs(fill = 'sub_grouping',
       title = 'Aztecan NL')

ggtree(otomanguean_native_tree_filtered) %<+% native_mapping +
  geom_tiplab(aes(fill = factor(MRCA_abbreviated)),
              color = 'black',
              geom = 'label', 
              label.padding = unit(0.15, 'lines'),
              label.size = 0,
              size = 2) + 
  xlim(0,1)+ 
  labs(fill = 'sub_grouping',
       title = 'Otomanguean NL')

### mapping polygons with accessions ################################################
nativelands_lang<-readOGR(dsn=path.expand("../data/indigenousLanguages_shp"))
nativelands_lang = processNativeLandsPolygons(nativelands_lang)

haynie_lang<-readOGR(dsn=path.expand("../data/Languages_NAm"))
haynie_lang = processHayniePolygons(haynie_lang)

## read accessions from master files
{
  otomanguean_NL_accessions = projectMaizeCoordinates(vroom('../results/NL/otomangueanMaster.csv'))
  otomanguean_haynie_accessions = projectMaizeCoordinates(vroom('../results/Haynie/otomangueanHaynieMaster.csv'))
  aztecan_NL_accessions = projectMaizeCoordinates(vroom('../results/NL/aztecanMaster.txt'))
  aztecan_haynie_accessions = projectMaizeCoordinates(vroom('../results/Haynie/aztecanHaynieMaster.csv'))
  mayan_NL_accessions = projectMaizeCoordinates(vroom('../results/NL/mayanMaster.csv'))
  mayan_haynie_accessions = projectMaizeCoordinates(vroom('../results/Haynie/mayanHaynieMaster.csv'))
}

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

## make pdf of all maps
pdf('../plots/all_maps.pdf', width = 9, height = 9)

## otomanguean_maps
ggmap(otomanguean_terrain_kernels_transformed) +
  coord_sf(crs = st_crs(3857)) +
  geom_sf(aes(fill = Name),
          data = st_transform(nativelands_lang %>% 
                                filter(Name %in% otomanguean_nativelands_mapping$`Polygon Name`), 
                              3857),
          inherit.aes = F,
          alpha = 0.4,
          size = 0.8) +
  geom_sf(data = st_transform(otomanguean_NL_accessions, 3857), 
          inherit.aes = F) +
  theme(legend.position = "bottom")+
  ggtitle('Otomanguean NL')

ggmap(otomanguean_terrain_kernels_transformed) +
  coord_sf(crs = st_crs(3857)) +
  geom_sf(aes(fill = Name),
          data = st_transform(haynie_lang %>% 
                                filter(Name %in% otomanguean_haynie_mapping$`Polygon Name`), 
                              3857),
          inherit.aes = F,
          alpha = 0.4,
          size = 0.8) +
  geom_sf(data = st_transform(otomanguean_haynie_accessions, 3857), 
          inherit.aes = F) +
  theme(legend.position = "bottom") +
  ggtitle('Otomanguean Haynie')

## aztecan maps
ggmap(aztecan_terrain_kernels_transformed) +
  coord_sf(crs = st_crs(3857)) +
  geom_sf(aes(fill = Name),
          data = st_transform(nativelands_lang %>% 
                                filter(Name %in% aztecan_nativelands_mapping$`Polygon Name`), 
                              3857),
          inherit.aes = F,
          alpha = 0.4) +
  geom_sf(data = st_transform(aztecan_NL_accessions, 3857), 
          inherit.aes = F,
          size = .8) +
  theme(legend.position = "bottom") +
  ggtitle('Aztecan NL')

ggmap(aztecan_terrain_kernels_transformed) +
  coord_sf(crs = st_crs(3857)) +
  geom_sf(aes(fill = Name),
          data = st_transform(haynie_lang %>% 
                                filter(Name %in% aztecan_haynie_mapping$`Polygon Name`), 
                              3857),
          inherit.aes = F,
          alpha = 0.4) +
  geom_sf(data = st_transform(aztecan_haynie_accessions, 3857), 
          inherit.aes = F,
          size = .8) +
  theme(legend.position = "bottom")+
  ggtitle('Aztecan Haynie')

## mayan maps
ggmap(mayan_terrain_kernels_transformed) +
  coord_sf(crs = st_crs(3857)) +
  geom_sf(aes(fill = Name),
          data = st_transform(nativelands_lang %>% 
                                filter(Name %in% mayan_nativelands_mapping$`Polygon Name`), 
                              3857),
          inherit.aes = F,
          alpha = 0.4) +
  geom_sf(data = st_transform(mayan_NL_accessions, 3857), 
          inherit.aes = F,
          size = .8) +
  theme(legend.position = "bottom")+
  ggtitle('Mayan NL')

ggmap(mayan_terrain_kernels_transformed) +
  coord_sf(crs = st_crs(3857)) +
  geom_sf(aes(fill = Name),
          data = st_transform(haynie_lang %>% 
                                filter(Name %in% mayan_haynie_mapping$`Polygon Name`), 
                              3857),
          inherit.aes = F,
          alpha = 0.4) +
  geom_sf(data = st_transform(mayan_haynie_accessions, 3857), 
          inherit.aes = F,
          size = .8) +
  theme(legend.position = "bottom") +
  ggtitle('Mayan Haynie')

dev.off()

### create table of unique accessions for each ########################################################
master_files = c('aztecan', 'aztecanHaynie', 'mayan', 'mayanHaynie', 'otomanguean', 'otomangueanHaynie')
accessions_per_file = lapply(master_files, function(master_file) read.csv(sprintf('../results/%sMaster.csv', master_file))$Unique.ID)
lapply(accessions_per_file, function(x) length(accessions_per_file))



