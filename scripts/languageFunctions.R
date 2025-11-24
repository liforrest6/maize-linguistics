## languageFunctions.R
## Authors: Forrest Li, Luke Sparreo
## Functions associated with maize linguistics project
## Intended to be loaded via source()

library(dplyr)
library(vroom)
library(sp)
library(ggmap)
library(sf)
library(ggtree)
library(readxl)
library(ape)
library(treeio)
library(rgdal)
library(ggtext)
library(cowplot)
library(grid)
library(gridExtra)
library(png)
library(ggplotify)
library(stringr)
library(scattermore)
library(scales)

## edit bbox for projections on ggmap
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

## general function for projecting coordinates with Latitude and Longitude columns
projectCoordinates = function(coord) {
  coordinates(coord) = ~ longitude + latitude
  proj4string(coord) = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
  coord <- st_as_sf(coord,coords = 1:2)
  return(coord)
}

## given passport data of maize data, limit to central america countries, and keep only those with latitude and longitude
filterMaizeCoordinatesbyCountry = function(maize_seed_data) {
  #filter by countries to just central America and with complete country code data = 9478
  maize_seed_data <- maize_seed_data[maize_seed_data$countries_country_name %in% c('MEXICO','BELIZE','GUATEMALA','EL SALVADOR','HONDURAS','NICARAGUA','COSTA RICA'), ]
  maize_seed_data_filtered_by_countrycode <- maize_seed_data[complete.cases(maize_seed_data[ , c('locations_latitude', 'locations_longitude')]), ]
  
  #remove any accessions without GWAS data
  maize_seed_data_filtered_by_countrycode <- maize_seed_data_filtered_by_countrycode[maize_seed_data_filtered_by_countrycode$Sample.ID.of.DNA.from.single.plants.used.in.GWAS %like% "SEEDGWAS", ]
  
  #putting coordinates to languages
  filtered_maize_coord = data.frame(Longitude = maize_seed_data_filtered_by_countrycode$locations_longitude, Latitude = maize_seed_data_filtered_by_countrycode$locations_latitude, Elevation= maize_seed_data_filtered_by_countrycode$locations_elevation, Sample.ID.of.DNA.from.single.plants.used.in.GWAS= maize_seed_data_filtered_by_countrycode$Sample.ID.of.DNA.from.single.plants.used.in.GWAS)
  names = maize_seed_data_filtered_by_countrycode$bank_number
  filtered_maize_coord$names = names
  return(filtered_maize_coord)
}

## project coordinates to WGS84 from lat/long coordinates and set it to spatial
projectMaizeCoordinates = function(filtered_maize_coord) {
  coordinates(filtered_maize_coord) = ~ Longitude + Latitude
  proj4string(filtered_maize_coord) = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
  filtered_maize_coord <- st_as_sf(filtered_maize_coord,coords = 1:2)
  filtered_maize_coord <- filtered_maize_coord %>% 
    st_transform(crs = 4326)
  return(filtered_maize_coord)
}
# use terra::extract to grab languages associated with each polygon that overlaps at specific coordinate of maize
# extracted_langs = terra::extract(lang, filtered_maize_coord)
# do a join of the maize metadataframe with the names of the extracted languages (one maize to many languages) = 15271
# coord_with_language_data <- merge(filtered_maize_coord_df, extracted_langs, by.x = 0, by.y = 'id.y', all.x = T)

## process polygons from Native land so that they are spatial polygons
processNativeLandsPolygons = function(lang) {
  polys = attr(lang,'polygons')
  names(polys)<-lang$Name
  ## DELETE THE ALEUTIANS
  lang = lang[-c(180), ]
  
  
  proj4string(lang) = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
  lang <- st_as_sf(lang,coords = 1:2)
  lang <- lang %>% 
    st_transform(crs = 4326)
  lang = lang[st_is_valid(lang),]
  return(lang)
}

## process polygons from Haynie so that they are spatial polygons
processHayniePolygons = function(lang) {
  polys = attr(lang,'polygons')
  names(polys)<-lang$Name
  ## DELETE THE ALEUTIANS
  # lang = lang[-c(180), ]
  
  
  proj4string(lang) = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
  lang <- st_as_sf(lang,coords = 1:2)
  lang <- lang %>% 
    st_transform(crs = 4326)
  lang = lang[st_is_valid(lang),]
  lang$long_glott <- gsub("\r\n", "", lang$long_glott)
  return(lang)
}


## calculate distances between all languages in language_matrix_list using the Mayan tree
calculateBranchDistances = function(language_matrix_list, tree) {
  language_matrix_list = language_matrix_list[language_matrix_list %in% tree$tip.label]
  distances = sapply(language_matrix_list, function(x) {
    print(x)
    distance = distNodes(tree = tree, x, clus = 0)
    # append(distance_list, distance[, 2])
    distance[,2]
  })
  distances = distances[language_matrix_list, language_matrix_list]
  return(distances)
}

## generate matrix of most recent common ancestor between languages
makeMRCAMatrix = function(lang_distance, mrca_langs) {
  mrca_matrix = matrix(nrow = length(mrca_langs), ncol = length(mrca_langs))
  for(i in 1:length(mrca_langs)) {
    for(j in 1:length(mrca_langs)) {
      if(i == j) mrca_matrix[i, j] = 0
      else {
        comparison = lang_distance[mrca_langs[[i]], mrca_langs[[j]]]
        avg_avg = data.frame(matrix(unlist(comparison), ncol = length(mrca_langs[[j]]))) %>% rowMeans() %>% mean()
        mrca_matrix[i, j] = avg_avg
        mrca_matrix[j, i] = avg_avg
      }
    }
  }
  return(mrca_matrix)
}

## generate matrix between accessions based on accession list with associated language
generateAccessionMatrix = function(langlist) {
  language_distance_long = do.call(rbind, langlist)
  colnames(language_distance_long) = c('A', 'B', 'distance')
  language_long_averages = language_distance_long %>% summarize(distance_avg = mean(distance), .by = c(A, B))
  accession_language_matrix = language_long_averages %>% pivot_wider(names_from = B, id_cols = A, values_from = distance_avg)
  # accession_test = accession_language_matrix
  # rownames(accession_test) = accession_test$A
  accession_matrix = as.matrix(accession_language_matrix[, 2:ncol(accession_language_matrix)])
  upperTriangle(accession_matrix) = lowerTriangle(accession_matrix, byrow = T)
  rownames(accession_matrix) = colnames(accession_matrix)
  accession_matrix = as.data.frame(accession_matrix)
  return(accession_matrix)
}

## calculate distances between accessions using haversine based on lat/lon
calculateHaversineDistance = function(bankUpdated_latlon) {
  haversine_distance = matrix(nrow= nrow(bankUpdated_latlon), ncol= nrow(bankUpdated_latlon))
  for(i in 1:nrow(bankUpdated_latlon)) {
    for (j in 1:i) {
      x=
        distm(unlist(bankUpdated_latlon[i, c('Longitude', 'Latitude')]), unlist(bankUpdated_latlon[j, c('Longitude', 'Latitude')]), fun = distHaversine)
      haversine_distance[i,j] = x
      haversine_distance[j,i] = x
    }
  }
  haversine_distance_df <- data.frame(haversine_distance)
  colnames(haversine_distance_df) <- as.character(bankUpdated_latlon$"Unique.ID")
  rownames(haversine_distance_df) <- as.character(bankUpdated_latlon$"Unique.ID")
  return(haversine_distance_df)
}


## calculate distance between accessions based on elevation from passport data
calculateElevationDistance = function(bankUpdated_latlon) {
  elev_diffs = matrix(nrow= nrow(bankUpdated_latlon), ncol= nrow(bankUpdated_latlon))
  for (i in 1:nrow(bankUpdated_latlon)){
    for (j in 1:i){
      elev_diffs[i, j] = abs(bankUpdated_latlon[i, 'Elevation'] - bankUpdated_latlon[j, 'Elevation'])
      elev_diffs[j, i] = elev_diffs[i, j]
    }
  }
  elev_diffs_df = as.data.frame(elev_diffs)
  colnames(elev_diffs_df) <- as.character(bankUpdated_latlon$"Unique.ID")
  rownames(elev_diffs_df) <- as.character(bankUpdated_latlon$"Unique.ID")
  
  return(elev_diffs_df)
}

calculateElevationDistanceSigned = function(bankUpdated_latlon) {
  elev_diffs = matrix(nrow= nrow(bankUpdated_latlon), ncol= nrow(bankUpdated_latlon))
  for (i in 1:nrow(bankUpdated_latlon)){
    for (j in 1:i){
      elev_diffs[i, j] = bankUpdated_latlon[i, 'Elevation'] - bankUpdated_latlon[j, 'Elevation']
      elev_diffs[j, i] = -elev_diffs[i, j]
    }
  }
  elev_diffs_df = as.data.frame(elev_diffs)
  colnames(elev_diffs_df) <- as.character(bankUpdated_latlon$"Unique.ID")
  rownames(elev_diffs_df) <- as.character(bankUpdated_latlon$"Unique.ID")
  
  return(elev_diffs_df)
}


## mmrr function from Wang et al, adapted to report statistics
MMRR<-function(Y,X,nperm=999){
  #compute regression coefficients and test statistics
  nrowsY<-nrow(Y)
  y<-unfold(Y)
  if(is.null(names(X)))names(X)<-paste("X",1:length(X),sep="")
  Xmats<-sapply(X,unfold)
  fit<-lm(y~Xmats)
  coeffs<-fit$coefficients
  summ<-summary(fit)
  r.squared<-summ$r.squared
  tstat<-summ$coefficients[,"t value"]
  Fstat<-summ$fstatistic[1]
  tprob<-rep(1,length(tstat))
  Fprob<-1
  anova_table = anova(fit)
  sse = anova_table[2,2]
  ssr = anova_table[1,2]
  
  #perform permutations
  for(i in 1:nperm){
    rand<-sample(1:nrowsY)
    Yperm<-Y[rand,rand]
    yperm<-unfold(Yperm)
    fit<-lm(yperm~Xmats)
    summ<-summary(fit)
    Fprob<-Fprob+as.numeric(summ$fstatistic[1]>=Fstat)
    tprob<-tprob+as.numeric(abs(summ$coefficients[,"t value"])>=abs(tstat))
  }
  
  #return values
  tp<-tprob/(nperm+1)
  Fp<-Fprob/(nperm+1)
  names(r.squared)<-"r.squared"
  names(coeffs)<-c("Intercept",names(X))
  names(tstat)<-paste(c("Intercept",names(X)),"(t)",sep="")
  names(tp)<-paste(c("Intercept",names(X)),"(p)",sep="")
  names(Fstat)<-"F-statistic"
  names(Fp)<-"F p-value"
  names(sse) <- 'SSE'
  names(ssr) <- 'SSR'
  return(list(r.squared=r.squared,
              coefficients=coeffs,
              tstatistic=tstat,
              tpvalue=tp,
              Fstatistic=Fstat,
              Fpvalue=Fp,
              SSE = sse,
              SSR = ssr))
}

# unfold converts the lower diagonal elements of a matrix into a vector
# unfold is called by MMRR
unfold<-function(X){
  x<-vector()
  for(i in 2:nrow(X)) x<-c(x,X[i,1:i-1])
  return(x)
}

# unfold converts the lower diagonal elements of a matrix into a vector
# unfold is called by MMRR
gaussian_unfold<-function(X, bandwidth){
  sd_X = sd(X)
  x<-vector()
  for(i in 2:nrow(X)) x<-c(x,X[i,1:i-1])
  return(scale(exp(-x/sd_X/bandwidth)))
}

scale_variance = function(x) {
  return((x - min(x))/(max(x) - min(x)))
}

kinship2dist<-function(x){
  diag(x)<-NA
  minx<-min(x,na.rm=TRUE)
  d.x<-1-(x-minx)/(1-minx)
  diag(d.x)<-0
  return(d.x)
}

## calculate levenshtein distances between languages instead of using MRCA or branch length
obtainLevenshteinDistances = function(language_list) {
  levenshtein = vroom('../data/levenshteinLanguageDistances.csv')
  top_leven = levenshtein %>% filter(language1 %in% language_list & language2 %in% language_list)
  bottom_leven = setNames(top_leven, c('language2', 'language1', 'L-distance'))[c(2, 1, 3)]
  sub_leven = rbind(top_leven, bottom_leven)
  sub_leven = pivot_wider(sub_leven, id_cols = 'language1', names_from = 'language2', values_from = 'L-distance', values_fill = 0) %>% column_to_rownames('language1')
  lang_distances = sub_leven[order(rownames(sub_leven)),order(colnames(sub_leven))]
  return(lang_distances)
}

## based on average distances between languages for overlapping languages or tips in a branch, create a language distance matrix between accessions
createAccessionDistanceList = function(MasterGWAS, language_matrix) {
  langlist = c()
  for(i in 1:length(MasterGWAS$`Unique.ID`)) {
    for (j in 1:i) {
      x = language_matrix[MasterGWAS[[i, 'LangName']],
                                  MasterGWAS[[j, 'LangName']]]
      langlist[[length(langlist) + 1]] = data.frame(A = MasterGWAS[[i, 'Unique.ID']], B = MasterGWAS[[j, 'Unique.ID']], distance = x)
      # langlist[[length(langlist) + 1]] = x
      # language_distance_df[i,j] = x
      # language_distance_df[j,i] = x
    }
  }
  return(langlist)
}

generateLanguageMatrix = function(family_language_distance_matrix, family_mrca_language_list, mrca_mapping_sheet) {
  ## given language distance matrix and list of created mrca language names, 
  ## plus mrca mapping sheet, generate a full language matrix
  
  mrca_means = lapply(family_mrca_language_list, function(lang) family_language_distance_matrix[, mrca_mapping_sheet %>% 
                                                                                                  filter(`MRCA_abbreviated` == lang) %>% 
                                                                                                  pull(`Tip_abbreviated`)] %>% rowMeans())
  
  left_side = rbind(family_language_distance_matrix, data.frame(do.call(rbind, mrca_means), row.names = family_mrca_language_list))
  mrca_matrix = makeMRCAMatrix(family_language_distance_matrix, lapply(family_mrca_language_list, function(lang) mrca_mapping_sheet %>% 
                                                                         filter(`MRCA_abbreviated` == lang) %>% 
                                                                         pull(`Tip_abbreviated`)))
  
  rownames(mrca_matrix) = family_mrca_language_list
  right_side = rbind(do.call(cbind, mrca_means), mrca_matrix)
  colnames(right_side) = family_mrca_language_list
  results = cbind(left_side, right_side)
}

