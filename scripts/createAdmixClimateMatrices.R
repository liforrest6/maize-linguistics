
## prepare admixture proportion table by removing sample metadata from IDs

admixProp = read.table('../data/admixProp.GBSsamples.txt')
admixProp$Unique.ID = gsub('.MRG.4.', ':', admixProp$V1)
admixProp$Unique.ID = gsub('.D17PEACXX.3.', ':', admixProp$Unique.ID)


## read language matrix given language and method
readLanguageMatrix = function(language, method) {
  language_matrix = read.csv(sprintf('../results/%s_%s_language_distances.csv', language, method))
  return(language_matrix)
}

## calculate admixture distance
calculateAdmixDist = function(language_matrix, admixProp) {
  languageNames = language_matrix[,1]
  admixProp_sorted = admixProp[match(languageNames, admixProp$Unique.ID), ]
  result = as.data.frame(outer(c(1:length(languageNames)), c(1:length(languageNames)), 
                               FUN = function(a, b, x) x[a] - x[b],
                               x = admixProp_sorted$V2),
                         row.names = languageNames)
  colnames(result) = languageNames
  return(result)
}

for (language in c('aztecan', 'mayan', 'otomanguean', 'aztecanHaynie', 'mayanHaynie', 'otomangueanHaynie')) {
  print(language)
  admix_dist_matrix = calculateAdmixDist(readLanguageMatrix(language, 'levenshtein'), admixProp)
  if (grepl('Haynie', language)) {
    write.csv(admix_dist_matrix, sprintf('../results/Haynie/%s_admixProp_distances.csv', language), row.names = T, quote = F)
  } else 
    write.csv(admix_dist_matrix, sprintf('../results/NL/%s_admixProp_distances.csv', language), row.names = T, quote = F)
}

climate_vals = read.csv('../data/GEA-climate-invnormtransformed.csv')

calculateClimateDist = function(language_matrix, climateValues) {
  languageNames = language_matrix[,1]
  climateValues_sorted = climateValues[match(languageNames, climateValues$Unique.ID), ]
  result = as.data.frame(outer(c(1:length(languageNames)), c(1:length(languageNames)), 
                               FUN = function(a, b, x, y) {
                                 sqrt((x[a] - x[b])^2 + (y[a] - y[b])^2)
                                 },
                               x = climateValues_sorted$tmax,
                               y = climateValues_sorted$precipTot),
                         row.names = languageNames)
  colnames(result) = languageNames
  return(result)
}
 

for (language in c('aztecan', 'mayan', 'otomanguean', 'aztecanHaynie', 'mayanHaynie', 'otomangueanHaynie')) {
  print(language)
  climate_dist_matrix = calculateClimateDist(readLanguageMatrix(language, 'levenshtein'), climate_vals)
  if (grepl('Haynie', language)) {
    write.csv(climate_dist_matrix, sprintf('../results/Haynie/%s_climate_distances.csv', language), row.names = T, quote = F)
  }
  else
    write.csv(climate_dist_matrix, sprintf('../results/NL/%s_climate_distances.csv', language), row.names = T, quote = F)
}


elevation_test = read.table('../results/NL/otomanguean_elevationSigned_distances.csv', sep = ',', row.names = 1, header = T)
admix_test = read.table('../results/otomanguean_admixProp_distances.csv', row.names = 1, sep = ',', header = T)
geography_test = read.table('../results/NL/otomanguean_haversine_distances.csv', row.names = 1, sep = ',', header = T)
levenshtein_test = read.table('../results/NL/otomanguean_levenshtein_language_distances.csv', row.names = 1, sep = ',', header = T)
branch_test = read.table('../results/NL/otomanguean_branch_language_distances.csv', row.names = 1, sep = ',', header = T)
climate_test = read.table('../results/NL/otomanguean_climate_distances.csv', row.names = 1, sep = ',', header = T)

cor(elevation_test[upper.tri(elevation_test)], admix_test[upper.tri(admix_test)])
cor(geography_test[upper.tri(geography_test)], admix_test[upper.tri(admix_test)])
cor(geography_test[upper.tri(geography_test)], elevation_test[upper.tri(elevation_test)])
cor(levenshtein_test[upper.tri(levenshtein_test)], branch_test[upper.tri(branch_test)])
cor(levenshtein_test[upper.tri(levenshtein_test)], elevation_test[upper.tri(elevation_test)])
cor(branch_test[upper.tri(levenshtein_test)], elevation_test[upper.tri(elevation_test)])
cor(geography_test[upper.tri(geography_test)], levenshtein_test[upper.tri(elevation_test)])
cor(climate_test[upper.tri(climate_test)], levenshtein_test[upper.tri(elevation_test)])
cor(climate_test[upper.tri(climate_test)], geography_test[upper.tri(geography_test)])
cor(climate_test[upper.tri(climate_test)], elevation_test[upper.tri(elevation_test)])

library(GGally)

#  make ggpairs correlation matrices between all variable distance matrices
for(language in c('aztecan', 'mayan', 'otomanguean', 'aztecanHaynie', 'mayanHaynie', 'otomangueanHaynie')) {
  if (grepl('Haynie', language)) {
    elevation_mat = read.table(sprintf('../results/Haynie/%s_elevationSigned_distances.csv', language), sep = ',', row.names = 1, header = T)
    admix_mat = read.table(sprintf('../results/Haynie/%s_admixProp_distances.csv', language), row.names = 1, sep = ',', header = T)
    geography_mat = read.table(sprintf('../results/Haynie/%s_haversine_distances.csv', language), row.names = 1, sep = ',', header = T)
    levenshtein_mat = read.table(sprintf('../results/Haynie/%s_levenshtein_language_distances.csv', language), row.names = 1, sep = ',', header = T)
    branch_mat = read.table(sprintf('../results/Haynie/%s_branch_language_distances.csv', language), row.names = 1, sep = ',', header = T)
    climate_mat = read.table(sprintf('../results/Haynie/%s_climate_distances.csv', language), row.names = 1, sep = ',', header = T)
  } else {
    elevation_mat = read.table(sprintf('../results/NL/%s_elevationSigned_distances.csv', language), sep = ',', row.names = 1, header = T)
    admix_mat = read.table(sprintf('../results/NL/%s_admixProp_distances.csv', language), row.names = 1, sep = ',', header = T)
    geography_mat = read.table(sprintf('../results/NL/%s_haversine_distances.csv', language), row.names = 1, sep = ',', header = T)
    levenshtein_mat = read.table(sprintf('../results/NL/%s_levenshtein_language_distances.csv', language), row.names = 1, sep = ',', header = T)
    branch_mat = read.table(sprintf('../results/NL/%s_branch_language_distances.csv', language), row.names = 1, sep = ',', header = T)
    climate_mat = read.table(sprintf('../results/NL/%s_climate_distances.csv', language), row.names = 1, sep = ',', header = T)
  }
  

  
  variables = data.frame(levenshtein = levenshtein_mat[upper.tri(levenshtein_mat)],
                         branch = branch_mat[upper.tri(branch_mat)],
                         geography = geography_mat[upper.tri(geography_mat)],
                         elevation = elevation_mat[upper.tri(elevation_mat)],
                         admixture = admix_mat[upper.tri(admix_mat)],
                         climate = climate_mat[upper.tri(climate_mat)]
                         )
  
  variable_correlations = ggpairs(variables, lower = list(continuous = 'density'),
                                  upper = list(continuous = wrap('cor', size = 2))) 
  ggsave(sprintf('variable_correlations_%s.png', language),
         plot = variable_correlations,
         device = 'png',
         '../plots/',
         width = 5,
         height = 5,
         units = 'in')
  
}
