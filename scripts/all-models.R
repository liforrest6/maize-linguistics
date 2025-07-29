library(data.table)
# library(PopGenReport)
library(dplyr)
library(ggplot2)
library(stringr)
# library(hierfstat)
library(tidyr)
library(tidyverse)
library(randomForest)

## for parallelization
library(foreach)
library(doParallel)
registerDoParallel(6)

run = commandArgs(t=T)
model = run[1]
print(model)

source('languageFunctions.R')


K = as.matrix(fread('../data/K_allChr.csv', header = T), rownames = 1)
language_families = c('mayan', 'aztecan', 'otomanguean', 'mayanHaynie', 'aztecanHaynie', 'otomangueanHaynie')


nperm = 200

if(model == 'levenshtein') {
#### for levenshtein distances #########################################################################################

  all_variances_levenshtein = foreach(language = language_families,.combine = bind_rows, .errorhandling = 'stop') %dopar% {
    print(language)

    if (grepl('Haynie', language)) {
      elevation_mat = read.table(sprintf('../results/Haynie/%s_elevationSigned_distances.csv', language), sep = ',', row.names = 1, header = T)
      admix_mat = read.table(sprintf('../results/Haynie/%s_admixProp_distances.csv', language), row.names = 1, sep = ',', header = T)
      geographic_mat = read.table(sprintf('../results/Haynie/%s_haversine_distances.csv', language), row.names = 1, sep = ',', header = T)
      language_mat = read.table(sprintf('../results/Haynie/%s_levenshtein_language_distances.csv', language), row.names = 1, sep = ',', header = T)
      # branch_mat = read.table(sprintf('../results/Haynie/%s_branch_language_distances.csv', language), row.names = 1, sep = ',', header = T)
      # climate_mat = read.table(sprintf('../results/Haynie/%s_climate_distances.csv', language), row.names = 1, sep = ',', header = T)
    } else {
      elevation_mat = fread(sprintf('../results/NL/%s_elevationSigned_distances.csv', language)) %>% remove_rownames %>% column_to_rownames(var="V1")
      admix_mat = fread(sprintf('../results/NL/%s_admixProp_distances.csv', language)) %>% remove_rownames %>% column_to_rownames(var="V1")
      geographic_mat = fread(sprintf('../results/NL/%s_haversine_distances.csv', language)) %>% remove_rownames %>% column_to_rownames(var="V1")
      language_mat = fread(sprintf('../results/NL/%s_levenshtein_language_distances.csv', language)) %>% remove_rownames %>% column_to_rownames(var="V1")
      # branch_mat = fread(sprintf('../results/NL/%s_branch_language_distances.csv', language)) %>% remove_rownames %>% column_to_rownames(var="V1")
      # climate_mat = fread(sprintf('../results/NL/%s_climate_distances.csv', language)) %>% remove_rownames %>% column_to_rownames(var="V1")
    }
    
    language_mat = as.matrix(language_mat)
    elevation_mat = abs(as.matrix(elevation_mat))
    geographic_mat = as.matrix(geographic_mat)
    admix_mat = abs(as.matrix(admix_mat))

    print('imported matrices')
    
    K2 = K[rownames(language_mat), rownames(language_mat)]
    K_distance = kinship2dist(K2)
    genetic_mat = as.matrix(K_distance)
    print('obtained genetic mat')
    
    ## do MMRR on all variables
    total_MMRR = MMRR(genetic_mat, list(elevation_mat, geographic_mat, language_mat, admix_mat), nperm = nperm)
    total_variance = total_MMRR$r.squared
    print('done with total_MMRR')

    ## do MMRR on three variables
    language_MMRR = MMRR(genetic_mat, list(geographic_mat, elevation_mat, admix_mat), nperm = 5)
    elevation_MMRR = MMRR(genetic_mat, list(geographic_mat, language_mat, admix_mat), nperm = 5)
    geographic_MMRR = MMRR(genetic_mat, list(elevation_mat, language_mat, admix_mat), nperm = 5)
    admix_MMRR = MMRR(genetic_mat, list(geographic_mat, elevation_mat, language_mat), nperm = 5)
    print('done with partial MMRR')

    ## calculate partial R2 for each variable based on SSE_reduced - SSE_total / SSE_reduced
    language_variance = (language_MMRR$SSE - total_MMRR$SSE) / language_MMRR$SSE
    elevation_variance = (elevation_MMRR$SSE - total_MMRR$SSE) / elevation_MMRR$SSE
    geographic_variance = (geographic_MMRR$SSE - total_MMRR$SSE) / geographic_MMRR$SSE
    admix_variance = (admix_MMRR$SSE - total_MMRR$SSE) / admix_MMRR$SSE
    # total_variance = MMRR(genetic_mat, list(elevation_mat, geographic_mat, language_mat, admix_mat), nperm = nperm)$r.squared
    # geographic_variance = total_variance - MMRR(genetic_mat, list(elevation_mat, language_mat, admix_mat), nperm = nperm)$r.squared
    # elevation_variance = total_variance - MMRR(genetic_mat, list(geographic_mat, language_mat, admix_mat), nperm = nperm)$r.squared
    # language_variance = total_variance - MMRR(genetic_mat, list(geographic_mat, elevation_mat, admix_mat), nperm = nperm)$r.squared
    # admix_variance = total_variance - MMRR(genetic_mat, list(geographic_mat, elevation_mat, language_mat), nperm = nperm)$r.squared
    
    language_summary = data.frame(language_family = language,
                                  total_variance = total_variance, 
                                  geographic_variance = geographic_variance, 
                                  elevation_variance = elevation_variance,
                                  language_variance = language_variance,
                                  admix_variance = admix_variance,
                                  Fstat = total_MMRR$Fstatistic,
                                  Fpvalue = total_MMRR$Fpvalue) 
    
    # all_variances_levenshtein = rbind(all_variances_levenshtein, language_summary)
    language_summary
  }


  ### write results to csv files
  write.csv(all_variances_levenshtein, '../results/all_variances_levenshtein.csv', quote = F)
  print('written csv')

  ### generate plots

  all_variances_levenshtein_long = pivot_longer(all_variances_levenshtein, 
                                                cols = c('elevation_variance', 'language_variance', 'geographic_variance', 
                                                         'admix_variance'), 
                                                names_to = 'distance_matrix',
                                                values_to = 'variance') %>% 
    mutate(variable = str_split_i(distance_matrix, '_', 1),
           variance_percent = variance / total_variance)

  ### color palette
  variance_palette = c("#D1235D","#9ECE27","#25B2F8","#bce7fd","#c492b1")

  levenshtein_variance_plot = ggplot(all_variances_levenshtein_long, aes(y = language_family, fill = variable, x = variance)) +
    geom_col(position = position_stack(reverse = T)) + 
    scale_fill_manual(values = variance_palette) +
    ylab('Language family') + 
    xlab('Proportion of total genetic variance explained') +
    ggtitle('Levenshtein distance')

  ggsave('levenshtein_variance.png',
       plot = levenshtein_variance_plot,
       device = 'png',
       '../plots/',
       width = 6,
       height = 4,
       units = 'in')

} else if (model == 'branch') {


#### for average branch distances #########################################################################################

  all_variances_average_branch = foreach(language = language_families,.combine = bind_rows, .errorhandling = 'stop') %dopar% {  
    
    print(language)
    if (grepl('Haynie', language)) {
      elevation_mat = read.table(sprintf('../results/Haynie/%s_elevationSigned_distances.csv', language), sep = ',', row.names = 1, header = T)
      admix_mat = read.table(sprintf('../results/Haynie/%s_admixProp_distances.csv', language), row.names = 1, sep = ',', header = T)
      geographic_mat = read.table(sprintf('../results/Haynie/%s_haversine_distances.csv', language), row.names = 1, sep = ',', header = T)
      language_mat = read.table(sprintf('../results/Haynie/%s_branch_language_distances.csv', language), row.names = 1, sep = ',', header = T)
      # climate_mat = read.table(sprintf('../results/Haynie/%s_climate_distances.csv', language), row.names = 1, sep = ',', header = T)
    } else {
      elevation_mat = fread(sprintf('../results/NL/%s_elevationSigned_distances.csv', language)) %>% remove_rownames %>% column_to_rownames(var="V1")
      admix_mat = fread(sprintf('../results/NL/%s_admixProp_distances.csv', language)) %>% remove_rownames %>% column_to_rownames(var="V1")
      geographic_mat = fread(sprintf('../results/NL/%s_haversine_distances.csv', language)) %>% remove_rownames %>% column_to_rownames(var="V1")
      language_mat = fread(sprintf('../results/NL/%s_branch_language_distances.csv', language)) %>% remove_rownames %>% column_to_rownames(var="V1")
      # climate_mat = fread(sprintf('../results/NL/%s_climate_distances.csv', language)) %>% remove_rownames %>% column_to_rownames(var="V1")
    }
    
    language_mat = as.matrix(language_mat)
    elevation_mat = abs(as.matrix(elevation_mat))
    geographic_mat = as.matrix(geographic_mat)
    admix_mat = abs(as.matrix(admix_mat))
    print('imported matrices')
    
    K2 = K[rownames(language_mat), rownames(language_mat)]
    K_distance = kinship2dist(K2)
    genetic_mat = as.matrix(K_distance)
    
    print('obtained genetic mat')

    total_MMRR = MMRR(genetic_mat, list(elevation_mat, geographic_mat, language_mat, admix_mat), nperm = nperm)
    total_variance = total_MMRR$r.squared
      print('done with total_MMRR')

    language_MMRR = MMRR(genetic_mat, list(geographic_mat, elevation_mat, admix_mat), nperm = 5)
    elevation_MMRR = MMRR(genetic_mat, list(geographic_mat, language_mat, admix_mat), nperm = 5)
    geographic_MMRR = MMRR(genetic_mat, list(elevation_mat, language_mat, admix_mat), nperm = 5)
    admix_MMRR = MMRR(genetic_mat, list(geographic_mat, elevation_mat, language_mat), nperm = 5)
    print('done with partial MMRR')

    language_variance = (language_MMRR$SSE - total_MMRR$SSE) / language_MMRR$SSE
    elevation_variance = (elevation_MMRR$SSE - total_MMRR$SSE) / elevation_MMRR$SSE
    geographic_variance = (geographic_MMRR$SSE - total_MMRR$SSE) / geographic_MMRR$SSE
    admix_variance = (admix_MMRR$SSE - total_MMRR$SSE) / admix_MMRR$SSE
    # total_variance = MMRR(genetic_mat, list(elevation_mat, geographic_mat, language_mat, admix_mat), nperm = nperm)$r.squared
    # geographic_variance = total_variance - MMRR(genetic_mat, list(elevation_mat, language_mat, admix_mat), nperm = nperm)$r.squared
    # elevation_variance = total_variance - MMRR(genetic_mat, list(geographic_mat, language_mat, admix_mat), nperm = nperm)$r.squared
    # language_variance = total_variance - MMRR(genetic_mat, list(geographic_mat, elevation_mat, admix_mat), nperm = nperm)$r.squared
    # admix_variance = total_variance - MMRR(genetic_mat, list(geographic_mat, elevation_mat, language_mat), nperm = nperm)$r.squared
    
    language_summary = data.frame(language_family = language,
                                  total_variance = total_variance, 
                                  geographic_variance = geographic_variance, 
                                  elevation_variance = elevation_variance,
                                  language_variance = language_variance,
                                  admix_variance = admix_variance,
                                  Fstat = total_MMRR$Fstatistic,
                                  Fpvalue = total_MMRR$Fpvalue) 
    
    # all_variances_levenshtein = rbind(all_variances_levenshtein, language_summary)
    language_summary
    
  }

  ### write results to csv files


  write.csv(all_variances_average_branch, '../results/all_variances_branch_length.csv', quote = F)
  print('written csv')

  ### generate plots


  all_variances_average_branch_long = pivot_longer(all_variances_average_branch, 
                                                   cols = c('elevation_variance', 'language_variance', 'geographic_variance', 
                                                            'admix_variance'), 
                                                   names_to = 'distance_matrix',
                                                   values_to = 'variance') %>% 
    mutate(variable = str_split_i(distance_matrix, '_', 1),
           variance_percent = variance / total_variance)

  ### color palette
  variance_palette = c("#D1235D","#9ECE27","#25B2F8","#bce7fd","#c492b1")

  average_branch_variance_plot = ggplot(all_variances_average_branch_long, aes(y = language_family, fill = variable, x = variance)) +
    geom_col(position = position_stack(reverse = T)) + 
    scale_fill_manual(values = variance_palette) +
    ylab('Language family') + 
    xlab('Proportion of total genetic variance explained') +
    ggtitle('Average branch length distance')


  ggsave('average_branch_variance.png',
         plot = average_branch_variance_plot,
         device = 'png',
         '../plots/',
         width = 6,
         height = 4,
         units = 'in')

} else if (model == 'RF') {

######### testing feature importance using random forests #########################################################################################


  all_node_purity = foreach(language = language_families,.combine = bind_rows, .errorhandling = 'stop') %dopar% {  
    print(language)
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
    
    K2 = K[rownames(levenshtein_mat), rownames(levenshtein_mat)]
    K_distance = kinship2dist(K2)
    genetic_mat = as.matrix(K_distance)

    print('imported matrices')
    
    
    variables = data.frame(genetic = genetic_mat[upper.tri(genetic_mat)],
                           levenshtein = levenshtein_mat[upper.tri(levenshtein_mat)],
                           branch = branch_mat[upper.tri(branch_mat)],
                           geography = geography_mat[upper.tri(geography_mat)],
                           elevation = abs(elevation_mat[upper.tri(elevation_mat)]),
                           admixture = abs(admix_mat[upper.tri(admix_mat)]),
                           climate = climate_mat[upper.tri(climate_mat)]
    )
    
    sample_size = floor(0.75 * nrow(variables))
    train_ind = sample(seq_len(nrow(variables)), size = sample_size)
    data_train = variables[train_ind,]
    data_test = variables[-train_ind,]
    
    rf_model_variables = randomForest(x = data_train %>% dplyr::select(levenshtein, geography, elevation, admixture), 
                                      y = data_train %>% dplyr::pull(genetic),
                                      ntree=100,
                                      importance=TRUE)

    print('finished training')

    # rf_model_predictions = predict(rf_model_variables, newdata = data_test %>% dplyr::select(levenshtein, geography, elevation, admixture, climate))
    # rf_model_variables_results = cor(data_test %>% pull(genetic), rf_model_predictions)
    
    this_importance = importance(rf_model_variables)
    t(this_importance)[2,]

    # print(all_node_purity)
  }

  write.csv(all_node_purity, '../results/randomForest_nodepurity.csv', quote = F)

}


print('done')