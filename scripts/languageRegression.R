library(data.table)
# library(PopGenReport)
library(dplyr)
library(ggplot2)
library(stringr)
# library(hierfstat)
library(tidyr)
library(tidyverse)
library(randomForest)
library(purrr)
library(vroom)

language_families = c('mayan', 'aztecan', 'otomanguean', 'mayanHaynie', 'aztecanHaynie', 'otomangueanHaynie')
languages = c('aztecanHaynie', 'aztecanNL', 'mayanNL', 'mayanHaynie', 'otomangueanHaynie', 'otomangueanNL')

#### regression using distance matrices using pca and pcoa #################
regression_pca_means = data.frame()
regression_pcoa_means = data.frame()
for(language in languages) {
  regression_results = vroom(sprintf('../results/regression_%s.csv', language))
  regression_pca = regression_results %>% mutate(geography_only = full.r.squared - no_geography.r.squared,
                                                 geography_only_adj = full.adj.r.squared - no_geography.adj.r.squared,
                                                 elevation_only = full.r.squared - no_elevation.r.squared,
                                                 elevation_only_adj = full.adj.r.squared - no_elevation.adj.r.squared,
                                                 admix_only = full.r.squared - no_admix.r.squared,
                                                 admix_only_adj = full.adj.r.squared - no_admix.adj.r.squared,
                                                 language_only = full.r.squared - no_language.r.squared,
                                                 language_only_adj = full.adj.r.squared - no_language.adj.r.squared) 
  regression_pcoa = regression_results %>% mutate(geography_only = full_pcoa.r.squared - no_geography_pcoa.r.squared,
                                                  geography_only_adj = full_pcoa.adj.r.squared - no_geography_pcoa.adj.r.squared,
                                                  elevation_only = full_pcoa.r.squared - no_elevation_pcoa.r.squared,
                                                  elevation_only_adj = full_pcoa.adj.r.squared - no_elevation_pcoa.adj.r.squared,
                                                  admix_only = full_pcoa.r.squared - no_admix_pcoa.r.squared,
                                                  admix_only_adj = full_pcoa.adj.r.squared - no_admix_pcoa.adj.r.squared,
                                                  language_only = full_pcoa.r.squared - no_language_pcoa.r.squared,
                                                  language_only_adj = full_pcoa.adj.r.squared - no_language._pcoaadj.r.squared)
  regression_pca_means = rbind(regression_pca_means, regression_pca %>% 
                                 dplyr::select(-c('snp')) %>% 
                                 dplyr::select(c(full.r.squared, full.adj.r.squared, 
                                                 geography_only, elevation_only, admix_only, language_only,
                                                 geography_only_adj, elevation_only_adj, 
                                                 admix_only_adj, language_only_adj)) %>% 
                                 colMeans())
  regression_pcoa_means = rbind(regression_pcoa_means, regression_pcoa %>% 
                                  dplyr::select(-c('snp')) %>% 
                                  dplyr::select(c(full_pcoa.r.squared, full_pcoa.adj.r.squared, 
                                                  geography_only, elevation_only, admix_only, language_only,
                                                  geography_only_adj, elevation_only_adj, 
                                                  admix_only_adj, language_only_adj)) %>% 
                                  colMeans())
}


colnames(regression_pca_means) = c('full', 'adj_full', 
                                   'geography', 'elevation', 'admix', 'language',
                                   'adj_geography', 'adj_elevation', 
                                   'adj_admix', 'adj_language')
colnames(regression_pcoa_means) = c('full', 'adj_full', 
                                    'geography', 'elevation', 'admix', 'language',
                                    'adj_geography', 'adj_elevation', 
                                    'adj_admix', 'adj_language')

regression_pca_means$language_set = c('aztecan', 'aztecan', 'mayan', 'mayan', 'otomanguean', 'otomanguean')
regression_pca_means$dataset = c('Haynie', 'NL', 'Haynie', 'NL', 'Haynie', 'NL')
regression_pcoa_means$language_set = c('aztecan', 'aztecan', 'mayan', 'mayan', 'otomanguean', 'otomanguean')
regression_pcoa_means$dataset= c('Haynie', 'NL', 'Haynie', 'NL', 'Haynie', 'NL')

regression_pca_means_pivot = regression_pca_means %>% pivot_longer(cols = -c(language_set, dataset),
                                                                   values_to = 'mean_r2',names_to = 'model')

ggplot(regression_pca_means_pivot, aes(x = model, y = mean_r2)) +
  facet_grid(rows = vars(dataset), cols = vars(language_set)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.2))


regression_pcoa_means_pivot = regression_pcoa_means %>% pivot_longer(cols = -c(language_set, dataset),
                                                                     values_to = 'mean_r2',names_to = 'model')

ggplot(regression_pcoa_means_pivot %>% filter(model %in% c('adj_language', 'adj_geography', 'adj_full', 'adj_elevation', 'adj_admix')), 
       aes(x = model, y = mean_r2 * 100)) +
  facet_grid(rows = vars(dataset), cols = vars(language_set)) +
  geom_col() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.2)) +
  coord_flip() +
  ggtitle('PCoA')

ggplot(regression_pca_means_pivot %>% filter(model %in% c('adj_language', 'adj_geography', 'adj_full', 'adj_elevation', 'adj_admix')), 
       aes(x = model, y = mean_r2 * 100)) +
  facet_grid(rows = vars(dataset), cols = vars(language_set)) +
  geom_col() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.2)) +
  coord_flip() +
  ggtitle('PCA')

ggplot(regression_pcoa_means_pivot %>% filter(model %in% c('language', 'geography', 'full', 'elevation', 'admix')), 
       aes(x = model, y = mean_r2 * 100)) +
  facet_grid(rows = vars(dataset), cols = vars(language_set)) +
  geom_col() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.2)) +
  coord_flip() +
  ggtitle('PCoA')

ggplot(regression_pca_means_pivot %>% filter(model %in% c('language', 'geography', 'full', 'elevation', 'admix')), 
       aes(x = model, y = mean_r2 * 100)) +
  facet_grid(rows = vars(dataset), cols = vars(language_set)) +
  geom_col() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.2)) +
  coord_flip() +
  ggtitle('PCA')

cbind(languages, regression_pcoa_means[, -c(11:12)] * 100)
cbind(languages, regression_pca_means[, -c(11:12)] * 100)


#### using point values instead of distance matrices ############################  

regression_point_means = data.frame()
languages = c('aztecanHaynie', 'aztecanNL', 'mayanNL', 'mayanHaynie', 'otomangueanHaynie', 'otomangueanNL')
for(language in languages) {
  regression_results = vroom(sprintf('../results/regression_%s_point.csv', language), show_col_types = F)
  regression_point_means = rbind(regression_point_means, regression_results %>% 
                                   dplyr::select(-c('snp')) %>% 
                                   colMeans())
  print(language)
  print(t.test(regression_results$language_point.adj.r.squared, 
               regression_results$geo_admix_point.adj.r.squared)$p.value)
  print(mean(regression_results$geo_admix_point.adj.r.squared) * 100)
  print(t.test(regression_results$language_point.adj.r.squared, 
               regression_results$geo_admix_point.adj.r.squared)$estimate)
  print(mean(regression_results$language_point.adj.r.squared) - mean(regression_point_means$geo_admix_point.adj.r.squared))
  
}

## read values for maize PVE from PCA #######
maize_PVE = data.frame()
languages_underscore = c('aztecan_Haynie', 'aztecan_NL', 'mayan_Haynie', 'mayan_NL', 'otomanguean_Haynie', 'otomanguean_NL')
for(language in languages_underscore) {
  this_PVE = vroom(sprintf('../results/%s_pca_eigens.txt', language))
  this_PVE = t(this_PVE[,-1])
  this_PVE =  as.data.frame(this_PVE)
  colnames(this_PVE) = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")
  this_PVE$language_set = str_to_title(str_split_1(language, '_')[1])
  this_PVE$dataset = str_split_1(language, '_')[2]
  maize_PVE = rbind(maize_PVE, this_PVE)
  
}

colnames(regression_point_means) = colnames(regression_results)[-c(1)]
## labeling of correct languages and datasets
regression_point_means$language_set = c('aztecan', 'aztecan', 'mayan', 'mayan', 'otomanguean', 'otomanguean')
regression_point_means$dataset = c('Haynie', 'NL', 'Haynie', 'NL', 'Haynie', 'NL')

## calculate difference
regression_point_means$language_difference.point = regression_point_means$language_point.adj.r.squared - regression_point_means$geo_admix_point.adj.r.squared
regression_point_main_results = regression_point_means %>% dplyr::select(language_set, dataset, geo_admix_point.adj.r.squared, language_point.adj.r.squared, 
                                                                         language_difference.point)
# write.table(regression_point_main_results, '../results/regression_point_main_results', sep = '\t', quote = F)


regression_point_means_pivot = regression_point_means %>% 
  dplyr::select(ends_with('point.adj.r.squared'), 
                ends_with('pcoa.adj.r.squared'), 
                language_set, dataset) %>% 
  pivot_longer(cols = -c(language_set, dataset),
               values_to = 'mean_r2',names_to = 'model') %>% 
  mutate(percent_r2 = mean_r2 * 100) %>% 
  separate(col = 'model', into = c('model', 'predictor_type'), sep = '_p', remove = T) %>% 
  mutate(predictor_type = recode(predictor_type, oint.adj.r.squared = 'point values',
                                 coa.adj.r.squared = 'distance matrices'),
         language_set = recode(language_set, aztecan = 'Aztecan', mayan = 'Mayan', otomanguean = 'Otomanguean'),
         model = recode(model, language = 'LLEA + language', geo_admix = 'LLEA'))

(regression_variance_plot = ggplot(regression_point_means_pivot, aes(x = model, y = percent_r2, fill = predictor_type)) +
    facet_grid(rows = vars(dataset), cols = vars(language_set)) +
    geom_col(position = 'dodge2') +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.2)) +
    coord_flip() + 
    ylab('Genetic variation explained by model (%)') + xlab('')
)


(regression_variance_plot = ggplot(regression_point_means_pivot %>% filter(predictor_type == 'point values'), 
                                   aes(x = language_set, y = percent_r2, fill = model)) +
    facet_grid(rows = vars(dataset)) +
    geom_bar(position = 'dodge2', stat = 'identity') +
    theme_bw() +
    theme(legend.position = 'top', 
          legend.margin = unit(c(.1, .1, .1, .1), 'cm'),
          legend.box.margin = unit(c(.1, .1, .1, .1), 'cm'),
          plot.margin = unit(c(.3,.3,.3,.3), 'cm'),
    ) +
    coord_flip() + 
    scale_fill_discrete(name = '') +
    ylab('Genetic variation explained by model (%)') + xlab('')
)

ggsave('regression_variance_plot.png', 
       plot = regression_variance_plot,
       device = 'png', 
       path = '../plots/', 
       width = 8, height = 5, units = 'in')
