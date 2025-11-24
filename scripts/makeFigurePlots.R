## figures for manuscript
library(grid)
library(gridExtra)
library(png)
library(cowplot)
library(ggplotify)
library(ggmap)
library(dplyr)
library(vroom)

source('languageFunctions.R')
source('apikey.R')
#### manuscript figure for all languages across all mesoamerica ################

### make map of all language families together ##############
mesoamerica_terrain_kernels = get_stadiamap(c(left = -114, bottom = 8, right = -82, top = 33),
                                            zoom = 6, maptype = 'stamen_toner_lite', color = 'bw', messaging = F)
mesoamerica_terrain_kernels_transformed = ggmap_bbox(mesoamerica_terrain_kernels)

all_haynie_accessions = rbind(otomanguean_haynie_accessions, aztecan_haynie_accessions, mayan_haynie_accessions)
all_NL_accessions = rbind(otomanguean_NL_accessions, aztecan_NL_accessions, mayan_NL_accessions)

(all_haynie_polygon_map = ggmap(mesoamerica_terrain_kernels_transformed) +
    coord_sf(crs = st_crs(3857)) +
    geom_sf(aes(fill = family_nam),
            data = st_transform(haynie_lang %>% 
                                  filter(Name %in% c(otomanguean_haynie_mapping$`Polygon Name`,
                                                     aztecan_haynie_mapping$`Polygon Name`,
                                                     mayan_haynie_mapping$`Polygon Name`)), 
                                3857),
            inherit.aes = F,
            alpha = 0.4) +
    geom_sf(data = st_transform(rbind(otomanguean_haynie_accessions, 
                                      aztecan_haynie_accessions,
                                      mayan_haynie_accessions), 3857), 
            inherit.aes = F,
            size = 0.1,
            alpha = 0.3) +
    scale_fill_discrete(name = '') +
    theme(legend.position = 'bottom',
          plot.margin = unit(c(.1,.1,.1,.1), 'cm')) +
    ggtitle('') +
    xlab('') +
    ylab('')
)

ggsave('all_haynie_polygon_map.png',
       plot = all_haynie_polygon_map, 
              device = 'png', 
              path = '../plots/', 
              width = 8, height = 6, units = 'in')

name_in_nativelands = function(x) {
  result = ''
  if(x %in% otomanguean_nativelands_mapping$`Polygon Name`)
    result = 'Otomanguean'
  else if(x %in% aztecan_nativelands_mapping$`Polygon Name`)
    result = 'Uto-Aztecan'
  else if(x %in% mayan_nativelands_mapping$`Polygon Name`)
    result = 'Mayan'
  return(result)
}

nativelands_lang$family_name = unlist(lapply(nativelands_lang$Name, name_in_nativelands))

(all_nativelands_polygon_map = ggmap(mesoamerica_terrain_kernels_transformed) +
    coord_sf(crs = st_crs(3857)) +
    geom_sf(aes(fill = family_name),
            data = st_transform(nativelands_lang %>% 
                                  filter(Name %in% c(otomanguean_nativelands_mapping$`Polygon Name`,
                                                     aztecan_nativelands_mapping$`Polygon Name`,
                                                     mayan_nativelands_mapping$`Polygon Name`)), 
                                3857),
            inherit.aes = F,
            alpha = 0.4) +
    geom_sf(data = st_transform(rbind(otomanguean_NL_accessions, 
                                      aztecan_NL_accessions,
                                      mayan_NL_accessions), 3857), 
            inherit.aes = F,
            size = 0.4,
            alpha = 0.5) +
    theme(legend.position = "bottom",
          plot.margin = unit(c(.1,.1,.1,.1), 'cm')) +
    scale_fill_discrete(name = 'language family') +
    ggtitle('') +
    xlab('') +
    ylab('')
)

ggsave('all_nativelands_polygon_map.png',
       plot = all_nativelands_polygon_map, 
       device = 'png', 
       path = '../plots/', 
       width = 8, height = 6, units = 'in')


#### make example language phylogenetic tree ############

mayan_haynie_tree_figure = mayan_haynie_tree_filtered
mayan_haynie_clades = data.frame(node = c(126, 161, 163, 166, 165, 107, 102, 98, 88, 101),
                                 clade = c('KICHE','CAKCHIQUEL', 'TZUTUJIL', 'KEKCHI', 'POQOM', 'MAM', 'TZOTZIL', 'CHOL', 'YUCATAN MAYAN', 'TZELTAL'))
{mayan_haynie_tree_figure$tip.label[77] = 'SOUTHERN_CAKCHIQUEL'
mayan_haynie_tree_figure$tip.label[19] = 'IXIL_CHAJUL'
mayan_haynie_tree_figure$tip.label[8] = 'CHONTAL_TABASCO'
mayan_haynie_tree_figure$tip.label[12] = 'CHORTI'
mayan_haynie_tree_figure$tip.label[11] = 'CHOLTI'
mayan_haynie_tree_figure$tip.label[7] = 'CHUJ'
mayan_haynie_tree_figure$tip.label[6] = 'TOJOLABAL'
mayan_haynie_tree_figure$tip.label[5] = 'QANJOBAL_SANTA_EULALIA'
mayan_haynie_tree_figure$tip.label[84] = 'HUASTEC_VERACRUZ'
mayan_haynie_tree_figure$tip.label[83] = 'CHICOMUCELTEC'
mayan_haynie_tree_figure$tip.label[85] = 'HUASTEC'}

mayan_haynie_tree_figure_transform = mayan_haynie_tree_figure 
# mayan_haynie_tree_figure_transform$edge.length[c(1,2)] = mayan_haynie_tree_figure$edge.length[c(1,2)] - 0.05
mayan_haynie_tree_figure_transform$edge.length[mayan_haynie_tree_figure_transform$edge.length > 0.12] = mayan_haynie_tree_figure_transform$edge.length[mayan_haynie_tree_figure_transform$edge.length > 0.12] - 0.05

(mayan_tree_figure = ggtree(mayan_haynie_tree_figure_transform) %>% 
    scaleClade(126, .07) %>% # scale kiche
    scaleClade(107, .1) %>% # scale mam
    scaleClade(102, .5) %>% # scale tzotzil
    scaleClade(88, .6) %>% # scale yucatan
    scaleClade(161, .6) %>% # scale cakchiquel
    scaleClade(163, .6) %>% # scale tzutujil
    scaleClade(98, .6) %>% # scale chol
    collapse(126, 'mixed', fill="darkgreen") %>% # kiche
    collapse(161, 'mixed', fill="salmon") %>%  #cakchiquel
    collapse(163, 'mixed', fill="magenta")  %>% #tzutujil
    collapse(166, 'mixed', fill="lightgreen") %>% #kekchi
    collapse(165, 'mixed', fill="skyblue") %>% #poqom
    collapse(107, 'mixed', fill="turquoise") %>% #mam
    collapse(102, 'mixed', fill = 'purple') %>% #tzotzil
    collapse(98, 'mixed', fill = 'orange') %>% #chol
    collapse(88, 'mixed', fill = 'red') %>% #yucatan, then tzeltal
    collapse(101, 'mixed', fill = 'blue') %<+% haynie_mapping +
    geom_tiplab(
      color = 'black',
      geom = 'label', 
      label.padding = unit(0.15, 'lines'),
      label.size = 0,
      size = 1.8) + 
    geom_cladelab(data = mayan_haynie_clades, 
                  mapping = aes(node = node, 
                                label = clade),
                  geom = 'text', 
                  align = F, fontsize = 1.8, offset = 0.05) + 
    xlim(0,3) +
    labs(fill = 'sub_grouping',
         title = 'Mayan Haynie') +
    theme(plot.margin = unit(c(.1,.1,.1,.1), 'cm')) +
    ggtitle('')
)

## grayscale tree
(mayan_tree_figure = ggtree(mayan_haynie_tree_figure_transform) %>% 
    scaleClade(126, .07) %>% # scale kiche
    scaleClade(107, .1) %>% # scale mam
    scaleClade(102, .5) %>% # scale tzotzil
    scaleClade(88, .6) %>% # scale yucatan
    scaleClade(161, .6) %>% # scale cakchiquel
    scaleClade(163, .6) %>% # scale tzutujil
    scaleClade(98, .6) %>% # scale chol
    collapse(126, 'mixed', fill="darkgray") %>% # kiche
    collapse(161, 'mixed', fill="darkgray") %>%  #cakchiquel
    collapse(163, 'mixed', fill="darkgray")  %>% #tzutujil
    collapse(166, 'mixed', fill="darkgray") %>% #kekchi
    collapse(165, 'mixed', fill="darkgray") %>% #poqom
    collapse(107, 'mixed', fill="darkgray") %>% #mam
    collapse(102, 'mixed', fill = 'darkgray') %>% #tzotzil
    collapse(98, 'mixed', fill = 'darkgray') %>% #chol
    collapse(88, 'mixed', fill = 'darkgray') %>% #yucatan, then tzeltal
    collapse(101, 'mixed', fill = 'darkgray') %<+% haynie_mapping +
    geom_tiplab(
      color = 'black',
      geom = 'label', 
      label.padding = unit(0.15, 'lines'),
      label.size = 0,
      size = 1.8) + 
    geom_cladelab(data = mayan_haynie_clades, 
                  mapping = aes(node = node, 
                                label = clade),
                  geom = 'text', 
                  align = F, fontsize = 1.8, offset = 0.05) + 
    xlim(0,3) +
    labs(fill = 'sub_grouping',
         title = 'Mayan Haynie') +
    theme(plot.margin = unit(c(.1,.1,.1,.1), 'cm')) +
    ggtitle('')
)

geom_text(aes(label = node), hjust = 2, size = 2.5) ## node numbers

#### stacked barplot for regression ################
fig3_PVE = regression_point_means_pivot %>% 
  filter(predictor_type == 'point values') %>% 
  pivot_wider(names_from = model, values_from = c('mean_r2', 'percent_r2')) %>% 
  mutate(language_difference = `percent_r2_LLEA + language` - percent_r2_LLEA,
         analysis = 'models') %>% 
  pivot_longer(cols = c('percent_r2_LLEA', 'language_difference'), names_to = 'model', values_to = 'percent') %>% 
  dplyr::select(language_set, dataset, model, percent, analysis)

PVE_maize10PCs = maize_PVE %>% 
  pivot_longer(cols = starts_with('PC'), values_to = 'percent', names_to = 'model') %>% 
  mutate(percent = percent * 100,
         analysis = 'top 10 maize PCs')

fig3_PVE_maize10PCs = rbind(fig3_PVE, PVE_maize10PCs)

ggplot(fig3_PVE_maize10PCs %>% filter(dataset == 'Haynie'), aes(x = c('models', 'top 10 maize PCs'))) +
  geom_bar(data = fig3_PVE_maize10PCs %>% filter(model %in% c('percent_r2_LLEA', 'language_difference')) %>% 
             filter(dataset == 'Haynie'),
           stat = 'identity', position = 'stack',
           aes(y = percent, x = 'models', fill = model)) +
  scale_fill_manual(values = c("percent_r2_LLEA"="#2C8C99",
                               "language_difference"="#42D9C8"),
                    labels = c('LLEA + language', 'LLEA alone')) +
  geom_bar(data = fig3_PVE_maize10PCs %>% filter(dataset == 'Haynie') %>%
             filter(str_detect(model, 'PC')) %>%
             mutate(PC = factor(model, levels = c('PC10', 'PC9', 'PC8', 'PC7', 'PC6', 'PC5', 'PC4', 'PC3', 'PC2', 'PC1'))),
           stat="identity",position="stack",
           fill="#999999",
           aes(y=percent, x= 'top 10 maize PCs', color = PC))+
  scale_color_viridis_d(option = "G") +
  facet_grid(dataset ~ language_set) +
  theme_bw() +
  guides(fill = guide_legend(title = "Model: ", position = "bottom"),
         color="none") +
  xlab("")+ylab("Genetic variation explained (%)") +
  theme(legend.position = 'bottom')

regression_plots_maize_PCs = ggplot(fig3_PVE_maize10PCs, aes(x = c('models', 'top 10 maize PCs'))) +
  geom_bar(data = fig3_PVE_maize10PCs %>% filter(model %in% c('percent_r2_LLEA', 'language_difference')),
           stat = 'identity', position = 'stack',
           aes(y = percent, x = 'models', fill = model)) +
  scale_fill_manual(values = c("percent_r2_LLEA"="#2C8C99",
                               "language_difference"="#42D9C8"),
                    labels = c('with language', 'LLEA alone')) +
  geom_bar(data = fig3_PVE_maize10PCs %>% 
             filter(str_detect(model, 'PC')) %>% 
             mutate(PC = factor(model, levels = c('PC10', 'PC9', 'PC8', 'PC7', 'PC6', 'PC5', 'PC4', 'PC3', 'PC2', 'PC1'))),
           stat="identity",position="stack",
           fill="#999999",
           aes(y=percent, x= 'top 10 maize PCs', color = PC))+
  scale_color_viridis_d(option = "G") +
  facet_grid(dataset ~ language_set) +
  theme_bw() +
  guides(fill = guide_legend(title = "Model: ", position = "bottom", reverse = T),
         color="none") +
  xlab("")+ylab("Genetic variation explained (%)") +
  theme(legend.position = 'top',
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 9),
        legend.key.height = unit(.25, "cm"),
        legend.spacing = unit(c(.1, .1, .1, .1), 'cm'),
        legend.box.margin = unit(c(.1, .1, .1, .1), 'cm'),
        plot.margin = unit(c(.1,.1,.1,.1), 'cm'),
        axis.text.x = element_text(angle = 30, hjust = 1))


fig3_regression_plots_maize_PCs = ggplot(fig3_PVE_maize10PCs %>% filter(dataset == 'Haynie'), aes(x = c('models', 'top 10 maize PCs'))) +
  geom_bar(data = fig3_PVE_maize10PCs %>% filter(dataset == 'Haynie') %>% filter(model %in% c('percent_r2_LLEA', 'language_difference')),
           stat = 'identity', position = 'stack',
           aes(y = percent, x = 'models', fill = model)) +
  scale_fill_manual(values = c("percent_r2_LLEA"="#2C8C99",
                               "language_difference"="#42D9C8"),
                    labels = c('with language', 'LLEA alone')) +
  geom_bar(data = fig3_PVE_maize10PCs %>% filter(dataset == 'Haynie') %>% 
             filter(str_detect(model, 'PC')) %>% 
             mutate(PC = factor(model, levels = c('PC10', 'PC9', 'PC8', 'PC7', 'PC6', 'PC5', 'PC4', 'PC3', 'PC2', 'PC1'))),
           stat="identity",position="stack",
           fill="#999999",
           aes(y=percent, x= 'top 10 maize PCs', color = PC))+
  scale_color_viridis_d(option = "G") +
  facet_grid(~ language_set) +
  theme_bw() +
  guides(fill = guide_legend(title = "Model: ", position = "bottom", reverse = T),
         color="none", ) +
  xlab("")+ylab("Genetic variation explained (%)") +
  theme(legend.position = 'top',
        legend.spacing = unit(c(.05, .05, .05, .05), 'cm'),
        legend.box.margin = unit(c(.05, .05, .05, .05), 'cm'),
        plot.margin = unit(c(.05,.05,.05,.05), 'cm'),
        axis.text.x = element_text(angle = 30, hjust = 1))

all_haynie_polygon_map +
  theme(legend.position = 'top',
        legend.box.margin = unit(c(0, 0, 0, 0), 'cm'),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 5),
        legend.key.height = unit(.05, "cm"),
        plot.margin = unit(c(t = -0.05, r= -0.05, b = -0.10, l= -0.05), 'cm'),
        axis.ticks = element_blank(),
        axis.text = element_blank())

## no mayan only languages map
# (figure_3 = plot_grid(plot_grid(all_haynie_polygon_map +
#                                   theme(legend.position = 'top',
#                                         legend.box.margin = unit(c(0, 0, 0, 0), 'cm'),
#                                         legend.key.height = unit(.05, "cm"),
#                                         legend.title = element_text(size = 6),
#                                         legend.text = element_text(size = 6),
#                                         plot.margin = unit(c(t = -0.10, 
#                                                              r= -0.05, 
#                                                              b = -0.10, 
#                                                              l= -0.05), 'cm'),
#                                         axis.ticks = element_blank(),
#                                         axis.text = element_blank()),
#                                 haynie_fst_plot + ggtitle('') + ylab('Fst'),
#                                 ncol = 2, labels = c('A', 'B')), 
#                      plot_grid(mayan_tree_figure + xlim(0, 0.6) , 
#                                fig3_regression_plots_maize_PCs + 
#                                  theme(legend.title = element_text(size = 7),
#                                        legend.text = element_text(size = 7),
#                                        legend.key.height = unit(.15, "cm")), 
#                                ncol = 2, labels = c('C', 'D'), rel_widths = c(1, 1.1)),
#                      nrow = 2, 
#                      rel_heights = c(1, 1.2),
#                      labels = c('', ''))
# )

mayan_haynie_map = ggmap(mayan_terrain_kernels_transformed) +
  coord_sf(crs = st_crs(3857)) +
  geom_sf(aes(fill = Name),
          data = st_transform(haynie_lang %>% 
                                filter(Name %in% mayan_haynie_mapping$`Polygon Name`), 
                              3857),
          inherit.aes = F,
          alpha = 0.4) +
  geom_sf(data = st_transform(mayan_haynie_accessions, 3857), 
          inherit.aes = F,
          size = .4, 
          alpha = 0.7) +
  theme(legend.position = "none",
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        title = element_text(size = 8)) +
  ggtitle('Haynie Mayan languages')

## with mayan languages map
(figure_3 = plot_grid(plot_grid(all_haynie_polygon_map +
                                  theme(legend.position = 'top',
                                        legend.box.margin = unit(c(0, 0, 0, 0), 'cm'),
                                        legend.key.height = unit(.05, "cm"),
                                        legend.title = element_text(size = 6),
                                        legend.text = element_text(size = 6),
                                        plot.margin = unit(c(t = -0.10, 
                                                             r= 0, 
                                                             b = -0.10, 
                                                             l= 0), 'cm'),
                                        axis.ticks = element_blank(),
                                        axis.text = element_blank(),
                                        axis.title = element_blank()),
                                mayan_haynie_map,
                                mayan_tree_figure + xlim(0, 0.7),
                                ncol = 3, labels = c('A', 'B', 'C')), 
                      plot_grid(haynie_fst_plot + ggtitle('') + ylab('Fst'), 
                                fig3_regression_plots_maize_PCs + 
                                  theme(legend.title = element_text(size = 8),
                                        legend.text = element_text(size = 8),
                                        legend.key.height = unit(.20, "cm"),
                                        axis.title.y = element_text(size = 8),
                                        axis.text.x = element_text(size = 7)), 
                                ncol = 2, labels = c('D', 'E'), rel_widths = c(1, 1)),
                      nrow = 2, 
                      rel_heights = c(1.1, 1),
                      labels = c('', ''))
)
  
ggsave('fig3_prototype_111825.pdf',
       plot = figure_3, 
       device = 'pdf', 
       path = '/Users/liforrest/Documents/Projects/maize_co-evolution/', 
       width = 18, height = 15, units = 'cm', bg = '#ffffff')

ggsave('regression_plots_maize_PCs.pdf',
       plot = regression_plots_maize_PCs, 
       device = 'pdf', 
       path = '/Users/liforrest/Documents/Projects/maize-linguistics/plots/', 
       width = 10, height = 10, units = 'cm', bg = '#ffffff')
  
NL_Haynie_fig = plot_grid(Fst_plot, regression_plots_maize_PCs, labels = c('A', 'B'), 
                          ncol = 2)

ggsave('NL_Haynie.pdf',
       plot = NL_Haynie_fig,
       device = 'pdf',
       path = '/Users/liforrest/Documents/Projects/maize_co-evolution/', 
       width = 18, height = 12, units = 'cm', bg = '#ffffff')

# use regression_variance_plot from scratch.R
  
# pcoa_regression =  ggplot(regression_pcoa_means_pivot %>% filter(model %in% c('adj_language', 'adj_geography', 'adj_full', 'adj_elevation', 'adj_admix')), 
#          aes(x = model, y = mean_r2 * 100)) +
#     facet_grid(rows = vars(dataset), cols = vars(language_set)) +
#     geom_col() +
#     theme(axis.text.x = element_text(angle = 90, vjust = 0.2)) +
#     coord_flip() +
#     ggtitle('PCoA')


  
  
  
fig2_plots = lapply(list('../plots/Yucatec_manhattan.png', '../plots/Otomi_manhattan.png', '../plots/feems_haynie_10res_HaynieNL.png'), 
                    function(x) {
                      img = as.raster(readPNG(x))
                      rasterGrob(img, interpolate = F)
                    }
)
  
(figure_2 = plot_grid(plot_grid(fig2_plots[[1]], fig2_plots[[2]], yucatec_plot, otomi_04_plot, otomi_09_plot, cora_04_plot, nrow = 5, rel_heights = c(1.5, 1.5, 1, 1, 1)), 
                         plot_grid(as_grob(pcoa_regression), fig2_plots[[3]], ncol = 1, nrow = 2), ncol = 2)
)
plot(figure_2)
  


marrangeGrob(fig2_plots, nrow = 3, ncol =1)

ggplot(as.ggplot(rasterGrob(as.raster(readPNG('../plots/feems_haynie_10res_HaynieNL.png')))))

### create plot for human association loci ########################################################
# human_gwas_results<-read_tsv("../results/21km_raw_LDprune_MAF0.01_LatLongElevAdmixPC1_10_lm_results.tsv")
# human_gwas_pos = read_tsv("../results/LDprune_MAF0.01_21km_MaizeGBS_SNP_pos.tsv")
# human_gwas_results = cbind(human_gwas_pos, human_gwas_results)
# colnames(human_gwas_results) = c('CHROM', 'POS', 'SNP', 'ID', 'r.squared', 'adj.r.squared')
# 
# human_max_bp = human_gwas_results %>% group_by(CHROM) %>% 
#   summarise(max_bp = max(POS), center = floor(mean(POS))) %>% 
#   mutate(bp_add = lag(cumsum(max_bp), default = 0),
#          center_add = bp_add + center)
# human_gwas_results = human_gwas_results %>% inner_join(human_max_bp, by = 'CHROM') %>% mutate(cumulative_bp = POS + bp_add)
# 
# human_gwas_plot = ggplot(human_gwas_results, aes(x = cumulative_bp, y = adj.r.squared, color = as.factor(CHROM))) +
#   geom_point(alpha = 0.3) +
#   scale_x_continuous(
#     label = max_bp %>% pull(CHROM),
#     breaks = max_bp %>% pull(center_add)) +
#   scale_color_manual(values = rep(
#     c("grey4", "grey30"),
#     unique(length(max_bp %>% pull(CHROM)))
#   )) + 
#   labs(x = 'Chromosome') +
#   theme_bw() +
#   theme(legend.position = 'none',
#         axis.text = element_text(size = 14)) +
#   ggtitle('R-squared outliers from human PCs')
# 
# 
# (figure_4 = plot_grid(all_gwas_plots, 
#                         plot_grid(yucatec_01_plot, otomi_04_plot, otomi_09_plot, cora_04_plot, nrow = 4, labels = c('D', 'E', 'F', 'G')),
#                         ncol = 2)
# )
# 
# ggsave('fig4_prototype_101725.pdf',
#        plot = figure_4, 
#        device = 'pdf', 
#        path = '/Users/liforrest/Documents/Projects/maize_co-evolution/', 
#        width = 18, height = 16, units = 'cm')

#### full americas for both human and maize ###########

americas_terrain_kernels = get_stadiamap(c(left = -120, bottom = -40, right = -35, top = 33),
                                            zoom = 5, maptype = 'stamen_toner_lite', color = 'bw', messaging = F)
americas_terrain_kernels_transformed = ggmap_bbox(americas_terrain_kernels)

maize_seed_data <- read.delim("../data/SeeDs_passport_fulll_2017.04.25_16.08.12.txt")
bankUpdated = read.csv('../data/selected_genotypeIDs.csv')
colnames(bankUpdated) = c('Unique.ID', 'Sample')

maize_coord_coevolution = data.frame(Longitude = maize_seed_data$locations_longitude, 
                         Latitude = maize_seed_data$locations_latitude, 
                         Elevation= maize_seed_data$locations_elevation, 
                         Sample.ID.of.DNA.from.single.plants.used.in.GWAS= maize_seed_data$Sample.ID.of.DNA.from.single.plants.used.in.GWAS)

# make a copy of the dataframe for left joining/merging later
maize_coord_coevolution$Unique.ID = bankUpdated[match(maize_coord_coevolution$Sample.ID.of.DNA.from.single.plants.used.in.GWAS,
                                                                  bankUpdated$Sample), 'Unique.ID']

human_maize_accessions = read.table('../data/human-maize-accessions.txt', header = T)
human_maize_accessions_gwas = merge(human_maize_accessions, bankUpdated, 
                                    by.x = 'Sample_ID_of_DNA_from_single_plants_used_in_GWAS.maize',
                                    by.y = 'Sample') %>% unique()

all_language_ID = c(all_haynie_accessions$Unique.ID, all_NL_accessions$Unique.ID)

maize_coord_coevolution_final = maize_coord_coevolution %>% filter(Unique.ID %in% all_language_ID | 
                                     Unique.ID %in% human_maize_accessions_gwas$Unique.ID) %>% 
  mutate('Only Language' = (Unique.ID %in% all_language_ID) & 
           !(Unique.ID %in% human_maize_accessions_gwas$Unique.ID),
         'Only Human' = !(Unique.ID %in% all_language_ID) & 
           (Unique.ID %in% human_maize_accessions_gwas$Unique.ID),
         'Both' = (Unique.ID %in% all_language_ID) &
           (Unique.ID %in% human_maize_accessions_gwas$Unique.ID) )

# maize_coord_coevolution = projectMaizeCoordinates(maize_coord_coevolution)
coordinates(maize_coord_coevolution_final) = ~ Longitude + Latitude
proj4string(maize_coord_coevolution_final) = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
maize_coord_coevolution_final <- st_as_sf(maize_coord_coevolution_final,coords = 1:2)
maize_coord_coevolution_final <- maize_coord_coevolution_final %>% 
  st_transform(crs = 3857)

(all_accessions_map = ggmap(americas_terrain_kernels_transformed) +
    geom_sf(data = maize_coord_coevolution_final %>% filter(`Only Human`), 
            # aes(fill = `Only Language`,
            #     shape = `Only Human`),
            # aes(x = Longitude, y = Latitude),
            inherit.aes = F,
            size = 0.5,
            alpha = 0.7) +
    theme(legend.position = "bottom") +
    ggtitle('')
)

### fig4 construction ###################
yucatan_genotypes_transform = yucatan_genotypes_transform %>% 
  mutate(S1_164490966_dosage = case_when(
    S1_164490966 == "0|0" ~ "Ref",
    S1_164490966 == "0|1" ~ "Het",
    S1_164490966 == "1|0" ~ "Het",
    S1_164490966 == "1|1" ~ "Alt"),
    S2_197685637_dosage = case_when(
      S2_197685637 == "0|0" ~ "Ref",
      S2_197685637 == "0|1" ~ "Het",
      S2_197685637 == "1|0" ~ "Het",
      S2_197685637 == "1|1" ~ "Alt"))

## yucatan map 
(yucatan_alleles_fig4 = ggmap(mayan_terrain_kernels) +
    geom_sf(data = yucatan_genotypes_transform, aes(color = S1_164490966_dosage),
            size = 0.5,
            alpha = 0.7,
            inherit.aes = F) +
    guides(color = guide_legend(title = "S1_164490966 genotype") ) +
    scale_color_manual(values = c('red', 'purple', 'blue'))+
    theme(legend.position = "right",
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 6),
          axis.text.y = element_text(size = 6),
          legend.text=element_text(size=8),
          legend.title = element_text(size = 8),
          plot.title = element_text(size = 12)) +
    ggtitle('Yucatan GWAS peak allele') +
    xlab('') +
    ylab('')
)

yucatec_gwas_sig_SNPs = yucatec_gwas_results %>% filter(gwas_sig) %>% pull(SNP)

(figure_4 = plot_grid(plot_grid(manual_manhattan(yucatec_gwas_results, yucatec_gwas_sig_SNPs, 
                                                 sig = 1e-08, chr_col = 'CHR', bp_col = 'BP', p_col = 'P', title = '') + 
                                  theme(axis.text.x = element_blank()),
                                yucatec_chr_plot + 
                                  labs(title = NULL) +
                                  theme(panel.grid.major.x = element_blank(),
                                        panel.grid.minor.x = element_blank())+
                                  scale_y_reverse(),
                                  ncol = 1, nrow = 2, labels = c('A', 'B')), 
                      yucatan_alleles_fig4 + 
                        theme(legend.position = 'top') + 
                        ggtitle(''),
                      ncol = 2, labels = c('', 'C')))

ggsave('fig4_prototype_111025.pdf',
       plot = figure_4, 
       device = 'pdf', 
       path = '/Users/liforrest/Documents/Projects/maize_co-evolution/', 
       width = 18, height = 10, units = 'cm', bg = '#ffffff')

ggsave('fig4_prototype_111025.png',
       plot = figure_4, 
       device = 'png', 
       path = '/Users/liforrest/Documents/Projects/maize_co-evolution/', 
       width = 18, height = 10, units = 'cm', bg = '#ffffff')
