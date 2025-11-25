source('apikey.R')
source('languageFunctions.R')
source('gwasAnalysis.R')

### create plot for gwas ########################################################

## make changeable ggplots for linguistic gwas
yucatec_gwas_plot = manual_manhattan(yucatec_gwas_results, yucatec_highlight, 
                                     sig = 1e-08, chr_col = 'CHR', bp_col = 'BP', p_col = 'P', title = 'Yucatan Mayan GWAS')
otomi_gwas_plot = manual_manhattan(otomi_gwas_results, otomi_highlight, 
                                   sig = 1e-08, chr_col = 'CHR', bp_col = 'BP', p_col = 'P', title = 'Otomi GWAS')
cora_gwas_plot = manual_manhattan(cora_gwas_results, cora_highlight, 
                                  sig = 1e-08, chr_col = 'CHR', bp_col = 'BP', p_col = 'P', title = 'Cora/Huichol GWAS')
chol_gwas_plot = manual_manhattan(chol_gwas_results, chol_highlight, 
                                  sig = 1e-08, chr_col = 'CHR', bp_col = 'BP', p_col = 'P', title = 'Chol/Chorti/Cholti GWAS')
pipil_gwas_plot = manual_manhattan(pipil_gwas_results, pipil_highlight, 
                                   sig = 1e-08, chr_col = 'CHR', bp_col = 'BP', p_col = 'P', title = 'Pipil and Eastern/Northern Puebla Nahuatl GWAS')
matlatzinca_gwas_plot = manual_manhattan(matlatzinca_gwas_results, matlatzinca_highlight, 
                                  sig = 1e-08, chr_col = 'CHR', bp_col = 'BP', p_col = 'P', title = 'Matlatzinca/Mazahua GWAS')
mangue_gwas_plot = manual_manhattan(mangue_gwas_results, mangue_highlight, 
                                  sig = 1e-08, chr_col = 'CHR', bp_col = 'BP', p_col = 'P', title = 'Mangue/Chiapenec/Mephaa GWAS')


(all_gwas_plots_2 = plot_grid(yucatec_gwas_plot, otomi_gwas_plot, cora_gwas_plot, pipil_gwas_plot, chol_gwas_plot, matlatzinca_gwas_plot, mangue_gwas_plot,
                            labels = 'AUTO', ncol = 2, nrow = 4))

ggsave('all_gwas_plots_2.png', 
       plot = all_gwas_plots_2,
       device = 'png', 
       path = '../plots/', 
       width = 10, height = 10, units = 'in')

ggsave('all_gwas_plots_2.pdf', 
       plot = all_gwas_plots_2,
       device = 'pdf', 
       path = '../plots/', 
       width = 10, height = 10, units = 'in')


# create plots for SweeD ----------------------------------------------------------------------
#### sweed plots for yucatec full genome for multiple chromosomes and for plotting in manuscript

# Yucatec individual chromosomes
(yucatec_01_plot = plotSweeD(sweed_yucatec %>% filter(Chr == 1), 
                             chr = 1, pos_intervals = list(c(165928841 - 500000, 167679000 + 500000)), language = 'Yucatan Mayan'))
(yucatec_02_plot = plotSweeD(sweed_yucatec %>% filter(Chr == 2), 
                             chr = 2, pos_intervals = list(c(204174296 - 500000, 204174757 + 500000)), language = 'Yucatan Mayan'))
(yucatec_08_plot = plotSweeD(sweed_yucatec %>% filter(Chr == 8), 
                             chr = 8, pos_intervals = list(c(27475338 - 500000, 27475338 + 500000)), language = 'Yucatan Mayan'))

(cora_04_plot = plotSweeD(cora_04, chr = 4, pos_intervals = list(c(10384903 - 500000, 10384903 + 500000)), language = 'Cora/Huichol'))
(otomi_04_plot = plotSweeD(otomi_04, chr = 4, pos_intervals = list(c(178304566 - 500000, 178304566 + 500000)), language = 'Otomi'))
(otomi_05_plot = plotSweeD(otomi_05, chr = 5, pos_intervals = list(c(199268719 - 500000, 199268719 + 500000)), language = 'Otomi'))
(otomi_09_plot = plotSweeD(otomi_09, chr = 9, pos_intervals = list(c(38991589 - 500000, 38991589 + 500000)), language = 'Otomi'))
(matlatzinca_08_plot = plotSweeD(matlatzinca_08, chr = 8, pos_intervals = list(c(5366101 - 500000, 5366101 + 500000), 
                                                                                 c(29095461 - 500000, 29095461 + 500000)), language = 'Matlatzinca/Mazahua'))
# (matlatzinca_08_plot_b = plotSweeD(matlatzinca_08, chr = 8, pos_intervals = list(c(29095461 - 500000, 29095461 + 500000)), language = 'Matlatzinca'))
(mangue_01_plot = plotSweeD(mangue_01, chr = 1, pos_intervals = list(c(28224363 - 500000, 28224363 + 500000)), language = 'Mangue/Chiapenec/Mephaa'))
(mangue_03_plot = plotSweeD(mangue_03, chr = 3, pos_intervals = list(c(181819647 - 500000, 181819647 + 500000)), language = 'Mangue/Chiapenec/Mephaa'))
(pipil_09_plot = plotSweeD(pipil_09, chr = 9, pos_intervals = list(c(148360542 - 500000, 148360542 + 500000)), language = 'Pipil and NE Puebla Nahuatl'))
(chol_05_plot = plotSweeD(chol_05, chr = 5, pos_intervals = list(c(213908202 - 500000, 213908202 + 500000)), language = 'Chol/Chorti/Cholti'))
(chol_07_plot = plotSweeD(chol_07, chr = 7, pos_intervals = list(c(158526307 - 500000, 158526307 + 500000)), language = 'Chol/Chorti/Cholti'))

## do not include mangue_01, maybe not yucatec_08, otomi_09
(all_sweed_plots_2 = plot_grid(yucatec_01_plot, yucatec_02_plot, otomi_04_plot, otomi_05_plot, 
                               mangue_03_plot, matlatzinca_08_plot, cora_04_plot, pipil_09_plot, chol_05_plot, chol_07_plot,
                             labels = 'AUTO', ncol= 3))
ggsave('all_sweed_plots.png', 
       plot = all_sweed_plots_2,
       device = 'png', 
       path = '../plots/', 
       width = 18, height = 13, units = 'cm')

ggsave('all_sweed_plots.pdf', 
       plot = all_sweed_plots_2,
       device = 'pdf', 
       path = '../plots/', 
       width = 18, height = 13, units = 'cm')



# create plot for Fst -------------------------------------------------------------------------

NL_fst_plot = ggplot(NL_fst_results, aes(x = cumulative_bp, y = WEIR_AND_COCKERHAM_FST, color = as.factor(CHROM))) +
  geom_scattermore(alpha = 0.3, pixels = c(768, 512), pointsize = 3.2) +
  scale_x_continuous(
    label = max_bp %>% pull(CHROM),
    breaks = max_bp %>% pull(center_add)) +
  scale_color_manual(values = rep(
    c("grey4", "grey30"),
    unique(length(max_bp %>% pull(CHROM)))
  )) + 
  labs(x = 'Chromosome', y = 'Fst') +
  theme_bw() +
  theme(legend.position = 'none',
        axis.text = element_text(size = 8)) +
  ggtitle('NL Fst between language families')

haynie_fst_plot = ggplot(haynie_fst_results, aes(x = cumulative_bp, y = WEIR_AND_COCKERHAM_FST, color = as.factor(CHROM))) +
  geom_scattermore(alpha = 0.30, pixels = c(768, 512), pointsize = 3.2) +
  # geom_point(alpha = 0.7)+
  scale_x_continuous(
    label = max_bp %>% pull(CHROM),
    breaks = max_bp %>% pull(center_add)) +
  scale_color_manual(values = rep(
    c("grey4", "grey30"),
    unique(length(max_bp %>% pull(CHROM)))
  )) + 
  labs(x = 'Chromosome', y = 'Fst') +
  theme_bw() +
  theme(legend.position = 'none',
        axis.text = element_text(size = 8)) +
  ggtitle('Haynie Fst between language families')


(Fst_plot = plot_grid(NL_fst_plot, haynie_fst_plot, nrow = 2))

ggsave('Fst_plots.png', 
       plot = Fst_plot,
       device = 'png', 
       path = '../plots/', 
       width = 8, height = 6, units = 'in')

ggsave('Fst_plots.pdf', 
       plot = Fst_plot,
       device = 'pdf', 
       path = '../plots/', 
       width = 8, height = 6, units = 'in')


ggplot(haynie_fst_results %>% filter(CHROM == 9), aes(x = cumulative_bp, 
                                                      y = WEIR_AND_COCKERHAM_FST, 
                                                      color = POS > 1.480e8 & POS < 1.490e8)) +
  geom_point(alpha = 0.3) +
  scale_x_continuous(
    label = max_bp %>% pull(CHROM),
    breaks = max_bp %>% pull(center_add)) +
  scale_color_manual(values = rep(
    c("grey4", "red"),
    unique(length(max_bp %>% pull(CHROM)))
  )) + 
  labs(x = 'Chromosome') +
  theme_bw() +
  theme(legend.position = 'none',
        axis.text = element_text(size = 14)) +
  ggtitle('Haynie Fst between language families')
          


