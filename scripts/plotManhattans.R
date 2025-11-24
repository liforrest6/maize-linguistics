source('apikey.R')
source('languageFunctions.R')


manual_manhattan = function(results, highlight_SNP_list, chr_filter = NA, sig = 1e-5,
                            chr_col = 'CHR', bp_col = 'BP', p_col = 'P', 
                            title = 'Multivariate GEA'){
  if(!is.na(chr_filter)) {
    results = results %>% filter((!!sym(chr_col)) == chr_filter)
  }
  
  max_BP <- results |>
    group_by((!!sym(chr_col))) |>
    summarise(max_bp = max((!!sym(bp_col)))) |>
    mutate(bp_add = lag(cumsum(max_bp), default = 0)) |>
    dplyr::select((!!sym(chr_col)), bp_add)
  
  results_add <- results |>
    inner_join(max_BP, by = chr_col) |>
    mutate(bp_cum = (!!sym(bp_col)) + bp_add)
  
  axis_set <- results_add |>
    group_by((!!sym(chr_col))) |>
    summarize(center = mean(bp_cum))
  
  ylim <- results_add |>
    filter((!!sym(p_col)) == min((!!sym(p_col)))) |>
    mutate(ylim = abs(floor(log10((!!sym(p_col))))) + 2) |>
    pull(ylim)
  
  # sig <- 1e-5
  
  (manhplot <- ggplot(results_add, aes(
    x = bp_cum, y = -log10((!!sym(p_col))),
    color = as.factor((!!sym(chr_col))), size = -log10((!!sym(p_col)))
  )) +
      geom_hline(
        yintercept = -log10(sig), color = "blue",
        linetype = "dashed"
      ) +
      geom_scattermore(alpha = 1, pointsize = 4.2, pixels = c(1280, 1024)) +
      # geom_point(alpha = 0.75) +
      scale_x_continuous(
        label = axis_set %>% pull((!!sym(chr_col))),
        breaks = axis_set$center
      ) +
      scale_y_continuous(expand = c(0, 0), limits = c(0, ylim)) +
      scale_color_manual(values = rep(
        c("grey4", "grey30"),
        unique(length(axis_set %>% pull(!!sym(chr_col))))
      )) +
      # geom_scattermore(data = results_add[results_add$SNP %in% highlight_SNP_list,],
      #                  alpha = 1, pointsize = 4.2, pixels = c(1280, 1024), color = 'red') + 
      geom_point(data = results_add[results_add$SNP %in% highlight_SNP_list,],
                       alpha = 1, size = 0.5, color = 'red') + 
      # scale_size_continuous(range = c(0.5, 3)) +
      labs(
        x = NULL,
        y = "-log<sub>10</sub>(p)"
      ) +
      theme_bw() +
      theme(
        legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title.y = element_markdown(),
        title = element_text(size = 9),
        axis.text.x = element_text(angle = 0, size = 8, vjust = 0.5)
      ) +
      ggtitle(title)
  )
  
}

plotSweeD = function(results, chr, language, pos_intervals = NA,
                            title = 'Multivariate GEA'){
  highlight_snps = c()
  if(all(!is.na(pos_intervals))) {
    for(interval in pos_intervals) {
      # highlight_snps = c(highlight_snps, results %>% filter(Position > pos_intervals[[1]][1] &
      #                                                         Position < pos_intervals[[1]][2]) %>% pull(SNP))
      highlight_snps = c(highlight_snps, results %>% filter(Position > interval[1] &
                                                              Position < interval[2]) %>% pull(SNP))
    }
    # highlight_snps = results %>% filter(Position > pos_intervals[[1]][1] &
    #                                       Position < pos_intervals[[1]][2]) %>% 
    #   pull(SNP)
    # results$highlighted = results[results$SNP %in% highlight_snps]
  }
  ggplot(results, aes(x = Position, y = Likelihood, 
                       color = SNP %in% highlight_snps)) +
    # geom_point(alpha = 1) + 
    geom_scattermore(alpha = 1, pointsize = 4.2, pixels = c(1280, 1024)) + 
    scale_color_manual(values = c('black', 'red')) +
    scale_x_continuous(n.breaks = 4, labels = unit_format(unit = "M", scale = 1e-6)) +
    labs(color='presence in GWAS') +
    ggtitle(sprintf('%s, chr %d', language, chr)) + 
    theme_bw() +
    theme(legend.position="none",
          axis.text = element_text(size = 5),
          title = element_text(size = 6),
          axis.title = element_text(size = 6))
}

### read results from gwas ########################################################
yucatec_gwas_results = vroom('../results/mayan_gemma_output/Yucatec_GWAS_results.txt')
otomi_gwas_results = vroom('../results/otomanguean_gemma_output/Otomi_GWAS_results.txt')
cora_gwas_results = vroom('../results/aztecan_gemma_output/Cora_GWAS_results.txt')
chol_gwas_results = vroom('../results/mayan_gemma_output/Chol_GWAS_results.txt')
matlatzinca_gwas_results = vroom('../results/otomanguean_gemma_output/Matlatzinca_GWAS_results.txt')
mangue_gwas_results = vroom('../results/otomanguean_gemma_output/Mangue_GWAS_results.txt')
pipil_gwas_results = vroom('../results/aztecan_gemma_output/Pipil_GWAS_results.txt')

yucatec_highlight = yucatec_gwas_results %>% filter(P < 1e-8)
otomi_highlight = otomi_gwas_results %>% filter(P < 1e-8)
cora_highlight = cora_gwas_results %>% filter(P < 1e-8)
chol_highlight = chol_gwas_results %>% filter(P < 1e-8)
matlatzinca_highlight = matlatzinca_gwas_results %>% filter(P < 1e-8)
mangue_highlight = mangue_gwas_results %>% filter(P < 1e-8)
pipil_highlight = pipil_gwas_results %>% filter(P < 1e-8)

write.table(yucatec_highlight$SNP, '../results/mayan_gemma_output/yucatec_highlight.txt', sep = '\t', quote = F, row.names = F, col.names =F)
write.table(otomi_highlight$SNP, '../results/otomanguean_gemma_output/otomi_highlight.txt', sep = '\t', quote = F, row.names = F, col.names =F)
write.table(cora_highlight$SNP, '../results/aztecan_gemma_output/cora_highlight.txt', sep = '\t', quote = F, row.names = F, col.names =F)
write.table(chol_highlight$SNP, '../results/mayan_gemma_output/chol_highlight', sep = '\t', quote = F, row.names = F, col.names =F)
write.table(matlatzinca_highlight$SNP, '../results/otomanguean_gemma_output/matlatzinca_highlight', sep = '\t', quote = F, row.names = F, col.names =F)
write.table(mangue_highlight$SNP, '../results/otomanguean_gemma_output/mangue_highlight', sep = '\t', quote = F, row.names = F, col.names =F)
write.table(pipil_highlight$SNP, '../results/aztecan_gemma_output/pipil_highlight', sep = '\t', quote = F, row.names = F, col.names =F)

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

# (all_gwas_plots = plot_grid(yucatec_gwas_plot, otomi_gwas_plot, cora_gwas_plot, 
#                             labels = 'AUTO', nrow = 3))
# 
# ggsave('all_gwas_plots.png', 
#        plot = all_gwas_plots,
#        device = 'png', 
#        path = '../plots/', 
#        width = 8, height = 6, units = 'in')
# 
# ggsave('all_gwas_plots.pdf', 
#        plot = all_gwas_plots,
#        device = 'pdf', 
#        path = '../plots/', 
#        width = 8, height = 6, units = 'in')

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

### create plot for sweed ########################################################
#### sweed plots for yucatec, full genome for multiple chromosomes and for plotting in manuscript
readSweeD = function(chr) {
  sweed_file = data.frame(vroom(sprintf('../results/SweeD_selection/Sweed_Report.Yucatec-%02d', chr), skip = 2))
  sweed_file$Chr = chr
  sweed_file$SNP = paste(sweed_file$Chr, sweed_file$Position,sep = '_')
  sweed_file
}
sweed_files = lapply(c(1:10), readSweeD)
sweed_yucatec = bind_rows(sweed_files)
max_bp = sweed_yucatec %>% group_by(Chr) %>% 
  summarise(max_bp = max(Position), center = floor(mean(Position))) %>% 
  mutate(bp_add = lag(cumsum(max_bp), default = 0),
         center_add = bp_add + center)
sweed_yucatec_results = sweed_yucatec %>% inner_join(max_bp, by = 'Chr') %>% mutate(cumulative_bp = Position + bp_add)

yucatec_chr1_interval = c(165928841 - 500000, 167679000 + 500000)
yucatec_chr2_interval = c(204174296 - 500000, 204174757 + 500000)
yucatec_chr8_interval = c(27475338 - 500000, 27475338 + 500000)

yucatec_gwas_results$gwas_sig = (yucatec_gwas_results$CHR == 1 & yucatec_gwas_results$BP > yucatec_chr1_interval[1] & yucatec_gwas_results$BP < yucatec_chr1_interval[2]) |
  (yucatec_gwas_results$CHR == 2 & yucatec_gwas_results$BP > yucatec_chr2_interval[1] & yucatec_gwas_results$BP < yucatec_chr2_interval[2]) |
  (yucatec_gwas_results$CHR == 8 & yucatec_gwas_results$BP > yucatec_chr8_interval[1] & yucatec_gwas_results$BP < yucatec_chr8_interval[2])

(yucatec_chr_plot = ggplot(sweed_yucatec_results, aes(x = cumulative_bp, y = Likelihood, color = as.factor(Chr))) +
  # geom_scattermore(alpha = 1, pixels = c(768, 512), pointsize = 3.2) +
  geom_scattermore(alpha = 1, pointsize = 4.2, pixels = c(1280, 1024)) +
  scale_x_continuous(
    label = max_bp %>% pull(Chr),
    breaks = max_bp %>% pull(center_add)) +
  scale_color_manual(values = rep(
    c("grey4", "grey30"),
    unique(length(max_bp %>% pull(Chr)))
  )) +
    geom_scattermore(
      data = subset(sweed_yucatec_results, gwas_sig == TRUE),
      color = "red", alpha = 1, pointsize = 4.2, pixels = c(1280, 1024)
    ) +
  labs(x = 'Chromosome') +
  theme_bw() +
  theme(legend.position = 'none',
        axis.text = element_text(size = 8)) +
  ggtitle('Yucatan Mayan SweeD results')
)

# Otomi ---------------------------------------------------------------------------------------

otomi_04 = vroom('../results/SweeD_selection/SweeD_Report.Otomi-04', skip = 2)
otomi_04$Chr = 4
otomi_04$SNP = paste(otomi_04$Chr, otomi_04$Position, sep = '_')

otomi_05 = vroom('../results/SweeD_selection/SweeD_Report.Otomi-05', skip = 2)
otomi_05$Chr = 5
otomi_05$SNP = paste(otomi_05$Chr, otomi_05$Position, sep = '_')

otomi_09 = vroom('../results/SweeD_selection/SweeD_Report.Otomi-09', skip = 2)
otomi_09$Chr = 9
otomi_09$SNP = paste(otomi_09$Chr, otomi_09$Position, sep = '_')

# Cora/Huichol --------------------------------------------------------------------------------
cora_04 = vroom('../results/SweeD_selection/SweeD_Report.Cora-Huichol-04', skip = 2)
cora_04$Chr = 4
cora_04$SNP = paste(cora_04$Chr, cora_04$Position,sep = '_')

# Chol/Chontal/Chorti/Cholti ----------------------------------------------------------------------------------------
chol_05 = vroom('../results/SweeD_selection/SweeD_Report.Chol-05', skip = 2)
chol_05$Chr = 5
chol_05$SNP = paste(chol_05$Chr, chol_05$Position,sep = '_')

chol_07 = vroom('../results/SweeD_selection/SweeD_Report.Chol-07', skip = 2)
chol_07$Chr = 7
chol_07$SNP = paste(chol_07$Chr, chol_07$Position,sep = '_')

# Pipil/Eastern + Northern Puebla Nahuatl ---------------------------------------------------------------------------------------
pipil_09 = vroom('../results/SweeD_selection/SweeD_Report.Pipil-09', skip = 2)
pipil_09$Chr = 9
pipil_09$SNP = paste(pipil_09$Chr, pipil_09$Position,sep = '_')

# Mangue/Chiapenec/Mephaa --------------------------------------------------------------------------------------
mangue_01 = vroom('../results/SweeD_selection/SweeD_Report.Mangue-01', skip = 2)
mangue_01$Chr = 1
mangue_01$SNP = paste(mangue_01$Chr, mangue_01$Position,sep = '_')

mangue_03 = vroom('../results/SweeD_selection/SweeD_Report.Mangue-03', skip = 2)
mangue_03$Chr = 3
mangue_03$SNP = paste(mangue_03$Chr, mangue_03$Position,sep = '_')

# matlatzinca/Mazahua ----------------------------------------------------------------------------------
matlatzinca_08 = vroom('../results/SweeD_selection/SweeD_Report.matlazinca-08', skip = 2)
matlatzinca_08$Chr = 8
matlatzinca_08$SNP = paste(matlatzinca_08$Chr, matlatzinca_08$Position,sep = '_')

#Yucatec individual chromosomes --------------------------------------------------------------
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
# figure plots --------------------------------------------------------------------------------
(yucatec_01_plot = ggplot(sweed_yucatec, aes(x = Position, y = Likelihood,
                                             color = Position < max(yucatec_highlight %>% 
                                                                      filter(CHR == 1) %>% 
                                                                      pull(BP)) + 1000000 & 
                                               Position > min(yucatec_highlight %>% 
                                                                filter(CHR == 1) %>% 
                                                                pull(BP)) - 1000000)
) +
    geom_scattermore(alpha = 1, pointsize = 3.2, pixels = c(2056, 1024)) + 
    scale_color_manual(values = c('black', 'red'), breaks = scales::pretty_breaks(n = 4)) +
    scale_x_continuous(n.breaks = 4, labels = unit_format(unit = "M", scale = 1e-6)) +
    labs(color='presence in GWAS') +
    ggtitle('Yucatan Mayan, chr 1') + 
    theme_bw() +
    theme(legend.position="none", 
          title = element_text(size = 8),
          axis.text.x = element_text(size = 7))
)

(otomi_04_plot = ggplot(otomi_04, aes(x = Position, y = Likelihood,
                                      color = Position < max(otomi_highlight %>% 
                                                               filter(CHR == 4) %>% 
                                                               pull(BP)) + 1000000 & 
                                        Position > min(otomi_highlight %>% 
                                                         filter(CHR == 4) %>% 
                                                         pull(BP)) - 1000000)
)+
    geom_scattermore(alpha = 1, pointsize = 3.2, pixels = c(2056, 1024)) + 
    scale_color_manual(values = c('black', 'red')) +
    scale_x_continuous(n.breaks = 4, labels = unit_format(unit = "M", scale = 1e-6)) +
    labs(color='presence in GWAS') +
    ggtitle('Otomi, chr 4') + 
    theme_bw() +
    theme(legend.position="none", 
          title = element_text(size = 8),
          axis.text.x = element_text(size = 7))
)

(otomi_09_plot = ggplot(otomi_09, aes(x = Position, y = Likelihood,
                                      color = Position < max(otomi_highlight %>% 
                                                               filter(CHR == 9) %>% 
                                                               pull(BP)) + 1000000 & 
                                        Position > min(otomi_highlight %>% 
                                                         filter(CHR == 9) %>% 
                                                         pull(BP)) - 1000000)
)+
    geom_scattermore(alpha = 1, pointsize = 3.2, pixels = c(2056, 1024)) + 
    scale_x_continuous(n.breaks = 4, labels = unit_format(unit = "M", scale = 1e-6)) +
    scale_color_manual(values = c('black', 'red')) +
    labs(color='presence in GWAS') +
    ggtitle('Otomi, chr 9') + 
    theme_bw() +
    theme(legend.position="none", 
          title = element_text(size = 8),
          axis.text.x = element_text(size = 7))
) 

(cora_04_plot = ggplot(cora_04, aes(x = Position, y = Likelihood,
                                    color = Position < min(cora_highlight %>% 
                                                             filter(CHR == 4) %>% 
                                                             pull(BP)) + 1000000 & 
                                      Position > min(cora_highlight %>% 
                                                       filter(CHR == 4) %>% 
                                                       pull(BP)) - 1000000)
)+
    geom_scattermore(alpha = 1, pointsize = 3.2, pixels = c(2056, 1024)) + 
    scale_color_manual(values = c('black', 'red'), breaks = scales::pretty_breaks(n = 5)) +
    scale_x_continuous(n.breaks = 4, labels = unit_format(unit = "M", scale = 1e-6)) +
    labs(color='presence in GWAS') +
    ggtitle('Cora/Huichol, chr 4') + 
    theme_bw() +
    theme(legend.position="none", 
          title = element_text(size = 8),
          axis.text.x = element_text(size = 7))
)




(all_sweed_plots = plot_grid(yucatec_01_plot + theme(legend.position="none"), otomi_04_plot, otomi_09_plot, cora_04_plot,
                             labels = 'AUTO', nrow = 4))
ggsave('all_sweed_plots.png', 
       plot = all_sweed_plots,
       device = 'png', 
       path = '../plots/', 
       width = 8, height = 8, units = 'in')

ggsave('all_sweed_plots.pdf', 
       plot = all_sweed_plots,
       device = 'pdf', 
       path = '../plots/', 
       width = 8, height = 8, units = 'in')

### create plot for Fst ########################################################
NL_fst = vroom('../results/NL/NL.weir.weir.fst')
haynie_fst = vroom('../results/Haynie/Haynie.weir.weir.fst')
NL_fst$SNP = paste(NL_fst$CHROM, NL_fst$POS,sep = '_')
haynie_fst$SNP = paste(haynie_fst$CHROM, haynie_fst$POS,sep = '_')

# manual_manhattan(NL_fst, chr_col = 'CHROM', sig = 0.01, 
#                  p_col = 'WEIR_AND_COCKERHAM_FST', bp_col = 'POS', title = 'FST', highlight_SNP_list = c())
# 
# manual_manhattan(haynie_fst, chr_col = 'CHROM', sig = 0.01, 
#                  p_col = 'WEIR_AND_COCKERHAM_FST', bp_col = 'POS', title = 'FST', highlight_SNP_list = c())

NL_fst$WEIR_AND_COCKERHAM_FST %>% mean()
haynie_fst$WEIR_AND_COCKERHAM_FST %>% mean()

max_bp = NL_fst %>% group_by(CHROM) %>% 
  summarise(max_bp = max(POS), center = floor(mean(POS))) %>% 
  mutate(bp_add = lag(cumsum(max_bp), default = 0),
         center_add = bp_add + center)
NL_fst_results = NL_fst %>% inner_join(max_bp, by = 'CHROM') %>% mutate(cumulative_bp = POS + bp_add)

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

haynie_fst_results = haynie_fst %>% inner_join(max_bp, by = 'CHROM') %>% mutate(cumulative_bp = POS + bp_add)
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


### check outliers #####################
haynie_fst %>% filter(WEIR_AND_COCKERHAM_FST > 0.25)
yucatec_gwas_results %>% filter(P < 1e-8 & CHR == 8) %>% pull(BP) %>% min()
## yucatec gwas chr 1 min at 165928841, max at 167679000
## yucatec gwas chr 2 min at 204174296, max at 204174757
## yucatec gwas chr 8 min at 27475338, max at 27475338

otomi_gwas_results %>% filter(P < 1e-8 & CHR == 4) %>% pull(BP) %>% min()
## otomi gwas chr 9 min at 38991555, max at 39215106
## otomi gwas chr 4 min at 178304566, max at 178304566

yucatec_highlight_plot = ggplot(sweed_yucatec, aes(x = Position, y = Likelihood, 
                                                   color = SNP %in% (sweed_yucatec %>% 
                                                                       filter(Position > 165928841 &
                                                                                Position < 167679000) %>% 
                                                                       pull(SNP)))) +
  geom_point() + 
  ggtitle('Yucatan Mayan, chr 1')

gene_ids = c('Zm00001eb030400,Zm00001eb030410,Zm00001eb030420,Zm00001eb030430,Zm00001eb030440,Zm00001eb030450,Zm00001eb030460,Zm00001eb030470,Zm00001eb030480,Zm00001eb030490,Zm00001eb030500,Zm00001eb030510,Zm00001eb030520,Zm00001eb030530,Zm00001eb030540,Zm00001eb030550,Zm00001eb030560,Zm00001eb030570,Zm00001eb030580,Zm00001eb030590,Zm00001eb030600,Zm00001eb030610,Zm00001eb030620,Zm00001eb030630,Zm00001eb030640,Zm00001eb030650,Zm00001eb030660,Zm00001eb030670,Zm00001eb030680,Zm00001eb030690,Zm00001eb030700,Zm00001eb030710,Zm00001eb030720,Zm00001eb030740,Zm00001eb030750'
             )

go_list = vroom('../results/GO_list.tsv')
go_list %>% select(`Gene > Gene ID`,
                   `Gene > GO Annotation > Ontology Term . Name`) %>% 
  count(`Gene > GO Annotation > Ontology Term . Name`) %>% 
  arrange(desc(n)) %>% print(n = 30)

go_list %>% filter(`Gene > GO Annotation > Ontology Term . Name` %in% c('zinc ion binding', 'metal ion binding')) %>% 
  select(`Gene > Gene ID`, `Gene > GO Annotation > Ontology Term . Name`) 

sweed_yucatec %>% arrange(desc(Likelihood)) %>% filter(abs(Position - 165849563) < 1000000)
yucatec_gwas_results %>% filter(CHR == 1 & P < 1e-7) %>% pull(BP)
sweed_yucatec %>% arrange(desc(Likelihood)) %>% View()

yucatec_01_plot + geom_point(data = sweed_yucatec %>% filter(Position == 166116864), aes(color = SNP)) + scale_color_manual(values = c('blue', 'green', 'orange', 'yellow')
)
          

getTopPercentageCLR = function(gwas_results, CHR, min_pos, max_pos, sweed_results) {
  top_gwas_results = gwas_results %>% filter(CHR == CHR & P < 1e-8) %>% pull(BP)
  sweed_results = sweed_results %>% filter(Chr == CHR)
  top_sweed_results = sweed_results %>% 
    arrange(desc(Likelihood))
  top_sweed_result = top_sweed_results %>%
    filter(Position > min_pos & 
             Position < max_pos) %>% 
    pull(Position) %>% 
    first()
  print(top_sweed_result)
  CLR_ind = which(top_sweed_results$Position == top_sweed_result)
  print(CLR_ind)
  return(CLR_ind / 40000 * 100)  
}          


# Sweed results CLR ranking -------------------------------------------------------------------
## Yucatec
getTopPercentageCLR(yucatec_gwas_results, 1, 165928841 - 500000, 167679000 + 500000, sweed_yucatec)
getTopPercentageCLR(yucatec_gwas_results, 2, 204174296 - 500000, 204174757 + 500000, sweed_yucatec)
getTopPercentageCLR(yucatec_gwas_results, 8, 27475338 - 500000, 27475338 + 500000, sweed_yucatec)
##Otomi
getTopPercentageCLR(otomi_gwas_results, 4, 178304566 - 500000, 178304566 + 500000, otomi_04)
getTopPercentageCLR(otomi_gwas_results, 5, 199268719 - 500000, 199268719 + 500000, otomi_05)
getTopPercentageCLR(otomi_gwas_results, 9, 38991589 - 500000, 38991589 + 500000, otomi_09)
## Cora/Huichol
# 10384903 - 500000, 10384903 + 500000
# getTopPercentageCLR(cora_gwas_results, 4, 2591046 - 50000, 2591046 + 50000, cora_04)
getTopPercentageCLR(cora_gwas_results, 4, 10384903 - 500000, 10384903 + 500000, cora_04)

## Matlatzinca
getTopPercentageCLR(matlatzinca_gwas_results, 8, 5366101 - 500000, 5366101 + 500000, matlatzinca_08)
getTopPercentageCLR(matlatzinca_gwas_results, 8, 29095461 - 500000, 29095461 + 500000, matlatzinca_08)

## Mangue
getTopPercentageCLR(mangue_gwas_results, 1, 28224363 - 500000, 28224363 + 500000, mangue_01)
getTopPercentageCLR(mangue_gwas_results, 3, 181819647 - 500000, 181819647 + 500000, mangue_03)

## Pipil
getTopPercentageCLR(pipil_gwas_results, 9, 148360542 - 500000, 148360542 + 500000, pipil_09)

## Chol
getTopPercentageCLR(chol_gwas_results, 5, 213908202 - 500000, 213908202 + 500000, chol_05)
getTopPercentageCLR(chol_gwas_results, 7, 158526307 - 500000, 158526307 + 500000, chol_07)

