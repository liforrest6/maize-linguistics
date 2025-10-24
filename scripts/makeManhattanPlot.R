library(ggplot2)
library(data.table)
library(dplyr)
library(stringr)
library(ggtext)

## author: Forrest Li
## script to generate automated manhattan plots for each linguistic GWAS


file_name = as.character(commandArgs(t=T))[1]

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
      geom_point(alpha = 0.75, size = 0.5) +
      scale_x_continuous(
        label = axis_set %>% pull((!!sym(chr_col))),
        breaks = axis_set$center
      ) +
      scale_y_continuous(expand = c(0, 0), limits = c(0, ylim)) +
      scale_color_manual(values = rep(
        c("grey4", "grey30"),
        unique(length(axis_set %>% pull(!!sym(chr_col))))
      )) +
      geom_point(data = results_add[results_add$SNP %in% highlight_SNP_list,],
                 alpha = 0.75, size = 0.75, color = 'red') + 
      scale_size_continuous(range = c(0.5, 3)) +
      labs(
        x = NULL,
        y = "-log<sub>10</sub>(p)"
      ) +
      theme_minimal() +
      theme(
        legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title.y = element_markdown(),
        axis.text.x = element_text(angle = 0, size = 10, vjust = 0.5)
      ) +
      ggtitle(title)
  )
  
}

file_path = str_split_1(file_name, '/')
grouping = str_split_1(file_path[length(file_path)], '_')[1]
results = read.table(file_name, sep = '\t', header = T)
p = manual_manhattan(results, highlight_SNP_list = list(), title = grouping)

# do.call(file.path, as.list(file_path[-c(length(file_path))]), grouping, '_manhattan.png')

png(paste0(do.call(file.path, as.list(file_path[-c(length(file_path))])), '/', grouping, '_manhattan.png'), width = 6, height = 3, units = 'in', res = 400)
plot(p)
dev.off()


