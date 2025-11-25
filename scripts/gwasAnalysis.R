# GWAS and selection scan analysis

source('languageFunctions.R')


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

# write.table(yucatec_highlight$SNP, '../results/mayan_gemma_output/yucatec_highlight.txt', sep = '\t', quote = F, row.names = F, col.names =F)
# write.table(otomi_highlight$SNP, '../results/otomanguean_gemma_output/otomi_highlight.txt', sep = '\t', quote = F, row.names = F, col.names =F)
# write.table(cora_highlight$SNP, '../results/aztecan_gemma_output/cora_highlight.txt', sep = '\t', quote = F, row.names = F, col.names =F)
# write.table(chol_highlight$SNP, '../results/mayan_gemma_output/chol_highlight', sep = '\t', quote = F, row.names = F, col.names =F)
# write.table(matlatzinca_highlight$SNP, '../results/otomanguean_gemma_output/matlatzinca_highlight', sep = '\t', quote = F, row.names = F, col.names =F)
# write.table(mangue_highlight$SNP, '../results/otomanguean_gemma_output/mangue_highlight', sep = '\t', quote = F, row.names = F, col.names =F)
# write.table(pipil_highlight$SNP, '../results/aztecan_gemma_output/pipil_highlight', sep = '\t', quote = F, row.names = F, col.names =F)



# read SweeD files ----------------------------------------------------------------------------

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

# Otomi 
otomi_04 = vroom('../results/SweeD_selection/SweeD_Report.Otomi-04', skip = 2)
otomi_04$Chr = 4
otomi_04$SNP = paste(otomi_04$Chr, otomi_04$Position, sep = '_')

otomi_05 = vroom('../results/SweeD_selection/SweeD_Report.Otomi-05', skip = 2)
otomi_05$Chr = 5
otomi_05$SNP = paste(otomi_05$Chr, otomi_05$Position, sep = '_')

otomi_09 = vroom('../results/SweeD_selection/SweeD_Report.Otomi-09', skip = 2)
otomi_09$Chr = 9
otomi_09$SNP = paste(otomi_09$Chr, otomi_09$Position, sep = '_')

# Cora/Huichol
cora_04 = vroom('../results/SweeD_selection/SweeD_Report.Cora-Huichol-04', skip = 2)
cora_04$Chr = 4
cora_04$SNP = paste(cora_04$Chr, cora_04$Position,sep = '_')

# Chol/Chontal/Chorti/Cholti
chol_05 = vroom('../results/SweeD_selection/SweeD_Report.Chol-05', skip = 2)
chol_05$Chr = 5
chol_05$SNP = paste(chol_05$Chr, chol_05$Position,sep = '_')

chol_07 = vroom('../results/SweeD_selection/SweeD_Report.Chol-07', skip = 2)
chol_07$Chr = 7
chol_07$SNP = paste(chol_07$Chr, chol_07$Position,sep = '_')

# Pipil/Eastern + Northern Puebla Nahuatl
pipil_09 = vroom('../results/SweeD_selection/SweeD_Report.Pipil-09', skip = 2)
pipil_09$Chr = 9
pipil_09$SNP = paste(pipil_09$Chr, pipil_09$Position,sep = '_')

# Mangue/Chiapenec/Mephaa 
mangue_01 = vroom('../results/SweeD_selection/SweeD_Report.Mangue-01', skip = 2)
mangue_01$Chr = 1
mangue_01$SNP = paste(mangue_01$Chr, mangue_01$Position,sep = '_')

mangue_03 = vroom('../results/SweeD_selection/SweeD_Report.Mangue-03', skip = 2)
mangue_03$Chr = 3
mangue_03$SNP = paste(mangue_03$Chr, mangue_03$Position,sep = '_')

# matlatzinca/Mazahua 
matlatzinca_08 = vroom('../results/SweeD_selection/SweeD_Report.matlazinca-08', skip = 2)
matlatzinca_08$Chr = 8
matlatzinca_08$SNP = paste(matlatzinca_08$Chr, matlatzinca_08$Position,sep = '_')

# Sweed results CLR ranking -------------------------------------------------------------------
## Yucatec
getTopPercentageCLR(yucatec_gwas_results, 1, 165928841 - 500000, 167679000 + 500000, sweed_yucatec)
getTopPercentageCLR(yucatec_gwas_results, 2, 204174296 - 500000, 204174757 + 500000, sweed_yucatec)
getTopPercentageCLR(yucatec_gwas_results, 8, 27475338 - 500000, 27475338 + 500000, sweed_yucatec)

## Otomi
getTopPercentageCLR(otomi_gwas_results, 4, 178304566 - 500000, 178304566 + 500000, otomi_04)
getTopPercentageCLR(otomi_gwas_results, 5, 199268719 - 500000, 199268719 + 500000, otomi_05)
getTopPercentageCLR(otomi_gwas_results, 9, 38991589 - 500000, 38991589 + 500000, otomi_09)

## Cora/Huichol
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



# Fst analysis --------------------------------------------------------------------------------

## read Fst files
NL_fst = vroom('../results/NL/NL.weir.weir.fst')
haynie_fst = vroom('../results/Haynie/Haynie.weir.weir.fst')
NL_fst$SNP = paste(NL_fst$CHROM, NL_fst$POS,sep = '_')
haynie_fst$SNP = paste(haynie_fst$CHROM, haynie_fst$POS,sep = '_')

## generate BP counter for manhattan plots
max_bp = NL_fst %>% group_by(CHROM) %>% 
  summarise(max_bp = max(POS), center = floor(mean(POS))) %>% 
  mutate(bp_add = lag(cumsum(max_bp), default = 0),
         center_add = bp_add + center)

## add column for BP counter
NL_fst_results = NL_fst %>% inner_join(max_bp, by = 'CHROM') %>% mutate(cumulative_bp = POS + bp_add)
haynie_fst_results = haynie_fst %>% inner_join(max_bp, by = 'CHROM') %>% mutate(cumulative_bp = POS + bp_add)

### check Fst outliers and enrichment #####################
haynie_fst %>% filter(WEIR_AND_COCKERHAM_FST > 0.25)
yucatec_gwas_results %>% filter(P < 1e-8 & CHR == 8) %>% pull(BP) %>% min()
## yucatec gwas chr 1 min at 165928841, max at 167679000
## yucatec gwas chr 2 min at 204174296, max at 204174757
## yucatec gwas chr 8 min at 27475338, max at 27475338

otomi_gwas_results %>% filter(P < 1e-8 & CHR == 4) %>% pull(BP) %>% min()
## otomi gwas chr 9 min at 38991555, max at 39215106
## otomi gwas chr 4 min at 178304566, max at 178304566

gene_ids = c('Zm00001eb030400,Zm00001eb030410,Zm00001eb030420,Zm00001eb030430,Zm00001eb030440,Zm00001eb030450,Zm00001eb030460,Zm00001eb030470,Zm00001eb030480,Zm00001eb030490,Zm00001eb030500,Zm00001eb030510,Zm00001eb030520,Zm00001eb030530,Zm00001eb030540,Zm00001eb030550,Zm00001eb030560,Zm00001eb030570,Zm00001eb030580,Zm00001eb030590,Zm00001eb030600,Zm00001eb030610,Zm00001eb030620,Zm00001eb030630,Zm00001eb030640,Zm00001eb030650,Zm00001eb030660,Zm00001eb030670,Zm00001eb030680,Zm00001eb030690,Zm00001eb030700,Zm00001eb030710,Zm00001eb030720,Zm00001eb030740,Zm00001eb030750'
)

go_list = vroom('../results/GO_list.tsv')

# assess enrichment
go_list %>% select(`Gene > Gene ID`,
                   `Gene > GO Annotation > Ontology Term . Name`) %>% 
  count(`Gene > GO Annotation > Ontology Term . Name`) %>% 
  arrange(desc(n)) %>% print(n = 30)

go_list %>% filter(`Gene > GO Annotation > Ontology Term . Name` %in% c('zinc ion binding', 'metal ion binding')) %>% 
  select(`Gene > Gene ID`, `Gene > GO Annotation > Ontology Term . Name`) 
