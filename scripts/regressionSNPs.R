library(tidyverse)
library(data.table)
library(foreach)
library(doParallel)
library(vegan)

cores = detectCores()
print(cores)
registerDoParallel(cores - 1)

args<-commandArgs(T)
language = args[1]
dataset = args[2]
# language = 'otomanguean'
# dataset = 'Haynie'

language_dataset = paste0(language, dataset)
print(language_dataset)

if (grepl('Haynie', language_dataset)) {
  elevation_mat = read.table(sprintf('../results/Haynie/%s_elevationSigned_distances.csv', language_dataset), sep = ',', row.names = 1, header = T)
  admix_mat = read.table(sprintf('../results/Haynie/%s_admixProp_distances.csv', language_dataset), row.names = 1, sep = ',', header = T)
  geographic_mat = read.table(sprintf('../results/Haynie/%s_haversine_distances.csv', language_dataset), row.names = 1, sep = ',', header = T)
  language_mat = read.table(sprintf('../results/Haynie/%s_levenshtein_language_distances.csv', language_dataset), row.names = 1, sep = ',', header = T)
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

## get unique order of accessions in language matrix
accession_order = rownames(language_mat)

language_mat = as.matrix(language_mat)
elevation_mat = abs(as.matrix(elevation_mat))
geographic_mat = as.matrix(geographic_mat)
admix_mat = abs(as.matrix(admix_mat))

print('imported matrices')

language_pca = prcomp(language_mat, rank = 5)
elevation_pca = prcomp(elevation_mat, rank = 5)
admix_pca = prcomp(admix_mat, rank = 5)
geographic_pca = prcomp(geographic_mat, rank = 5)
pca_predictors = cbind(geographic_pca$rotation, elevation_pca$rotation, admix_pca$rotation, language_pca$rotation)
colnames(pca_predictors) = outer(c('PC1', 'PC2', 'PC3', 'PC4', 'PC5'), c('geographic', 'elevation', 'admix', 'language'), paste, sep = '_')


print("computed PCA")

language_pcoa = wcmdscale(d = language_mat,
                          k = 5,
                          eig = TRUE,
                          add = FALSE,
                          x.ret = FALSE
)
elevation_pcoa = wcmdscale(d = elevation_mat,
                          k = 5,
                          eig = TRUE,
                          add = FALSE,
                          x.ret = FALSE
)
admix_pcoa = wcmdscale(d = admix_mat,
                       k = 5,
                       eig = TRUE,
                       add = FALSE,
                       x.ret = FALSE
)
geographic_pcoa  = wcmdscale(d = geographic_mat,
                            k = 5,
                            eig = TRUE,
                            add = FALSE,
                            x.ret = FALSE
)

print('computed PCoA')

pcoa_predictors = cbind(geographic_pcoa$points, 
	elevation_pcoa$points, 
	admix_pcoa$points, 
	language_pcoa$points)
colnames(pcoa_predictors) = c('geographic_PCo1', 'geographic_PCo2', 'geographic_PCo3', 'geographic_PCo4', 'geographic_PCo5',
	'elevation_PCo1',
	'admix_PCo1',
	'language_PCo1', 'language_PCo2', 'language_PCo3', 'language_PCo4', 'language_PCo5') 

## get point values for each accession instead of distances
admixProp = read.table('../data/admixProp.GBSsamples.txt')
admixProp$Unique.ID = gsub('.MRG.4.', ':', admixProp$V1)
admixProp$Unique.ID = gsub('.D17PEACXX.3.', ':', admixProp$Unique.ID)
## get order
admixProp_sorted = admixProp[match(accession_order, admixProp$Unique.ID), ]

coordinates = read.table('/group/jrigrp10/maize-linguistics/data//feems/all_linguistics_coordinates.txt', header = T)
coordinates_sorted = coordinates[match(accession_order, coordinates$V1), c('locations_latitude', 'locations_longitude', 'locations_elevation')]

combined_predictors = cbind(coordinates_sorted, admixProp_sorted$V2, language_pcoa$points)
rownames(combined_predictors) = accession_order
colnames(combined_predictors) = c('latitude', 'longitude', 'elevation', 'admixProp', 
	'language_PCo1', 'language_PCo2', 'language_PCo3', 'language_PCo4', 'language_PCo5')
combined_predictors = as.matrix(combined_predictors)


genotypes = fread(sprintf('../data/%s/%s.%s.LD-pruned.dosage.regression.vcf', language, language, dataset), sep = '\t', header = T)
genotypes_mat = as.matrix(genotypes)

print("Read in raw genotypes")

# ## for test
# test_genotypes = genotypes_mat[1:10, ]

snp_results = foreach(snp = c(1:nrow(genotypes_mat)), .combine = bind_rows, .errorhandling = 'pass') %dopar% {
		result = tryCatch({

			if(snp %% 1000 == 0) {
				print(snp)
			}

			genotype = genotypes_mat[snp, ]

			## LLEA + language models
			language_point_results = summary(lm(genotype ~ combined_predictors))
			language_pca_results = summary(lm(genotype ~ pca_predictors))
			language_pcoa_results = summary(lm(genotype ~ pcoa_predictors))

			## LLEA models
			no_language_point_results = summary(lm(genotype ~ combined_predictors[, -c(5:9)]))
			no_language_pca_results = summary(lm(genotype ~ pca_predictors[, -c(16:20)]))
			no_language_pcoa_results = summary(lm(genotype ~ pcoa_predictors[,-c(8:12)]))
			
			},
			warning = function(war) {
				result = c(snp = snp, 
					r.squared = NA, 
					adj.r.squared = NA)
				},
			error = function(err) {
				results = c(snp = snp, 
					r.squared = NA, 
					adj.r.squared = NA)
				},
			finally = {
				results = c(snp = snp, 
					geo_admix_pca.r.squared = no_language_pca_results$r.squared,
					geo_admix_pcoa.r.squared = no_language_pcoa_results$r.squared,
					geo_admix_point.r.squared = no_language_point_results$r.squared, 

					language_pca.r.squared = language_pca_results$r.squared,
					language_pcoa.r.squared = language_pcoa_results$r.squared,
					language_point.r.squared = language_point_results$r.squared,

					geo_admix_pca.adj.r.squared = no_language_pca_results$adj.r.squared,
					geo_admix_pcoa.adj.r.squared = no_language_pcoa_results$adj.r.squared,
					geo_admix_point.adj.r.squared = no_language_point_results$adj.r.squared, 

					language_pca.adj.r.squared = language_pca_results$adj.r.squared,
					language_pcoa.adj.r.squared = language_pcoa_results$adj.r.squared,
					language_point.adj.r.squared = language_point_results$adj.r.squared

					)
				})
		return(results)
}

print("Finished regression")

# snp_results = foreach(snp = c(1:nrow(genotypes_mat)), .combine = bind_rows, .errorhandling = 'pass') %dopar% {
# 		result = tryCatch({

# 			if(snp %% 1000 == 0) {
# 				print(snp)
# 			}

# 			genotype = genotypes_mat[snp, ]

# 			lm_results = lm(genotype ~ pca_predictors)
# 			lm.summary <-summary(lm_results)

# 			no_geography_results = summary(lm(genotype ~ pca_predictors[,6:20]))
# 			no_elevation_results = summary(lm(genotype ~ pca_predictors[,c(1:5, 11:20)]))
# 			no_admix_results = summary(lm(genotype ~ pca_predictors[,c(1:10, 16:20)]))
# 			no_language_results = summary(lm(genotype ~ pca_predictors[,c(1:15)]))

# 			lm_pcoa_results = lm(genotype ~ pcoa_predictors)
# 			lm.pcoa.summary <-summary(lm_pcoa_results)

# 			no_geography_pcoa_results = summary(lm(genotype ~ pcoa_predictors[,-c(1:5)]))
# 			no_elevation_pcoa_results = summary(lm(genotype ~ pcoa_predictors[,-c(6)]))
# 			no_admix_pcoa_results = summary(lm(genotype ~ pcoa_predictors[,-c(7)]))
# 			no_language_pcoa_results = summary(lm(genotype ~ pcoa_predictors[,-c(8:12)]))

			
# 			},
# 			warning = function(war) {
# 				result = c(snp = snp, 
# 					r.squared = NA, 
# 					adj.r.squared = NA)
# 				},
# 			error = function(err) {
# 				results = c(snp = snp, 
# 					r.squared = NA, 
# 					adj.r.squared = NA)
# 				},
# 			finally = {
# 				results = c(snp = snp, 
# 					full.r.squared = lm.summary$r.squared, 
# 					full.adj.r.squared = lm.summary$adj.r.squared,
# 					no_geography.r.squared = no_geography_results$r.squared, 
# 					no_geography.adj.r.squared = no_geography_results$adj.r.squared,
# 					no_elevation.r.squared = no_elevation_results$r.squared, 
# 					no_elevation.adj.r.squared = no_elevation_results$adj.r.squared,
# 					no_admix.r.squared = no_admix_results$r.squared, 
# 					no_admix.adj.r.squared = no_admix_results$adj.r.squared,
# 					no_language.r.squared = no_language_results$r.squared, 
# 					no_language.adj.r.squared = no_language_results$adj.r.squared,

# 					full_pcoa.r.squared = lm.pcoa.summary$r.squared, 
# 					full_pcoa.adj.r.squared = lm.pcoa.summary$adj.r.squared,
# 					no_geography_pcoa.r.squared = no_geography_pcoa_results$r.squared, 
# 					no_geography_pcoa.adj.r.squared = no_geography_pcoa_results$adj.r.squared,
# 					no_elevation_pcoa.r.squared = no_elevation_pcoa_results$r.squared, 
# 					no_elevation_pcoa.adj.r.squared = no_elevation_pcoa_results$adj.r.squared,
# 					no_admix_pcoa.r.squared = no_admix_pcoa_results$r.squared, 
# 					no_admix_pcoa.adj.r.squared = no_admix_pcoa_results$adj.r.squared,
# 					no_language_pcoa.r.squared = no_language_pcoa_results$r.squared, 
# 					no_language_pcoa.adj.r.squared = no_language_pcoa_results$adj.r.squared

# 					)
# 				})
# 		return(results)
# }


write.csv(snp_results, sprintf('../results/regression_%s_point.csv', language_dataset), quote = F, row.names = F)

print("Wrote results")

# maize_genotypes<-read_tsv(args[2])
# maize_genotypes<-select(maize_genotypes, ID, contains("SEEDGWAS"))
# #transpose the genotypes so the SNPs are columns and individuals are rows and conserve names
# print(paste("Transposing genotypes...", Sys.time()))
# t_maize_genotypes<-t(select(maize_genotypes, contains("GWAS")))  %>% as_tibble(rownames = NA) %>% rownames_to_column()
# colnames(t_maize_genotypes)<-c("sample_id",maize_genotypes$ID)
# print(paste("Joining transposed genotypes to PCs...", Sys.time()))
# t_maize_genotypes<-mutate(t_maize_genotypes, sample_id = str_split(sample_id, "\\.", simplify = T)[,1]) %>% left_join(., maize_wHumanPCcentroids, by="sample_id")


