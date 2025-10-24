library(tidyverse)
library(data.table)

## author: Forrest Li
## script to generate PCA based on SeeD GBS data - in this case, LD pruned SNPs used for equivalent linguistic regression analysis

args<-commandArgs(T)
language = args[1]
dataset = args[2]
# language = 'otomanguean'
# dataset = 'Haynie'

language_dataset = paste0(language, dataset)
print(language_dataset)


## use genotypes used for equivalent population of maize for specific language and dataset
genotypes = fread(sprintf('../data/%s/%s.%s.LD-pruned.dosage.regression.vcf', language, language, dataset), sep = '\t', header = T)
genotypes_mat = as.matrix(genotypes)

## same imputation method as in human PCs
impute_mean = function(X) {
  X_mean = colMeans(X,na.rm=T)
  X[is.na(X)] = X_mean[rep(1:ncol(X),colSums(is.na(X)))]
  X
}

## transpose

genotypes_t = t(genotypes_mat)

genotypes_impute = impute_mean(genotypes_t)
genotype.pca = prcomp(genotypes_impute, center = T, scale = F)
## finished pca

genotype.summary = summary(genotype.pca)

## get % of variance explained by top eigenvalues of prcomp
print(genotype.summary$importance[2,1:10])

write.table(enframe(genotype.summary$importance[2,1:10]) , 
	sprintf('../results/%s_%s_pca_eigens.txt', language, dataset),
	quote = F,
	sep = '\t',
	row.names = F)