library(tidyverse)
library(vroom)

language_families = c('mayan', 'aztecan', 'otomanguean')

for(language in language_families) {
	lan<-vroom(sprintf("/group/jrigrp10/maize-linguistics/results/%sHaynieMaster.csv", language), delim = ',') 
	lan_small<-select(lan,"Unique.ID", "Polygon Language Name")
	colnames(lan_small)<-c("Unique.ID","Polygon.Language")
	languages<-levels(factor(lan_small$Polygon.Language))
	new_dat<-data.frame(sapply(languages, function(i) as.numeric(lan_small$Polygon.Language==i)))
	new_new_dat<-cbind(lan_small,new_dat)
	write.table(new_new_dat,file=sprintf("/group/jrigrp10/maize-linguistics/data/%s/0_1_%s.txt", language, language),row.names=F,quote=F,col.names=T,sep="\t")

}

