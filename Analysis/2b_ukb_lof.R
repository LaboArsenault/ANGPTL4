#!/usr/bin/env Rscript
library(data.table)
library(tidyverse)
wd<-"/mnt/sda/gagelo01/Projects/small_MR_exploration/Test_potential_project/ANGPTL4/"
setwd(wd)
nrows<-Inf#Inf
ukbpath<-"/home/couchr02/Mendel_UKB/Source/Phenotype/February_2023_Update/ukb671338.tab"
columnname <- fread(ukbpath, nrows = 2)
columnname <- names(columnname)
data_pca <- fread("Data/Modified/data_pca.txt")

##############Create lof_individual##########
wd<-"/mnt/sda/boujer01/ELOI/ANGPTL4_PTVs/ANGPTL4_individual/"
k <- list.files(wd)
res_list <- map(1:length(k), function(i) { 
data <- fread(paste0(wd, k[i]))
data[,ptv:=1]
data <- merge(data_pca[!is.na(PCA1),.(f.eid)], data, by.x ="f.eid", by.y = "column",  all.x = TRUE)
data[is.na(ptv),ptv:=0]
setnames(data, "ptv", gsub("_sample.csv|_sample.miss.csv", "", k[i]))
return(data)
})
lof_individual <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "f.eid"),
       res_list)
colnom<- setdiff(colnames(lof_individual), "f.eid")
setnames(lof_individual, colnom, paste0("lof_", colnom))
fwrite(lof_individual, "Data/Raw/lof_individual.txt")

########create lof_cadd#######
k <-fread("/mnt/sda/boujer01/ELOI/ANGPTL4_PTVs/VEP/ANGPTL4_VEP_all.txt")
k <- k[SYMBOL=="ANGPTL4",.(Location, Allele, CADD_PHRED, SYMBOL, Existing_variation)] %>% distinct
k$rsid <- map_chr(strsplit(k$Existing_variation, split = ","), function(x) x[1])
cadd <- separate(k, col = "Location", into = c("chr", "posmin", "posmax"), sep = ":|-")

k<- colnames(lof_individual)
k <- setdiff(k, "f.eid")
k <- data.table(lof_name=k)
k<- separate(k, col = "lof_name", into = c("todump","chr", "pos", "nea", "ea"), sep = "_", remove = FALSE)
cadd_loff <- merge(k, cadd, by.x = c("chr", "pos"), by.y = c("chr",  "posmin"), all.x = TRUE)
cadd_loff <- cadd_loff[,.(lof_name, CADD_PHRED, SYMBOL, rsid)]
setnames(cadd_loff, c("CADD_PHRED", "SYMBOL"), c("cadd_score", "hgnc"))
fwrite(cadd_loff, "Data/Modified/cadd_loff.txt")
message("This script finished without errors")
