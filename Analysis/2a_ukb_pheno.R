#!/usr/bin/env Rscript
library(data.table)
library(tidyverse)
library(lubridate)
library(survminer)
library(survival)
library(furrr)
library(GagnonMR)
wd<-"/mnt/sda/gagelo01/Projects/small_MR_exploration/Test_potential_project/ANGPTL4/"
setwd(wd)
nrows<-Inf
ukbpath<-"/home/couchr02/Mendel_UKB/Source/Phenotype/February_2023_Update/ukb671338.tab"
columnname <- fread(ukbpath, nrows = 2)
columnname <- names(columnname)

###PCA###
pca<-fread("/home/couchr02/Mendel_UKB/Source/Phenotype/October_30_2019_Refresh/ukb38266.tab", nrows =2)
data_pca <- fread(input = "/home/couchr02/Mendel_UKB/Source/Phenotype/October_30_2019_Refresh/ukb38266.tab",
                  select = c(1, which(grepl("f.22009.", colnames(pca)))),
                  nrows = nrows)
setnames(data_pca, colnames(data_pca), gsub("f.22009.0.","PCA" ,colnames(data_pca), fixed = TRUE))
data_pca <- data_pca[,1:11]
fwrite(data_pca, "Data/Modified/data_pca.txt")
#Unique thing to change
fieldidok<-c( "I21" )
fidtoexclude<-character()
########################
#create the data_pheno data.table
fieldid <- c("f.34.0.0", "f.52.0.0", "f.53.0.0", "f.31.0.0",  "f.21001.0.0", "f.48.0.0",
             "f.21000.0.0", #ethnicity
             "f.189.0.0", #townsed deprivation index at recruitment
             "f.20116.0.0", #smoking status
             "f.1558.0.0", #alcohol intake frequency
             "f.30870.0", "f.30760.0", "f.30780.0", #Triglycerides, HDL, LDL
             "f.40000.0.0") #date of death
             
data_pheno <- fread(input = ukbpath,
                      select = c(1, which(grepl(paste(fieldid, collapse = "|"), columnname))),
                      nrows = nrows)
data_pheno <- data_pheno[2:.N,]

setnames(data_pheno, c("f.21001.0.0", "f.48.0.0", "f.21000.0.0", "f.31.0.0", "f.30870.0.0", "f.30760.0.0", "f.30780.0.0", "f.189.0.0",
                       "f.20116.0.0", "f.1558.0.0", "f.40000.0.0"),
         c("BMI", "WC", "eth", "sex", "tg", "hdl", "ldl", "townsend", "smoking", "alcohol", "date_death"))
data_pheno[eth %in% c(-3, -1), eth := NA] 
data_pheno[!is.na(eth),eth := substr(eth,1,1)]
data_pheno[, eth := factor(eth, levels = 1:6, labels = paste0("Ethnicity (", c("White", "Mixed", "Southeast_Asian", "Black", "Chinese", "Other"), ")"))] #1 = white
data_pheno[, smoking := smoking %>% ifelse(. %in% 1:2, 1, .) %>% ifelse(.%in%c(-3), NA, .)]
data_pheno[, alcohol :=  ifelse(alcohol %in% -3, NA, 6-alcohol)]
data_pheno[, date_birth := ymd(paste0(f.34.0.0, "-", f.52.0.0, "-01")) %>% as.IDate(.)]
data_pheno[, age_enrollment := f.53.0.0 - date_birth]
col <- c("smoking", "sex"); data_pheno[, (col):=lapply(.SD, as.logical), .SDcols = col]
colnom<-c("BMI", "WC","tg", "hdl", "ldl", "townsend")
data_pheno[,(colnom):=lapply(.SD, function(x) (x-mean(x, na.rm = TRUE))/sd(x, na.rm = TRUE)) ,.SDcols = colnom]
fwrite(data_pheno, "Data/Modified/data_pheno.txt")

message("This script finished without errors")
