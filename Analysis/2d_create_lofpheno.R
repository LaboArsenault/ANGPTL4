#!/usr/bin/env Rscript
#import instrument
library(data.table)
library(tidyverse)
library(furrr)
library(survival)
library(lubridate)
library(survminer)

setwd("/mnt/sda/gagelo01/Projects/small_MR_exploration/Test_potential_project/ANGPTL4/")

res_disease <- fread( "Data/Modified/res_disease.txt")
cadd_loff <- fread( "Data/Modified/cadd_loff.txt")
lof_individual <- fread( "Data/Raw/lof_individual.txt")
data_pheno <-fread( "Data/Modified/data_pheno.txt")
data_pca<-fread("Data/Modified/data_pca.txt")

dat <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "f.eid", all = TRUE), list(data_pheno, data_pca, res_disease, lof_individual))

dat[, date_birth := ymd(paste0(f.34.0.0, "-", f.52.0.0, "-01")) %>% as.IDate(.)]
dat[, age_enrollment := f.53.0.0 - date_birth]
colnom_date <- colnames(dat)[grepl("_date$", colnames(dat))]

date_end_followup <- dat[, lapply(.SD, function(x) max(x, na.rm = TRUE)), .SDcols = colnom_date]
date_end_followup <- round_date(as.IDate(date_end_followup[,apply(.SD, 1, max)]), unit = "month")+30 #Here is the date the follow up ends.

vec_arg <- gsub("_date", "", colnom_date)
start_followup<-"date_birth"

###toinclude
dat[,toinclude:=0]
dat[!is.na(PCA1),toinclude:=1]

#####test each variant individually with TG######
IV <- cadd_loff$lof_name
DV <- c("tg")
cov_inc<-""
# cov_inc <- paste(c("+ age_enrollment + sex", paste0("PCA", 1:10)), collapse = " + ")
arguments <- tidyr::expand_grid(IV, DV, cov_inc)
setDT(arguments)
res <- map(split(arguments, 1:arguments[,.N]), function(x) 
  lm(paste0(x$DV, " ~ ", x$IV, x$cov_inc), data = dat) %>%
    broom::tidy(.)) %>%
  rbindlist(., fill = TRUE)

lof_sign_tg <- res[grepl("^lof",term) & p.value<0.05, ]$term
fwrite(res, "Data/Modified/res_logit_tg.txt")
#######create the variable ANGPTL4_lof######
cadd_threshold <- 20
lof_cadd_thresh <- cadd_loff[cadd_score>cadd_threshold,lof_name]
tosel <- intersect(lof_sign_tg, lof_cadd_thresh)
dat[, E40K_carrier := apply(.SD, 1, function(x) any(x, na.rm = FALSE)),.SDcols = cadd_loff[rsid=="rs116843064", lof_name]]
dat[, has_lof := apply(.SD, 1, function(x) any(x, na.rm = FALSE)),.SDcols = tosel]
dat[,sum(has_lof, na.rm = TRUE)]
carrier_toremove <- setdiff(colnames(dat)[grepl("^lof", colnames(dat))], tosel)
dat[,  toinclude := ifelse(apply(.SD, 1, function(x) any(x, na.rm = FALSE))==TRUE, 0, toinclude),.SDcols = carrier_toremove]

######make sure that dat is ready for cox#######
map(1:length(vec_arg),  function(i) { 
  outcome_date <- paste0(vec_arg[i], "_date")
  censoring_date <- (paste0(vec_arg[i], "_censoring_date"))
  dat[,(censoring_date) := dplyr::if_else(is.na(date_death), ymd(date_end_followup) %>% as.IDate(.), date_death) %>%
        dplyr::if_else(!is.na(get(outcome_date)), get(outcome_date), .)]
  dat[, (paste0(vec_arg[i], "_censored")) := ifelse(is.na(get(outcome_date)), 0, 1)]
  
  if(start_followup=="date_birth"){dat[, (paste0(vec_arg[i], "_time")) := as.IDate(get(censoring_date)) - as.IDate(date_birth)]}
  if(start_followup=="date_assessement"){dat[, (paste0(vec_arg[i], "_time")) := as.IDate(get(censoring_date)) - as.IDate(f.53.0.0)]}
  dat[, (paste0(vec_arg[i], "_toinclude")):=toinclude]
  dat[get(paste0(vec_arg[i], "_toexclude"))==TRUE,(paste0(vec_arg[i], "_toinclude")):=0]
  dat[get(paste0(vec_arg[i], "_time"))<1,(paste0(vec_arg[i], "_toinclude")):=0]
})
#######
#####run the coxph######
include_E40K_ascarrier <- c(TRUE, FALSE)
IV <- "has_lof"
DV <- c(vec_arg)
cov_inc <- paste(c("+ age_enrollment + sex", paste0("PCA", 1:10)), collapse = " + ")
arguments <- tidyr::expand_grid(IV, DV, cov_inc, include_E40K_ascarrier)
setDT(arguments)
res <- map(split(arguments, 1:arguments[,.N]), function(x) {
  data<-data.table(dat)
  if(!(x$include_E40K_ascarrier)) {data[E40K_carrier==TRUE,(paste0(x$DV, "_toinclude")):=0]}
  res <- GagnonMR::run_coxph_wrapper(dat = data[get(paste0(x$DV, "_toinclude"))==1,], IV = x$IV, DV = x$DV, cov_inc = x$cov_inc)
  res[,include_E40K_ascarrier:=x$include_E40K_ascarrier]
  return(res)}) %>%  rbindlist(., fill = TRUE)

res[exposure=="has_lofTRUE", ]
fwrite(res, "Data/Modified/cox_res.R")
fwrite(dat, "Data/Modified/cox_data.R")

message("This script finished without errors")