#!/usr/bin/env Rscript
#import instrument
library(data.table)
library(tidyverse)
library(furrr)
library(survival)
library(lubridate)
library(survminer)
setwd("/mnt/sda/gagelo01/Projects/small_MR_exploration/Test_potential_project/ANGPTL4/")
my.lowerthan <- function(k1, k2) {
  totransformNA<- is.na(k1)&is.na(k2)
  k1[is.na(k1)]<- Inf
  k2[is.na(k2)]<- Inf
  res <- k1 < k2
  res[totransformNA]<-NA
  return(res)
}
my.min <- function(x) ifelse( !all(is.na(x)), min(x, na.rm=T), NA)


#CAD
####
arg_cad <- data.table(code_inclusion_icd10 = list(c("I21", "I22", "I23", "I241", "I252", "I240", "I240", "I248", "I249",
                                                    paste0("I25", c(0:2, 5:9)))),
                      code_inclusion_opcs = list(c(paste0("K40", 1:4), paste0("K41", 1:4), paste0("K45", 1:5), paste0("K49", c(1:2,8:9), "K502", paste0("K75", c(1:4, 8, 9))))),
                      code_exclusion_icd10 = list(NULL),
                      code_exclusion_opcs = list(NULL),
                      outcome_name = "CAD")

arg_as <- data.table(code_inclusion_icd10 = list(c("I350", "I352")),
                     code_inclusion_opcs = list(NULL),
                     code_exclusion_icd10 = list(paste0("I0", 1:9)),
                     code_exclusion_opcs = list(NULL),
                     outcome_name = "AS")

arg_t2d <- data.table(code_inclusion_icd10 = list(paste0("E1",1:4)),
                     code_inclusion_opcs = list(NULL),
                     code_exclusion_icd10 = list(NULL),
                     code_exclusion_opcs = list(NULL),
                     outcome_name = "T2D")
arguments <- rbindlist(list(arg_cad, arg_as, arg_t2d))
arguments[,nrows:=Inf]

list_res <- map(split(arguments, 1:arguments[,.N]), function(x) {
  
  k_icd10 <- GagnonMR::ukb_format_HAR_OPCS(code_inclusion = x$code_inclusion_icd10[[1]], code_exclusion = x$code_exclusion_icd10[[1]],
                                           nrows = x$nrows)
  
  if(!is.null(x$code_inclusion_opcs[[1]])) {
    k_opcs <-  GagnonMR::ukb_format_HAR_OPCS(code_inclusion = x$code_inclusion_opcs[[1]], code_exclusion = x$code_exclusion_opcs[[1]],
                                             nrows = x$nrows, fieldid_date = "f.41282.0.0", fieldid_code = "f.41272.0.0")
    k_icd10 <- merge(k_icd10, k_opcs, by = "f.eid", all = TRUE)
    k_icd10[,date := apply(.SD, 1, my.min), .SDcols = c("outcome_date.y", "outcome_date.x")]
    k_icd10[,toexclude := my.lowerthan(as.IDate(exclude_date.x), as.IDate(date))|my.lowerthan(as.IDate(exclude_date.y), as.IDate(date))]
  } else {setnames(k_icd10, "outcome_date", "date")}
  
  
  k <- k_icd10[,.(f.eid, date, toexclude)]
  setnames(k, c("date", "toexclude"), paste0(x$outcome_name, "_", c("date", "toexclude")))
  return(k)
})


res_disease <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "f.eid", all = TRUE),
                      list_res)

fwrite(res_disease, "Data/Modified/res_disease.txt")
message("This script finished without errors")
