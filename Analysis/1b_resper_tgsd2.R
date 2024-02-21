#!/usr/bin/env Rscript
setwd("/home/gagelo01/workspace/Projects/small_MR_exploration/Test_potential_project/ANGPTL4")
source("Analysis/tosource.R")
source("/mnt/sda/gagelo01/Projects/small_MR_exploration/F2_finngen/Analysis/func_F2.R")

out <- fread("Data/Modified/all_outcome.txt" )

inst <- out[id.outcome=="trait-16-4", ]
colnom<-colnames(inst)[grepl("outcome", colnames(inst))]
setnames(inst, colnom, gsub("outcome", "exposure", colnom))
# harm <- TwoSampleMR::harmonise_data(inst[SNP=="rs116843064",],
#                             outcome_dat = out[!grepl("dis-15", id.outcome),],
#                             action = 1)

res <- wald_quick(inst_all_sign_clump_arg = inst[SNP=="rs116843064",], all_outcome_arg = out[!grepl("dis-15", id.outcome),], 
                  id_exposure_name = "trait-16-4")

fwrite(res, "Data/Modified/res_tg_sd.txt")
message("This script finished without yours")
