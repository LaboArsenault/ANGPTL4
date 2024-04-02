#!/usr/bin/env Rscript
setwd("/mnt/sda/gagelo01/Projects/small_MR_exploration/Test_potential_project/ANGPTL4")
mysource <- "tosource.R"

library(GagnonMR)
myArray <- GagnonMR::list_masterscript("/mnt/sda/gagelo01/Projects/Pipelines/Robust_efficient_drug_pipeline/Analysis/")
myArray <- myArray[c(2:5, 8)]
GagnonMR::write_masterscript(myArray = myArray, mysource = mysource, output_name = "Analysis/drugpipeline.sh")
message("This script finished without errors")