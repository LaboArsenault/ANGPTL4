#!/usr/bin/env Rscript
library(tidyverse)
library(data.table)
library(GagnonMR)
library(openxlsx)
wd<-"/mnt/sda/gagelo01/Projects/small_MR_exploration/Test_potential_project/ANGPTL4/"
setwd(wd)
df_index <- data.table::fread("/mnt/sda/gagelo01/Vcffile/server_gwas_id.txt")
all_outcome <- fread("Data/Modified/all_outcome.txt" )

####align on the good reference allele
k <- all_outcome[id.outcome == "trait-16-4", ] %>% GagnonMR::convert_outcome_to_exposure(.)
harm <- TwoSampleMR::harmonise_data(k, all_outcome)
tosel <- colnames(harm)[!grepl("exposure", colnames(harm))]
data_res <- harm[,tosel] %>% as.data.table(.)
res_cophescan <- fread("Data/Modified/res_cophescan_hyperpriors.txt")
res_tg_sd <- fread("Data/Modified/res_tg_sd.txt")


####
dataset <- df_index[id %in% unique(c(data_res$id.outcome)), ]
dataset[, c("group_name", "author", "consortium", "unit", "nsnp", "initial_build", "category", "sd", "note" ) := NULL]

suptab2 <- data_res
suptab2 <- suptab2[SNP=="rs116843064", ]
suptab2[,c("mr_keep.outcome", "pval_origin.outcome", "action", "mr_keep", "remove", "palindromic", "ambiguous") := NULL]

suptab4<-res_tg_sd
suptab4 <- suptab4[method %in% c("Inverse variance weighted", "Wald ratio")]

list_supdat <- list( "Supplementary Table 1" = dataset,
                     "Supplementary Table 2" = suptab2,
                     "Supplementary Table 3" = res_cophescan,
                     "Supplementary Table 4" = suptab4)



#

dt_title <- data.table(title = paste0("ST", c(1:4)),
                       caption = c( "Description of the datasets used for Genetic analysis.",
                                    "Effect of E40K on the human phenome",
                                    "Results of the cophescan analyses",
                                    "Effect of LPL enhancement, LIPC enhancement and ANGPTL4 inhibition on cardiometabolic traits and diseases. (all associations are scaled to a 1-SD change in triglycerides)"))

###
col_description<- vector(mode = "list", length = length(list_supdat))
col_description[[1]] <- tribble(
  ~x, ~y,
  "id", "a unique id",
  "trait", "a unique name",
  "year", "the year the GWAS was published",
  "trait", "The author of the GWAS",
  "consortium", "the name of the consortium of the GWAS",
  "sex", "sex included",
  "population", "ancestry",
  "sample_size", "the sample_size",
  "pmid", "the pubmed ID",
  "ncase", "the number of cases",
  "ncontrol", "the number of controls",
  "adjustments", "what variables were included in the GWAS adjustment model",
  "clean_variable_name", "A publication ready name that can be used to plot figures."
) %>% as.data.table(.)

col_description[[2]] <- tribble(
  ~x, ~y,
  "id.outcome", "a unique id",
  "chr.outcome", "chromosome",
  "pos.outcome", "Position build 37",
  "other_allele.outcome", "The non-effect allele",
  "effect_allele.outcome", "The effect allele",
  "beta.outcome", "beta effect estimate",
  "se.outcome", "standard error of the estimate",
  "pval.outcome", "p-valueor of the estimate",
  "eaf.outcome", "effect allele frequency",
  "samplesize.outcome", "sample size",
  "ncase.outcome", "number of cases",
  "SNP", "rsid",
  "ncontrol.outcome", "number of controls",
  "outcome", "A unique name for the outcome",
  "clean_variable_name", "A publication ready name that can be used to plot figures."
) %>% as.data.table(.)

col_description[[3]] <- tribble(
  ~x, ~y,
  "PP.Hn", "posterio probability that there is no causal SNP",
  "PP.Ha", "posterior probability that there is a causal SNP distinct from the query variant",
  "PP.Hc", "posterior probability that the query variant is causal",
  "nsnps", "the number of SNPs",
  "querysnp", "The query SNP (in our case always rs116843064)",
  "querytrait", "the unique id of the trait",
  "typeBF", "the type of test either SuSiE or ABF using coloc",
  "grp", "the unique group",
  "cophe.hyp.call", "the hypothesis prioritised by cophescan",
) %>% as.data.table(.)

col_description[[4]] <- tribble(
  ~x, ~y,
  "id.outcome", "a unique id for the outcome",
  "id.exposure", "a unique id for the exposure",
  "outcome", "a unique name for the outcome",
  "exposure", "a unique name for twosample MR",
  "method", "The method used for the cochran's Q",
  "nsnp", "the number of SNPs used as IVs",
  "b", "beta",
  "se", "standard error",
  "pval", "p-value",
  "lci", "lower confidence interval",
  "uci", "upper confidence interval",
  "hgnc", "the hugo gene name nomenclature of the exposure",
) %>% as.data.table(.)

bold_st <- createStyle(textDecoration = "Bold")
wb <- createWorkbook()
for(i in 1:length(list_supdat)) {
  addWorksheet(wb, sheetName =  dt_title[i, title])

  title <- dt_title[i,paste0(title, " : ", caption)]
  writeData(wb, sheet = i, x = title, startCol = 1, startRow = 1)
  addStyle(wb, sheet = i, style = bold_st, rows = 1, cols = 1:2)
  writeData(wb, sheet = i, x = col_description[[i]], startCol = 1, startRow = 2)
  addStyle(wb, sheet = i, style = bold_st, rows = 2:col_description[[i]][,.N+2], cols = 1)
  deleteData(wb, sheet = i, rows = 2, cols = 1:2, gridExpand = TRUE)
  writeData(wb, sheet = i, x = list_supdat[[i]], startCol = 1, startRow = col_description[[i]][,.N+4])
  addStyle(wb, sheet = i, style = bold_st, rows = col_description[[i]][,.N+4], cols = 1:length(colnames(list_supdat[[i]])), gridExpand = TRUE, stack = TRUE)
}
saveWorkbook(wb, file = "Results/supplementary_tables.xlsx", overwrite = TRUE)

message("This script finished without errors")

