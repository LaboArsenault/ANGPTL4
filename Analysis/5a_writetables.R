#!/usr/bin/env Rscript
library(tidyverse)
library(data.table)
library(GagnonMR)
library(openxlsx)
library(gtsummary)
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
res_multicis <- fread("Data/Modified/res_multicis_independent.txt")
cadd_loff <- fread("Data/Modified/cadd_loff.txt")
#######Baseline characteristics######
cox_data <- fread("Data/Modified/cox_data.R")
cox_data[,age:=age_enrollment/365.25]
baseline_char <- cox_data[toinclude==1, .(WC, BMI, alcohol, ldl, hdl, tg, sex, age,
                                          smoking, has_lof,CAD_censored, T2D_censored)]
baseline_char[,alcohol:=as.numeric(alcohol)]
baseline_char$alcohol %>% class
baseline_char[,has_lof := ifelse(has_lof, "ANGPTL4 LOF carriers", "ANGPTL4 LOF non-carriers")]
tbl <- tbl_summary(baseline_char, by = "has_lof", label = list(PRS_dis_12_1 = "Pancreatitis GRS, (SD)",sex = "Sex (Men)", 
                                                               age = "age, (years)", BMI="Body Mass Index, (Kg/m²)", 
                                                               WC="Waist circumference, (cm)", alcohol = "Alcohol intake, (Likert 0-5)", smoking = "Smoking status, yes",
                                                               sbp = "Systolic blood pressure, (mm Hg)", 
                                                               dbp = "Diastolic blood pressure, (mm Hg)", Cholesterol="Total-Cholesterol, (mmol/L)", ldl="LDL-Cholesterol, (mmol/L)",
                                                               hdl="HDL-Cholesterol, (mmol/L)", tg="Triglycerides, (mmol/L)",
                                                               blood_glucose="Blood glucose, (mmol/L)", 
                                                               hba1c="Glycated hemoglobin, (mmol/mol)", 
                                                               apob="apoB, (mmol/L)",
                                                               alt = "Alanine aminotransferase, (U/L)", c_reactive_protein = "C-reactive protein, (mg/L)",
                                                               cad="Coronary Artery Disease (CAD) before baseline"),
                   type= list(alcohol="continuous"), missing="no") %>% 
  add_p()

baseline_xlsx <- tbl %>% gtsummary::as_tibble()
##############
dataset <- df_index[id %in% unique(c(data_res$id.outcome)), ]
dataset[, c("group_name", "author", "unit", "nsnp", "initial_build", "category", "sd", "note", "rawfile_wget", "rawfile_path_server" ) := NULL]
dataset[id=="dis-2-1", rawfile_path_web := "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90091001-GCST90092000/GCST90091033/GCST90091033_buildGRCh37.tsv.gz"]
dataset[id=="trait-1-1", rawfile_path_web := "https://zenodo.org/records/1251813/files/bmi.giant-ukbb.meta-analysis.combined.23May2018.txt.gz?download=1"]
dataset[id=="trait-10-1", rawfile_path_web := "https://zenodo.org/records/1251813/files/whradjbmi.giant-ukbb.meta-analysis.combined.23May2018.txt.gz?download=1"]
dataset[id=="trait-10-3", rawfile_path_web := "https://zenodo.org/records/1251813/files/whr.giant-ukbb.meta-analysis.combined.23May2018.txt.gz?download=1"]
dataset[id=="trait-6-1", rawfile_path_web := "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST008001-GCST009000/GCST008598/Results_90th_percentile.txt.gz"]
dataset[grepl("trait-13-", id), rawfile_path_web := "http://ldsc.broadinstitute.org/ldhub/"]

  
suptab2 <- data_res
suptab2 <- suptab2[SNP=="rs116843064", ]
suptab2[,c("mr_keep.outcome", "pval_origin.outcome", "action", "mr_keep", "remove", "palindromic", "ambiguous") := NULL]

res_tg_sd <- fread("Data/Modified/res_tg_sd.txt")
res_tg_sd[, hgnc:= "E40K"]
dt_pan <- res_multicis[grepl("trait-16-4", id.exposure) & method == "Inverse variance weighted", ]
k1<- separate(dt_pan[grepl("trait-16-4", id.exposure), ], col = "id.exposure", into = c("hgnc", "id.exposure"), sep = "_", remove = TRUE)
k1 <- rbind(k1, res_tg_sd, fill = TRUE)
k1[, c("type_of_test", "lci", "uci","fstat","rsq", "steiger_dir", "steiger_pval") := NULL]
suptab3 <- k1
scatter <- dcast(k1[method%in%c("Inverse variance weighted",  "Wald ratio") & !grepl("dis-15-", id.outcome),], id.outcome ~ hgnc, value.var = c("b", "se"))
colnom<-colnames(scatter)[grepl("^b_", colnames(scatter))]
suptab4 <- cor(scatter[,.SD,.SDcols = colnom], method = c("pearson"), use="complete.obs")

cadd_loff



list_supdat <- list( "Supplementary Table 1" = dataset,
                     "Supplementary Table 2" = suptab2,
                     "Supplementary Table 3" = suptab3,
                     "Supplementary Table 4" = suptab4,
                     "Supplementary Table 5" = cadd_loff,
                     "Supplementary Table 6" = baseline_xlsx,
                     "Supplementary Table 7" = res_cophescan)



#

dt_title <- data.table(title = paste0("ST", c(1:7)),
                       caption = c( "Description of the datasets used for Genetic analysis.",
                                    "Effect of E40K on the human phenome",
                                    "Effect of LPL enhancement, LIPC enhancement and ANGPTL4 inhibition on cardiometabolic traits and diseases. (all associations are scaled to a 1-SD change in triglycerides)",
                                    "Correlation matrix between E40K, ANGPTL4 inhibition, LIPC enhancement, and LPL enhancement",
                                    "LOF variants included",
                                    "Baseline characteristics of the UKB cohort",
                                    "Results of the cophescan analyses"))

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

col_description[[4]] <- tribble(
  ~x, ~y,
  "b_ANGPTL4", "ANGPTL4 inhibition proxied using 22 SNPs in the ANGPTL4 gene region (window 300Kb) associated with TG",
  "b_E40K", "ANGPTL4 inhibition proxied using E40K",
  "b_LIPC", "LIPC enhancement proxied with 71 independent SNPs in the LPL region +/- 300 Kb significantly associated with triglycerides",
  "b_HL", "HL enhancement was proxied with 23 SNPs in the region of the gene encoding HL significantly associated with triglyceride levels"
) %>% as.data.table(.)

col_description[[5]] <- tribble(
  ~x, ~y,
  "lof_name", "lof name",
  "cadd_score", "CADD score phred",
  "hgnc", "hgnc",
  "rsid", "rsid"
) %>% as.data.table(.)

col_description[[6]] <- tribble(
  ~x, ~y,
  "", "",
) %>% as.data.table(.)

col_description[[7]] <- tribble(
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

