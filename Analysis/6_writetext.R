#!/usr/bin/env Rscript
library(tidyverse)
library(data.table)
library(GagnonMR)
wd<-"/mnt/sda/gagelo01/Projects/small_MR_exploration/Test_potential_project/ANGPTL4/"
setwd(wd)
source("Analysis/tosource.R")

all_outcome <- fread("Data/Modified/all_outcome.txt" )
df_index <- fread("/mnt/sda/gagelo01/Vcffile/server_gwas_id.txt")

#####
k <- all_outcome[id.outcome == "trait-16-4", ] %>% GagnonMR::convert_outcome_to_exposure(.)
harm <- TwoSampleMR::harmonise_data(k, all_outcome)
tosel <- colnames(harm)[!grepl("exposure", colnames(harm))]
data_res <- harm[,tosel] %>% as.data.table(.)
data_res[, lci.outcome := beta.outcome-1.96*se.outcome]
data_res[, uci.outcome := beta.outcome+1.96*se.outcome]
data_res <- merge(data_res, df_index[,.(id,clean_variable_name)], by.x = "id.outcome", by.y = "id")
cox_data <- fread( "Data/Modified/cox_data.R")
cox_res <- fread("Data/Modified/cox_res.R")
cadd_loff <- fread( "Data/Modified/cadd_loff.txt")
res_logit_tg <- fread( "Data/Modified/res_logit_tg.txt")
res_tg_sd <- fread("Data/Modified/res_tg_sd.txt")
res_tg_sd[, hgnc:= "E40K"]
dt_pan <- fread( "Data/Modified/res_multicis_independent.txt")
k1<- separate(dt_pan[grepl("trait-16-4", id.exposure), ], col = "id.exposure", into = c("hgnc", "id.exposure"), sep = "_", remove = TRUE)
k1 <- rbind(k1, res_tg_sd, fill = TRUE)
scatter <- dcast(k1[method%in%c("Inverse variance weighted",  "Wald ratio") & !grepl("dis-15-", id.outcome),], id.outcome ~ hgnc, value.var = c("b", "se"))
res_cophescan <- fread("Data/Modified/res_cophescan_hyperpriors.txt")
res_hyprcoloc <-fread("Data/Modified/res_hyprcoloc.txt")
# inst_tgsd <- fread("Data/Modified/inst_tgsd.txt")
######
return_format_data<-function(data) {
  k <- data[, paste0("OR = ", format(round(exp(beta.outcome), digits = 2), nsmall = 2), ", 95% CI=", format(round(exp(lci.outcome), digits = 2), nsmall = 2), " to ",  format(round(exp(uci.outcome), digits = 2), nsmall = 2), ", p=",pval.outcome %>% formatC(., format = "e", digits = 1))]
  names(k) <- data[,paste0(outcome)]
  return(k)
}
return_format_data_noexp <-function(data) {
  k<-data[, paste0(format(round(beta.outcome, digits = 2), nsmall = 2) , " 95% CI=", format(round(lci.outcome, digits = 2), nsmall = 2), " to ",  format(round(uci.outcome, digits = 2), nsmall = 2), ", p=",pval.outcome %>% formatC(., format = "e", digits = 1))]
  names(k) <- data[,paste0(outcome)]
  return(k)
}

return_format_data_noexp_HR <-function(data) {
  k<-data[, paste0("HR = ", format(round(HR, digits = 2), nsmall = 2) , " 95% CI=", format(round(lci, digits = 2), nsmall = 2), " to ",  format(round(uci, digits = 2), nsmall = 2), ", p=",pval %>% formatC(., format = "e", digits = 1))]
  names(k) <- data[,paste0(outcome)]
  return(k)
}

return_format_fstat <-function(data) {
k <- data[, paste0(N, " SNPs (r2 = ", round(rsq*100, digits =2), "%; F-statistics = ",  round(fstat, digits = 0), ")")]
names(k) <- data$hgnc
return(k)
  }

summary_logical <- function(x) {
  k1 <- sum(x)
  k2<- length(x)
  paste0(k1,"/", k2, " (", round(k1/k2, digits = 2)*100, "%)")
}

####Abstract####
data_res[!(grepl("dis-15-", id.outcome)), length(unique(id.outcome))]
data_res[(grepl("dis-15-", id.outcome)), length(unique(id.outcome))]
df_index[grepl("dis-15-", id), max(as.numeric(sample_size))]
k <- data_res[!(grepl("dis-15-", id.outcome)) & pval.outcome < 0.05/1589, ]
data_res[SNP == "rs116843064" & id.outcome %in% c("dis-13-1", "dis-19-1", "dis-23-2"), ] %>% return_format_data
data_res[(grepl("dis-15-", id.outcome)) & pval.outcome < 0.005 & beta.outcome > 0,]
data_res[(grepl("dis-15-", id.outcome)) & pval.outcome < 0.05/1589 & beta.outcome < 0,]
(pearson_cor <- cor(scatter$b_E40K, scatter$b_LPL, use="complete.obs", method = "pearson"))
cor(scatter$b_E40K, scatter$b_LIPC, use="complete.obs", method = "pearson")
####Introduction########
data_res[!(grepl("dis-15-", id.outcome)), id.outcome %>% unique %>% length]

####Results#####
#para2
data_res[SNP == "rs116843064" & id.outcome %in% c("dis-13-1", "dis-19-1", "dis-23-2"), ] %>% return_format_data
data_res[SNP == "rs116843064" & id.outcome %in% c("trait-7-2"), ]
data_res[SNP == "rs116843064" & !grepl("dis-15-", id.outcome) & grepl("dis-", id.outcome) & pval.outcome < 0.05 & beta.outcome > 0]
#para3
cox_data[toinclude==1&E40K_carrier!=TRUE, sum(has_lof)]
cox_data[toinclude==1, sum(has_lof)]
# IV <- "has_lof"
# DV <- c("ldl", "tg", "hdl")
# cov_inc <- paste(c("+ age_enrollment + sex", paste0("PCA", 1:10)), collapse = " + ")
# arguments <- tidyr::expand_grid(IV, DV, cov_inc) %>% as.data.table(.)
# 
# res <- map(split(arguments, 1:arguments[,.N]), function(x) {
#   res <- lm(paste0(x$DV, " ~ ", x$IV, x$cov_inc), data = cox_data[toinclude==1, ]) %>% broom::tidy(.)
#   res$DV <- x$DV
#   return(res)}) %>%
#   rbindlist(., fill = TRUE)
# res[term == "has_lofTRUE", ]
cox_res[exposure=="has_lofTRUE"]
cox_data[CAD_toinclude==1, summary_logical(has_lof)]
cox_res[exposure=="has_lofTRUE" & outcome == "CAD" & include_E40K_ascarrier==FALSE] %>% return_format_data_noexp_HR
cox_res[exposure=="has_lofTRUE" & outcome == "AS" & include_E40K_ascarrier==FALSE] %>% return_format_data_noexp_HR
cox_data[AS_toinclude==1 & AS_censored == 1,has_lof] %>% summary_logical

###para 4
res_cophescan[querytrait%in% c("dis-13-1", "dis-19-1", "dis-23-2", "trait-16-4", "trait-16-1"), ]
res_hyprcoloc$regional_prob
res_hyprcoloc$posterior_explained_by_snp
##### Para 5 #####
data_res[grepl("dis-15-", id.outcome) & SNP == "rs116843064" & beta.outcome > 0 & pval.outcome<0.005, ]
# Association with lymphadenitis?
data_res[grepl("lymphade|ascite|periton", tolower(clean_variable_name)) & pval.outcome < 0.05, ] #lymphadenopathy, ascites, and peritonitis
k<-data_res[grepl("dis-15-", id.outcome) & grepl("steno", tolower(clean_variable_name)) & clean_variable_name == "Calcific aortic valvular stenosis" & pval.outcome < 0.05, ] #lymphadenopathy, ascites, and peritonitis
k%>%return_format_data

####para 6####
inst_tgsd <-  inst[grepl("trait-16-4", id.exposure), ] %>% TwoSampleMR::add_rsq()
setDT(inst_tgsd)
res_rsq_fstat <- map(unique(inst_tgsd$id.exposure), function(x) {
k <- inst_tgsd[id.exposure==x, ]
res<-data.table(id.exposure=x)
res$fstat<-GagnonMR::fstat_fromdat(list(k))
res$rsq <- sum(k$rsq.exposure)
res$N <- k[,.N]
return(res)
}) %>% rbindlist
return_format_fstat(res_rsq_fstat)
cor(scatter$b_E40K, scatter$b_LPL, use="complete.obs", method = "pearson")
cor(scatter$b_E40K, scatter$b_LIPC, use="complete.obs", method = "pearson")
### Discussion ####
#para 3
(pearson_cor <- cor(scatter$b_E40K, scatter$b_LPL, use="complete.obs", method = "pearson"))

# para 5
dt_united <- GagnonMR::get_tpm_for_genes_on_all_tissues(c("ANGPTL4", "LPL", "LIPC"))
dt_united[hgnc=="LPL"][order(-TPM),]
dt_united[tissue=="Adipose Visceral Omentum" & hgnc == "LPL", ]$TPM / dt_united[tissue=="Liver" & hgnc == "LPL", ]$TPM
dt_united[hgnc=="LIPC"][order(-TPM),]
dt_united[tissue=="Liver" & hgnc == "LIPC", ]$TPM / dt_united[tissue=="Adipose Visceral Omentum" & hgnc == "LIPC", ]$TPM


#####Methods######
#para 1
data_res[!(grepl("dis-15-", id.outcome)) & grepl("dis-", id.outcome), length(unique(id.outcome))]
data_res[!(grepl("dis-15-", id.outcome)) & grepl("trait-", id.outcome), length(unique(id.outcome))]

#para 2
data_res[(grepl("dis-15-", id.outcome)), length(unique(id.outcome))]
df_index[grepl("dis-15-", id), max(as.numeric(sample_size))]

#lof mutation para 1
cox_data[!is.na(PCA1), .N]

#lof mutation para 2
cadd_loff[rsid=="rs116843064", cadd_score]
res_logit_tg[term == cadd_loff[rsid=="rs116843064",lof_name], p.value]

#lof mutation para 4
colnom_date <- colnames(cox_data)[grepl("_date$", colnames(cox_data))]
date_end_followup <- cox_data[, lapply(.SD, function(x) max(x, na.rm = TRUE)), .SDcols = colnom_date]
date_end_followup[,apply(.SD, 1, max)]


