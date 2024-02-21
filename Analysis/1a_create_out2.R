#!/usr/bin/env Rscript
setwd("/home/gagelo01/workspace/Projects/small_MR_exploration/Test_potential_project/ANGPTL4")
source("Analysis/tosource.R")

#####create out#####
out_server <- c(df_index[grepl("dis-15-", id)&ncase>1000, paste0("/mnt/sda/gagelo01/Vcffile/Server_vcf/", id, "/", id, ".vcf.gz")], out_server)
all_outcome <- fread("Data/Modified/all_outcome.txt")

out2 <- future_map(out_server, function(x, rsiid = "rs116843064") {
  gwasvcf::set_bcftools()
  gwasvcf::set_plink()
  outr<-GagnonMR::extract_outcome_variant(snps = rsiid, outcomes = x, rsq = 0.8, parameters = parameters)
  return(outr)
}, .options = furrr_options(seed = TRUE)) %>% rbindlist(.,fill = TRUE)


out<-rbind(all_outcome, out2, fill = TRUE) %>% distinct(.)
out<- distinct(out)
fwrite(out, "Data/Modified/all_outcome.txt" )

message("This script finished without errors")
