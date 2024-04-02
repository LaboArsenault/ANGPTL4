#!/usr/bin/env Rscript
library(data.table)
library(tidyverse)
library(GagnonMR)
library(tictoc)
library(furrr)

wd<-"/mnt/sda/gagelo01/Projects/small_MR_exploration/Test_potential_project/ANGPTL4"
setwd(wd)
gwasvcf::set_bcftools()
gwasvcf::set_plink()
ldref<-"/home/couchr02/Mendel_Commun/Christian/LDlocal/big_EUR_rs"
x <- paste0("awk -F ' ' '{print $2}' ","/home/couchr02/Mendel_Commun/Christian/LDlocal/EUR_rs", ".bim")
dt_ref_snp <- fread(cmd = x, header = FALSE) #Those are SNP in 1000G with MAF>0.01
df_index <- data.table::fread("/mnt/sda/gagelo01/Vcffile/server_gwas_id.txt")
if(file.exists("Data/Modified/dt_gene_region.txt")){dt_gene_region<-fread("Data/Modified/dt_gene_region.txt")} else{
dt_gene_region <- fread("/mnt/sda/gagelo01/Projects/small_MR_exploration/Test_potential_project/Tribal/Data/Modified/dt_gene_region.txt")}
ncores <- 40
options(future.globals.maxSize= 5e9)
plan(multicore, workers = ncores, gc = TRUE)

######change parameters#########
######
parameters <- GagnonMR::default_param()
parameters$path <- c(GagnonMR::default_param()$path, paste0(wd, "/Data/Modified/Gtex_vcf/"))
parameters$uni_cis_minpval <- 1
parameters$ldref <- ldref
parameters$snp_bim <- dt_ref_snp$V1 #I only include SNP in 1000G with maf > 0.01
parameters$multicis_clumping["clump_r2"]<-0.6


#####select inst#####
holistic_selection<-FALSE
study_to_selectinst <- list(c(NULL))
study_noinst_butexposure <- NULL
typeof_sel = c("lead_snp", "multicis_independent_clumping")
should_skip_homogenous = TRUE
should_select_pan<-FALSE
should_QTL_instrument_selection <- FALSE
############Choose the gene you wish to include and the outcome you wish to include  ###########
gene_toinclude <-  c("ANGPTL4", "LIPC", "LPL")
# vec_tissue_gtex <- c("Liver", "Adipose_Subcutaneous", "Adipose_Visceral_Omentum")
thewindow <- 3e5
arguments_exposure_proxies <- dt_gene_region[hgnc%in%gene_toinclude,]
arguments_exposure_proxies

arguments_exposure_proxies[,gene_region:=GagnonMR::from_genecard_to_generegion(hgnc, 6e5)[hgnc]]
arguments_exposure_proxies[,vcffile_inst := "/mnt/sda/gagelo01/Vcffile/Server_vcf/trait-16-4/trait-16-4.vcf.gz"]
############Choose the outcome you wish to include  ###########
# ID_server_out <- c("dis-14-7", "trait-31-4", "trait-31-2")
ID_server_out <- c("dis-13-1",
                   # df_index[grepl("dis-15-",id) & ncase > 1000, id],
                   df_index[grepl("trait-16", id) & grepl("^HDL|^LDL|^logTG", trait) & population %in% c("European") & sex == "Males and Females", id],
                   "trait-19-5",
                   "trait-7-2", "trait-6-1","dis-14-6", "dis-18-1", "dis-19-1",
                   "dis-2-1", "trait-14-8", "trait-14-10", paste0("trait-27-", 1:2), "dis-23-2",
                   "trait-36-2", "dis-7-1",  "trait-29-20", "dis-12-1",
                   "trait-13-1", "trait-13-2", "trait-1-1", "trait-10-1",
                   "dis-14-7", "trait-29-15", "trait-29-16", "trait-31-4", "trait-31-2") %>% unique#"trait-29-7"
## "trait-25-4", "trait-25-13", "trait-25-19", "trait-25-16", "trait-25-1", "dis-20-1",  "dis-4-1", "trait-2-2", "trait-2-4",
out_server <- paste0("/mnt/sda/gagelo01/Vcffile/Server_vcf/", ID_server_out, "/", ID_server_out, ".vcf.gz")
out_mrbase<-NULL

#####methods
split_outcome<-FALSE
all_mr_methods_short = TRUE
all_mr_methods_skip_presso = TRUE
tsmr_robust<-c("mr_weighted_mode", "mr_weighted_median")

make.dir <- function(fp) {
  if(!file.exists(fp)) {  # If the folder does not exist, create a new one
    dir.create(fp)
  } else {   # If it existed, delete and replace with a new one  
    unlink(fp, recursive = TRUE)
    dir.create(fp)
  }
}
