#!/usr/bin/env Rscript
library(cophescan)
library(data.table)
library(tidyverse)
library(GagnonMR)
wd<-"/mnt/sda/gagelo01/Projects/small_MR_exploration/Test_potential_project/ANGPTL4"
setwd(wd)
source("Analysis/tosource.R")
generegion  = GagnonMR::from_genecard_to_generegion("ANGPTL4", window = 6e5)
gwasvcf::set_bcftools()
gwasvcf::set_plink()
set.seed(1)
########
###
from_tsmr_to_cophescan_input <- function(dat) {
  type <- ifelse(all(is.na(dat$ncase.exposure) & 
                       is.na(dat$ncontrol.exposure)), "quant", "cc")
  out <- dat %>% {
    list(pvalues = .$pval.exposure, N = .$samplesize.exposure,
         MAF = .$eaf.exposure %>% ifelse(.<0.5, ., 1-.),
         beta = .$beta.exposure,
         varbeta = .$se.exposure^2, type = type, snp = .$SNP,
         z = .$beta.exposure/.$se.exposure, chr = .$chr.exposure,
         pos = .$pos.exposure, id = .$id.exposure)}
  if (type == "cc") {
    out$s <- mean(dat$ncase.exposure/dat$samplesize.exposure, na.rm = TRUE)
  }
  return(out)
}
#########
###We identify variants###
parameters$ldref<-"/home/couchr02/Mendel_Commun/Christian/LDlocal/EUR_rs" #This is better because LD will better match common variants. These are also variants that are generally well imputed.
ldref<-parameters$ldref
x <- paste0("awk -F ' ' '{print $1, $2, $4}' ",ldref, ".bim")
dt_ref_snp <- fread(cmd = x, header = FALSE) 
colnames(dt_ref_snp)<-c("chr", "rsid", "pos")
k <- strsplit(generegion, split = ":|-")[[1]]
variants <- dt_ref_snp[chr==as.integer(k[1]) & pos > as.numeric(k[2]) & pos<as.numeric(k[3]),]$rsid

#####we create the LD matrix
LD_matrix<- ieugwasr::ld_matrix_local(variants = variants, 
                                      plink_bin = genetics.binaRies::get_plink_binary(), bfile = parameters$ldref, with_alleles = TRUE)

###We extract the corresponding inst
plan(multicore, workers = 40, gc = TRUE)

res_cophescan <- future_map(unique(out_server), function(x) { 
message(paste0("****Cophescan on ", x, "****"))
dat <- GagnonMR:::intern_vcfpath_to_TwoSampelMR_region(vcf=x, chrompos = generegion) 

###We align both 
k <- TwoSampleMR::harmonise_ld_dat(x =dat, ld = LD_matrix)

cophescan_input <- from_tsmr_to_cophescan_input(k$x)
cophescan_input$N<-mean(cophescan_input$N) ##This line is necessary otherwise there will be complex bug further down the pipeline
cophescan_input$LD <- k$ld
if(!("rs116843064" %in% cophescan_input$snp)){return(data.table())} 
res.susie <- cophescan::cophe.susie(cophescan_input, querysnpid = "rs116843064", querytrait=dat$id.exposure[1])

return(summary(res.susie))
}, .options = furrr_options(seed = TRUE)) %>% rbindlist(., fill = TRUE)


#cophesscan susie with hierachical priors
covar=FALSE
covar_vec = NULL
cophe.hier.res <- run_metrop_priors(res_cophescan, avg_posterior=TRUE, avg_pik = TRUE, covar_vec = covar_vec, covar = covar, nits = 30000)
res.post.prob = cbind(cophe.hier.res$avg.posterior, cophe.hier.res$data)

res.hier.predict <- cophe.hyp.predict(as.data.frame(res.post.prob ))
col_disp <- c( "PP.Hn", "PP.Ha", "PP.Hc", "nsnps", "querysnp", "querytrait", "typeBF",  "grp", "cophe.hyp.call")
res_cophescan_hyperpriors <- res.hier.predict[, col_disp]

fwrite(res_cophescan_hyperpriors, "Data/Modified/res_cophescan_hyperpriors.txt")

message("This script finished without errors")
