#!/usr/bin/env Rscript
library(cophescan)
library(data.table)
library(tidyverse)
library(GagnonMR)
library(gassocplot)
wd<-"/mnt/sda/gagelo01/Projects/small_MR_exploration/Test_potential_project/ANGPTL4"
setwd(wd)
source("Analysis/tosource.R")
generegion  = GagnonMR::from_genecard_to_generegion("ANGPTL4", window = 6e5)
gwasvcf::set_bcftools()
gwasvcf::set_plink()
res_cophescan <- fread( "Data/Modified/res_cophescan_hyperpriors.txt")
res_cophescan[cophe.hyp.call=="Hc" & !grepl("dis-15", querytrait)]
df_index <- fread("/mnt/sda/gagelo01/Vcffile/server_gwas_id.txt")
#perform hyprcoloc analyses
ID <- c("trait-16-4", "trait-16-1","dis-13-1", "dis-23-2", "dis-19-1")
inst <- GagnonMR:::intern_vcfpath_to_TwoSampelMR_region(vcf = paste0("/mnt/sda/gagelo01/Vcffile/Server_vcf/", ID,"/", ID, ".vcf.gz"), chrompos = GagnonMR::from_genecard_to_generegion("ANGPTL4", window = 6e5), parameters = parameters)
inst <- merge(inst, df_index[, .(id, clean_variable_name)], by.x = "id.exposure", by.y = "id")
inst[,exposure:=clean_variable_name]
order_exposure <-  c("HDL cholesterol", "Triglyceride", "Coronary artery disease","Type 2 diabetes","Aortic stenosis" )
inst_mvmr <- GagnonMR::prepare_for_mvmr(inst, inst, should_clump = FALSE, harmonise_strictness = 1)
res <- GagnonMR::run_hypr_on_aligned(inst_mvmr)
inst_mvmr[, chr.exposure := as.integer(chr.exposure)]

fwrite(inst_mvmr, "Data/Modified/data_hyprcoloc.txt")
fwrite(res, "Data/Modified/res_hyprcoloc.txt")
A <- GagnonMR::stack_assoc_plot_wrapper(df_aligned = inst_mvmr,res_hypr1 = res, ldref = parameters$ldref,
                                   traits_inorder = order_exposure, build = 37)

res <- sensitivity.plot_wrapper(df_aligned = inst_mvmr,
                                traits_inorder = order_exposure)
B<-drawheatmap(res[[2]])
twopanel <-  cowplot::ggdraw() +
  cowplot::draw_plot(ggplotify::as.ggplot(A) + theme(text = element_text(size = 0.4)), x = 0.08, y =0, width = .6, height = 1) +
  cowplot::draw_plot(B, x = .66, y =0.1, width = .31, height = 0.7) +
  cowplot::draw_plot_label(label = c("", ""), size = 25,
                  x = c(0, 0.62), y = c(0.9, 0.9))


png("Results/figure4.png",  width=650/72,height=850/72,  unit = "in", res =300)
twopanel
dev.off()

# ggsave(plot = twopanel, filename = paste0("Results/figure4.png"),
#        width = 790/72,height = 683/72,units="in",scale=1, device = "png")

message("This script finished without errors")
