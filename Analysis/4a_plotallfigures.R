#!/usr/bin/env Rscript
library(tidyverse)
library(data.table)
library(GagnonMR)
library(gassocplot)
library(forestploter)
library("grid")
library("ggplotify")
library(cowplot)
library(survival)
library(survminer)

wd<-"/mnt/sda/gagelo01/Projects/small_MR_exploration/Test_potential_project/ANGPTL4"
setwd(wd)

all_outcome <- fread("Data/Modified/all_outcome.txt" )
all_outcome <- all_outcome %>% unique
df_index <- fread("/mnt/sda/gagelo01/Vcffile/server_gwas_id.txt")

####align on the good reference allele
k <- all_outcome[id.outcome == "trait-16-4", ] %>% GagnonMR::convert_outcome_to_exposure(.)
harm <- TwoSampleMR::harmonise_data(k, all_outcome)
tosel <- colnames(harm)[!grepl("exposure", colnames(harm))]
data_res <- harm[,tosel] %>% as.data.table(.)
cox_res <- fread("Data/Modified/cox_res.R")
cox_dat <- fread("Data/Modified/cox_data.R")
#####Phewas#####

phewas <- data_res[grepl("dis-15-", id.outcome), ]
info = fread("/mnt/sda/gobemi01/PheWAS/PheWAS_MR_ANGPTL3/data/FinnGenR7_infofile_mod.csv")
phewas  <- merge(phewas , info[,.(phenocode,category)], by.x="outcome",by.y="phenocode")
phewas  <- merge(phewas, df_index[,.(id,clean_variable_name)], by.x = "id.outcome", by.y = "id")
# phewas <- merge(phewas, dt_gene_region[,.(hgnc, UniProt, id, study, gene_region)], by.x = "id.exposure", by.y ="id")


ntest<- phewas[,length(unique(clean_variable_name))]
setnames(phewas, c("outcome", "pval.outcome") ,c("phenotype", "p"))
phewas[,groupnum:=as.factor(category)]
phewas[,direction := ifelse(beta.outcome>=0, "25", "24")]
phewas[,description := clean_variable_name]
phewas[, value:=-log10(as.numeric(p))]


signif_adj = -log10(0.05/ntest)
max_y = max(phewas$value)
if (max_y<=10) { y_interval = 1} else if (max_y>10 & max_y<=20) { y_interval = 2} else { y_interval = 5}
round_max_y = ceiling(max_y/y_interval)*y_interval + 1

phewas[ p<0.05/ntest, annotate := TRUE]
library(MetBrewer)
source("/mnt/sda/gagelo01/Projects/small_MR_exploration/F2_finngen/Analysis/my_phenotypePlot_function.R")

plot <-phenotypePlot(phewas[SNP == "rs116843064",.(groupnum, direction, p, description, phenotype, value, annotate)],
                     suggestive.line=-log10(0.05),significant.line=signif_adj,direction=T,color.palette=rep(met.brewer("Redon"),12),
                     x.group.labels=T,y.axis.interval=y_interval, annotate.size=3, size.x.labels=8, size.y.labels=14,
                     x.axis.label="", y.axis.label="-log(p-value)",
                     title="Phewas of E40K")
plot

ggsave(filename = paste0("Results/angptl4_figure5.png"),
       width = 491/72,height = 486/72,units="in",scale=1, device = "png") 

#####
#######Forest plot with all sorts of outcome
dt<- rbindlist(list(
  data.table(nom = c(paste0("trait-2-", 1:6), paste0("dis-6-", 1:3), "dis-23-2", "trait-29-15", "trait-29-16", "trait-31-2", "trait-31-4"), value = "Glucose homeostasis"),
  data.table(nom = c(paste0("dis-13-", 1:4), "dis-5-1", df_index[grepl("dis-14-", id), id], paste0("trait-13-", 1:2), "dis-18-1", "dis-19-1", "dis-1-1"), value = "Vascular"),
  data.table(nom = c("dis-2-2", "trait-14-8", paste0("trait-27-", 1:2)), value = "Liver"),
  data.table(nom = c("dis-7-1", "trait-12-2"), value = "Kidney"),
  data.table(nom = c(unique(df_index[grepl("trait-16-", id), id]), "trait-19-5", "trait-29-20") , value = "Lipids"),
  data.table(nom =   paste0("trait-28-", 1:2), value = "Bone"),
  data.table(nom = c(unique(df_index[grepl("met-", id), id]), "trait-20-8"), value = "Metabolites"),
  data.table(nom = c("trait-1-1", "trait-10-1", "trait-10-3", df_index[grepl("trait-25-", id), id]), value = "Anthropometric"),
  data.table(nom = c("trait-6-1", "trait-7-2"), value = "Lifespan"),
  data.table(nom = c("dis-12-1", "trait-14-10"), value = "Pancreas"),
  data.table(nom = c("dis-4-1"), value = "Brain")
))

dt_forest <- data_res
dt_forest <- rbind(dt_forest, data.table(id.outcome = c("trait-2-2", "trait-2-4"), SNP = "rs116843064"), fill = TRUE)
dt_forest <- merge(dt_forest, distinct(dt), by.x = "id.outcome", by.y = "nom")
dt_forest <- merge(dt_forest, df_index[,.(id,clean_variable_name)], by.x = "id.outcome", by.y = "id")
dt_forest[id.outcome=="dis-2-2", clean_variable_name:= "NAFLD"]
####
dt_forest[,lci := beta.outcome-1.96*se.outcome]
dt_forest[,uci := beta.outcome+1.96*se.outcome]
setnames(dt_forest, c("beta.outcome", "pval.outcome"), c("b", "pval"))
dt_forest <- merge(data.table(SNP = c("rs115849089", "rs116843064"), hgnc = c("LPL", "ANGPTL4")), dt_forest, by = "SNP")
dt_forest<- dt_forest[clean_variable_name!=""]
dt_forest <- dt_forest[order(value, clean_variable_name),]
dt_forest <- dt_forest[!(id.outcome %in% c("trait-16-36", "trait-16-37", "trait-16-39")),] #remove mixed in GLGC
dt_forest[, .N,by="clean_variable_name"][N==1,]
dt_forest <- dt_forest[!grepl("trait-2-", id.outcome), ]
dt_forest <- dt_forest[hgnc=="ANGPTL4",]
dt_forest[,panel :="ANGPTL4 E40K",]

p <- my_forestplotter_fancy(data = dt_forest[grepl("trait-", id.outcome) & id.outcome != "trait-6-1", ], 
                            col.below.header = "clean_variable_name", 
                            col.header = "value", 
                            col.header.heading = "Outcomes", 
                            effect.name = "beta (95% CI)", 
                            col.right = NULL, 
                            exponentiate = FALSE,
                            xlab = "Standard deviation (95% CI)")

ggsave(p, filename = paste0("Results/","angptl4_figure1", ".png"),
       width = 521/72,height = 546/72,units="in",scale=1, device = "png")

# tosel <- c("dis-23-2", "dis-7-1", "dis-2-2", "dis-13-1", "dis-19-1", "dis-14-6", "dis-1-1",  "dis-18-1")
# dt_forest <- dt_forest[!(id.outcome %in% c("dis-6-1","dis-13-2", "dis-5-1")), ]

p <- my_forestplotter_fancy(data = dt_forest[grepl("dis-", id.outcome)  | id.outcome == "trait-6-1", ], 
                            col.below.header = "clean_variable_name", 
                            col.header = "value", 
                            col.header.heading = "Outcomes", 
                            effect.name = "OR (95% CI)", 
                            col.right = NULL, 
                            exponentiate = TRUE,
                            xlab = "Odds ratio (95% CI)")
ggsave(p, filename = paste0("Results/","angptl4_figure2", ".png"),
       width = 521/72,height = 346/72,units="in",scale=1, device = "png")

#######Create the lof mutation forestplot######
###create the lof mutation with tg, ldl and HDL
#####test each variant individually with TG######
include_E40K_ascarrier <- c(TRUE, FALSE)
IV <- "has_lof"
DV <- c("tg", "ldl","hdl")
cov_inc <- paste(c("+ age_enrollment + sex", paste0("PCA", 1:10)), collapse = " + ")
arguments <- tidyr::expand_grid(IV, DV, cov_inc, include_E40K_ascarrier)
setDT(arguments)
res <- map(split(arguments, 1:arguments[,.N]), function(x) {
  if(x$include_E40K_ascarrier==TRUE){data<-cox_dat[toinclude==1,]}else{
    data <- cox_dat[toinclude==1 & E40K_carrier == FALSE,]}
  res<- lm(paste0(x$DV, " ~ ", x$IV, x$cov_inc), data = data) %>%
    broom::tidy(.) %>% as.data.table(.)
  res[,DV:=x$DV]
  res[,include_E40K_ascarrier:=x$include_E40K_ascarrier]
  return(res)}) %>%
  rbindlist(., fill = TRUE)

dt_forest <- res[term=="has_lofTRUE",]
setnames(dt_forest, c("estimate", "std.error", "p.value"), c("b", "se", "pval"))
dt_forest[,lci:=b-se*1.96]
dt_forest[,uci:=b+se*1.96]
dictionnary <- c(tg = "Triglycerides",  ldl = "LDL cholesterol", hdl = "HDL cholesterol")
dt_forest[,panel := "ANGPTL4 loss-of-function \nmutations carriers  "]
dt_forest[,outcome := dictionnary[DV]]
dt_forest <- dt_forest[order(outcome)]
dt_forest[,header:=ifelse(include_E40K_ascarrier, "Including E40K carriers", "Excluding E40K carriers")]

data<- data.table(dt_forest)
p1 <- my_forestplotter_fancy(data = data, 
                            col.below.header = "outcome", 
                            col.header = "header", 
                            col.header.heading = "Outcomes", 
                            col.left = NULL,
                            col.left.heading = NULL,
                            effect.name = "Beta (95% CI)", 
                            col.right = NULL, 
                            exponentiate = FALSE,
                            xlab = "1-SD deviation scale")
p1  

# ggsave(plot = p1, filename = paste0("Results/Presentation/","forest_lof_lipids", ".png"),
#        width = 450/72,height = 175/72,units="in",scale=1, device = "png")
####create the table prevalence#####
colnom <- apply(expand.grid(c("CAD", "AS", "T2D"), c("_toinclude", "_censored")), 1, paste, collapse="")

k<- cox_dat[, .SD,.SDcols = c("f.eid", "has_lof", "E40K_carrier", colnom) ]
k_long = melt(k, id.vars = c("f.eid", "has_lof", "E40K_carrier"))
k_long <- separate(k_long, col = "variable", into = c("outcome", "variable"))
k_wide <- dcast(k_long, f.eid + has_lof + outcome + E40K_carrier ~ variable, value.var = "value")

summary_logical <- function(x) {
  k1 <- sum(x)
  k2<- length(x)
  paste0(k1,"/", k2, " (", round(k1/k2, digits = 3)*100, "%)")
}

getdfnew<-function(k_wide){
  df_new <- k_wide[, summary_logical(.SD[toinclude==1, censored]), by = c("outcome", "toinclude", "has_lof")][toinclude==1,]
  df_new[, has_lof := ifelse(has_lof, "carriers", "noncarriers")]
  df_new <- dcast(df_new, outcome ~ has_lof, value.var = "V1")
  return(df_new)}
k1<-getdfnew(k_wide)
k1[ ,include_E40K_ascarrier:=TRUE]
k2<-getdfnew(k_wide[E40K_carrier!=TRUE])
k2[ ,include_E40K_ascarrier:=FALSE]
df_new <- rbind(k1,k2)
#####create the dt_forest#######
dt_forest <- merge(cox_res[exposure=="has_lofTRUE",], df_new, by = c("outcome","include_E40K_ascarrier"))
dt_forest[,b:=log(HR)]
dt_forest[,z_score := qnorm(1 - pval/2)*(b/abs(b))]
dt_forest[,se:=b/z_score]
dt_forest[,lci:=b-se*1.96]
dt_forest[,uci:=b+se*1.96]
# dt_forest[,panel:=stringr::str_pad("", 45, "right", " ")]
dt_forest[,header:=ifelse(include_E40K_ascarrier, "Including E40K carriers", "Excluding E40K carriers")]
dt_forest[,panel := paste0("ANGPTL4 loss-of-function \nmutation carriers")]
k<- c(AS = "Aortic stenosis", CAD = "Coronary artery disease", T2D = "Type 2 diabetes")
dt_forest[,outcome:=k[outcome]]
data<-data.table(dt_forest)
p2 <- my_forestplotter_fancy(data = data, 
                            col.below.header = "outcome", 
                            col.header = "header", 
                            col.header.heading = "Outcomes", 
                            col.left = c("carriers", "noncarriers"),
                            col.left.heading = c("Prevalence in \ncarriers", "Prevalence in \nnoncarriers"),
                            effect.name = "HR (95% CI)", 
                            col.right = NULL, 
                            exponentiate = TRUE,
                            xlab = "")

# ggsave(plot = p2, filename = paste0("Results/figure3", ".png"),
#        width = 700/72,height = 115/72,units="in",scale=1, device = "png")

####plot lof pane#####
p1gg<-as.ggplot(p1)
p2gg<-as.ggplot(p2)
# panel_3 <- plot_grid(p1gg, p2gg, nrow = 2, labels=c("A)", "B)"), align = "v")
# class(panel_3)
# ggsave(plot = panel_3, filename = paste0("Results/angptl4_figure3", ".png"),
#        width = 727/72,height = 435/72,units="in",scale=1, device = "png")

######create the coxplot########
cox_dat <- fread("Data/Modified/cox_data.R")
cox_plot <- data.table(cox_dat)
cox_plot[,toinclude:=NULL]
colnom1 <- c("f.eid", "E40K_carrier", "has_lof")
colnom2 <- colnames(cox_plot)[grep("toinclude$|censored$|time$", colnames(cox_plot))]
cox_plot <- cox_plot[, .SD, .SDcols =c(colnom1,colnom2)]
tobind <- expand_grid(f.eid = cox_plot$f.eid, include_E40K = c(TRUE, FALSE)) %>% as.data.table(.)
cox_plot <- merge(tobind, cox_plot, by = "f.eid", all.x = TRUE)
cox_plot <- melt(cox_plot, id.vars = c(colnom1, "include_E40K") , measure.vars = colnom2)
cox_plot[, c("disease", "forwide") := data.table::tstrsplit(variable, "_", fixed=TRUE)]
cox_plot[,variable:=NULL]
cox_plot <- dcast(cox_plot, 
                  formula = paste0(paste(c(colnom1, "include_E40K", "disease"),  collapse = " + "), "~ forwide"),
                  value.var = "value")
cox_plot[,time:=time/365.25]
cox_plot[include_E40K==FALSE & E40K_carrier==TRUE, toinclude := 0]

fit <- survfit( Surv(time, censored) ~ has_lof, data = cox_plot[toinclude==1,] )
cox_figure <- ggsurvplot_facet(fit, cox_plot[toinclude==1,], facet.by = c("disease", "include_E40K"),
                 palette = "jco", pval = TRUE, linetype = "solid", censor = FALSE,
                 pval.coord = c(0,0.93)) +
  xlab("Age (years)") +
  ylab("Disease-free probability") +
  ylim(c(0.7,1)) +
  guides(color=guide_legend(title="ANGPTL4 mutation carriers"))

panel_3 <- cowplot::plot_grid(p2gg, p1gg, labels = c('A)', 'B)'),
                              label_size = 12, ncol = 1)

ggsave(plot = panel_3, filename = paste0("Results/angptl4_figure3", ".png"),
       width = 727/72,height = 450/72,units="in",scale=1, device = "png")

# panel_3 <- cowplot::plot_grid(cox_figure, p2gg, p1gg, labels = c('A)', 'B)', "C)"),
#                               label_size = 12, ncol = 1, rel_heights = c(2, 1, 1))
# 
# ggsave(plot = panel_3, filename = paste0("Results/angptl4_figure3", ".png"),
#        width = 727/72,height = 900/72,units="in",scale=1, device = "png")
#####hyprcoloc plot########
inst_mvmr <- fread("Data/Modified/data_hyprcoloc.txt")
res <- fread("Data/Modified/res_hyprcoloc.txt")
order_exposure <-  c("HDL cholesterol", "Triglyceride", "Coronary artery disease","Type 2 diabetes","Aortic stenosis" )
source("Analysis/tosource.R")

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


png("Results/angptl4_figure4.png",  width=650/72,height=850/72,  unit = "in", res =300)
twopanel
dev.off()

#####Concordance with LPL enhancement######
res_tg_sd <- fread("Data/Modified/res_tg_sd.txt")
res_tg_sd[, hgnc:= "E40K"]
res_multicis <- fread( "Data/Modified/res_multicis_independent.txt")
dt_pan <- res_multicis[grepl("trait-16-4", id.exposure) & method == "Inverse variance weighted", ]
k1<- separate(dt_pan[grepl("trait-16-4", id.exposure), ], col = "id.exposure", into = c("hgnc", "id.exposure"), sep = "_", remove = TRUE)
k1 <- rbind(k1, res_tg_sd, fill = TRUE)
scatter <- dcast(k1[method%in%c("Inverse variance weighted",  "Wald ratio") & !grepl("dis-15-", id.outcome),], id.outcome ~ hgnc, value.var = c("b", "se"))
colnom<-colnames(scatter)[grepl("^b_", colnames(scatter))]
cor(scatter[,.SD,.SDcols = colnom], method = c("pearson"), use="complete.obs")
scatter <- merge(scatter, distinct(dt), by.x = "id.outcome", by.y = "nom")
data.scatter <- data.table(scatter)

plot_scatter <- function(scatter,
                         my_xlab = "Genetically proxied ANGPTL4 inhibition",
                         my_ylab = "Genetically proxied LPL enhancing",
                         title) {
  pearson_cor <- cor(scatter$b_x, scatter$b_y, use="complete.obs", method = "pearson")
  lm_res <- lm(formula = "b_x ~ b_y", data = scatter) %>%
    broom::tidy() %>%
    as.data.table()
  
  
  ggplot2::ggplot(data = scatter, ggplot2::aes(x = b_x,
                                               y = b_y, colour = value)) + ggplot2::geom_errorbarh(ggplot2::aes(xmin = b_x -
                                                                                                                  se_x, xmax = b_x + se_x),
                                                                                                   colour = "grey", height = 0) + ggplot2::geom_errorbar(ggplot2::aes(ymin = b_y -
                                                                                                                                                                        se_y, ymax = b_y + se_y),
                                                                                                                                                         colour = "grey", width = 0) + ggplot2::geom_point(size =2) +
    ggplot2::geom_abline(ggplot2::aes(intercept = lm_res[1,estimate], slope = lm_res[2, estimate]), linetype = "dashed") +
    
    ggplot2::scale_colour_manual(values = c(rgb(112, 54, 153, maxColorValue = 255),
                                            rgb(66, 94, 191, maxColorValue = 255),
                                            rgb(84,201,237, maxColorValue = 255),
                                            rgb(59,166,18, maxColorValue = 255),
                                            rgb(255,110,26, maxColorValue = 255),
                                            rgb(149,199,71, maxColorValue = 255),
                                            rgb(161,15,125, maxColorValue = 255),
                                            # rgb(249, 106, 27, maxColorValue = 255),
                                            # rgb(214,15,102, maxColorValue = 255),
                                            rgb(8, 161, 217, maxColorValue = 255),
                                            rgb(255,186,51, maxColorValue = 255),
                                            rgb(54, 150, 214, maxColorValue = 255))) +
    ggplot2::labs( x = my_xlab, y = my_ylab ) +
    ggplot2::theme(legend.position = "top",legend.direction = "vertical") + ggplot2::guides(colour = ggplot2::guide_legend(ncol = 2)) +
    expand_limits(x = 0, y = 0) +
    theme(
      legend.background = element_rect(fill = "white", size = 4, colour = "white"),
      axis.ticks = element_line(colour = "black", size = 0.4),
      axis.title=element_text(size=14,face="bold"),
      # panel.grid.major = element_line(colour = "grey70", size = 0.2),
      panel.grid.minor = element_blank(),
      legend.key=element_blank(),
      panel.background = element_rect(fill = 'white', colour = 'white'),
      # panel.border = element_rect(colour = "grey70", fill=NA, size=1),
      axis.line.x.bottom = element_line(color = "black", size = 0.2),
      axis.line.y.left   = element_line(color = "black", size = 0.2),
      axis.line.y.right  = element_blank(),
      axis.text.y.right  = element_blank(),
      panel.border       = element_blank(),
      legend.title= element_blank(),
      legend.margin=margin(c(1,5,5,5)),
      plot.caption = element_text(hjust=0.5),
      plot.title = element_text(hjust = 0.5, size=14,face="bold"),
      legend.position=c(0.3,0.85)) +
    labs(caption = "[Each measures are scaled to 1-SD lower triglyceride levels]", size = 4) +
    annotate(geom="text", x=0.8, y=-0.8, label= paste0("~italic(r)==", round(pearson_cor, digits = 2)), parse=TRUE, size=7) +
    ggtitle(title)
}

scatter <- data.table(data.scatter)
scatter <- scatter[,.(b_E40K, se_E40K, b_LPL, se_LPL, value)]
setnames(scatter, c("b_E40K", "se_E40K", "b_LPL", "se_LPL"), c("b_x", "se_x", "b_y", "se_y"))
fig6_a <- plot_scatter(scatter = scatter,
                       my_xlab = "Genetically proxied ANGPTL4 inhibition",
                       my_ylab = "Genetically proxied LPL enhancing",
                       title = "ANGPTL4 inhibition vs. LPL enhancement")

scatter <- data.table(data.scatter)
scatter <- scatter[,.(b_E40K, se_E40K, b_LIPC, se_LIPC, value)]
setnames(scatter, c("b_E40K", "se_E40K", "b_LIPC", "se_LIPC"), c("b_x", "se_x", "b_y", "se_y"))
fig6_b<- plot_scatter(scatter = scatter,
                      my_xlab = "Genetically proxied ANGPTL4 inhibition",
                      my_ylab = "Genetically proxied HL enhancing",
                      title = "ANGPTL4 inhibition vs. HL enhancement")

panel_6 <- cowplot::plot_grid(fig6_a, fig6_b, labels = c('A)', 'B)'),
                              label_size = 22, ncol = 2)
ggsave(plot = panel_6, filename = paste0("Results/angptl4_figure6", ".png"),
       width = 1000/72,height = 500/72,units="in",scale=1, device = "png")

message("This script finished without errors")