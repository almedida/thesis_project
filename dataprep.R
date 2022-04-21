rm(list = ls())
##############################################
#install.packages('pacman')
library(pacman, devtools)
p_load("tidyverse", "matrixTests", "gtools", "tmvtnorm")

#install BiocManager
#if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")

#install limma and qvalue packages
#BiocManager::install(c("limma", "qvalue"))
library(limma)
library(qvalue)
#############################################
kidney_CBA = read.table("kidneyCBA.csv", header= TRUE, sep=",")
kidney_Tech = read.table("kidneyTech.csv", header= TRUE, sep=",")

#load liver data
liver_CBA = read.table("liverCBA.csv", header= TRUE, sep=",")
liver_Tech = read.table("liverTech.csv", header= TRUE, sep=",")

###################student t test ######################################3

kidney_ttest = row_t_equalvar(kidney_CBA, kidney_Tech, alternative = "two.sided", 
                              mu = 0, conf.level = 0.95)
liver_ttest =  row_t_equalvar(liver_CBA, liver_Tech, alternative = "two.sided", 
                              mu = 0, conf.level = 0.95)

pval_raw = list(kidney_ttest$pvalue, liver_ttest$pvalue) %>% as.data.frame()
colnames(pval_raw) = c("kidney_pval", "liver_pval")

pvals1 = pval_raw["kidney_pval"] %>% as.matrix()
pvals2 = pval_raw["liver_pval"] %>% as.matrix()

adj_pval_kidney = p.adjust(data.matrix(pvals1, "BH")) %>% as.data.frame()
adj_pval_liver = p.adjust(data.matrix(pvals2, "BH")) %>% as.data.frame()
adj_pvals = list(adj_pval_kidney, adj_pval_liver) %>% as.data.frame()
###############################################################
#######################limma t test#######################################
#gene expression analysis for the 2 levels (Tech and CBA) for kidney
#load table on data frame kidneyExpr
kidney_expr = list(kidney_CBA, kidney_Tech) %>% as.data.frame()

#create a factor list for the differentially expressed genes with Tech set as first level
kidney = factor(
  x = c(rep("Tech",5), rep("CBA",5)),
  levels=c("Tech","CBA")            # Set Tech to be the first level
)

design_kidney = model.matrix(~kidney)          # Remove the zero

#Now we can run the differential expression pipeline
fit_kidney = lmFit(kidney_expr, design_kidney)
fit_kidney = eBayes(fit_kidney)
results_kidney = decideTests(fit_kidney)

kidney_limma_pval = fit_kidney$p.value[, 2] %>% as.matrix()
adj_kidney_limma_pval = p.adjust(data.matrix(kidney_limma_pval, "BH")) %>% as.data.frame()

#visualization of the results
#plotMD(fit_kidney, coef="kidneyCBA", status=results_kidney)

#####################################################################
#load table on data frame liverExpr
liver_expr = list(liver_CBA, liver_Tech) %>% as.data.frame()

#create a factor list for the differentially expressed genes with Tech set as first level
liver <- factor(
  x = c(rep("Tech",6), rep("CBA",6)),
  levels=c("Tech","CBA")            # Set Tech to be the first level
)

design <- model.matrix(~liver)          # Remove the zero

#Now we can run the differential expression pipeline
fit_liver <- lmFit(liver_expr, design)
fit_liver <- eBayes(fit_liver)
results_liver <- decideTests(fit_liver)

liver_limma_pval = fit_liver$p.value[, 2] %>% as.matrix()
adj_liver_limma_pval = p.adjust(data.matrix(liver_limma_pval, "BH")) %>% as.data.frame()

#visualization of the results
#plotMD(fit_liver, coef="liverCBA", status=results)_liver

adj_limmas = list(adj_kidney_limma_pval, adj_liver_limma_pval) %>% as.data.frame()
