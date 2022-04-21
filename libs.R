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
