exp1 = read.csv("../simulation_4/E_trep1_10k.csv")
exp2 = read.csv("../simulation_4/E_trep2_10k.csv")

#exp1 = read.csv("E_trep1_10k.csv")
#exp2 = read.csv("E_trep2_10k.csv")

# AN  = read.csv("Adjacent_non_tumor.csv", header = TRUE, sep = ",")
# LA = read.csv("Lung_adenocarcinoma.csv", header = TRUE, sep = ",")
# 
# #paired t test
# exp_paired = row_t_paired(AN, LA, alternative = "two.sided", mu = 0, conf.level = 0.95)
# head(exp_paired)
# exp_mean_diff = exp_paired$mean.diff
# # 
# mod_LA = LA + exp_mean_diff
# #mod_exp_paired = row_t_paired(AN, mod_exp2, alternative = "two.sided", mu = 0, conf.level = 0.95)
# 
# hist(mod_exp_paired$pvalue)
# m = 10000
# nn_genes = nrow(AN)
# selected_g = sample(nn_genes, m)
# 
# exp1 = AN[selected_g, ]
# exp2 = mod_LA[selected_g, ]
# 
# write.csv(exp1, "exp1_tumor_10k_new.csv", row.names = FALSE)
# write.csv(exp2, "exp2_tumor_10k_new.csv", row.names = FALSE)

#exp1 = read.csv("exp1_tumor_10k_new.csv", header = TRUE, sep=",")
#exp2 = read.csv("exp2_tumor_10k_new.csv", header = TRUE, sep=",")
###############################################################
#randomly selected genes of size m from the total genes
# #both experiments have the same genes
set.seed(1)
m= 10000
n_genes = nrow(exp1) #n_gene = no of genes (rows in the dataset)
selected_genes = sample(n_genes, m) #the selected_genes will be the same for both experiments
set.seed(NULL)
#set.seed(Sys.time())




