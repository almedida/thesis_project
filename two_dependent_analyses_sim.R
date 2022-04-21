source("libs.R")

##############################################
genej_sim_dataj = read.csv("genej_sim_data2.csv", header = TRUE, sep = ",")

############################################
m=10000
simulation = function(n_samples, n_mu, n_de){
  # n_de = 1000
  # m=10000
  # n_mu = 2
  # n_samples = 4
  N = ncol(genej_sim_dataj) #N = number columns (samples in the dataset)
  n_gene = nrow(genej_sim_dataj) #n_gene = no of genes (rows in the dataset)
  
  genej_sim_data = genej_sim_dataj
  genej_std = apply(genej_sim_data, 1, sd)  #calculate std for all genes
  
  randomly_selected_samples  = sample(N, 2*n_samples) #randomly select samples of size 2*n_trt_grp 
  #for treatment group
  treatment_group1 = randomly_selected_samples[1:n_samples]
  treatment_group2 = randomly_selected_samples[(n_samples + 1):(2*n_samples)]
  
  treatment_group1_data = genej_sim_data[treatment_group1]
  
  #group 2 data
  treatment_group2_data = genej_sim_data[treatment_group2]
  
  #generate treatment effects, j_effects=10000 for the m genes
  mu = n_mu*genej_std
  treatment_effect_j = rnorm(m, mu, genej_std)
  
  #add 10k treatments effects generated to group 2 data
  de_treatment_group2 = treatment_group2_data + treatment_effect_j
  
  n_treated_genes_grp2 = nrow(de_treatment_group2)
  
  #we extract m genes from DE_treatment_group2 (genes with treatment effects)
  de_genes_index = sample(n_treated_genes_grp2, n_de)
  
  #we copied untreated (genes without effects) genes from treatment_group2_data
  treated_genes_group2 = treatment_group2_data
  
  #we replaced the m untreated genes (without treatment effects) with the corresponding treated genes(with treatment effects)
  treated_genes_group2[de_genes_index, ] = de_treatment_group2[de_genes_index, ]
  
  group1_data = treatment_group1_data
  group2_data = treated_genes_group2
  
  #student t-test
  ttest = row_t_equalvar(group1_data, group2_data, alternative = "two.sided", mu = 0, conf.level = 0.95)
  
  pvalue = ttest$pvalue %>% as.data.frame()
  adj_pvalue = p.adjust(data.matrix(pvalue), "BH") %>% as.data.frame()
  
  #limma 
  limma_gene_exp = c(group1_data, group2_data) %>% as.data.frame()
  
  expr_list_10k = factor(
    x = c(rep("group1", n_samples), rep("group2", n_samples)),
    levels=c("group1","group2")            # Set group 1 to be the first level
  )
  
  design = model.matrix(~expr_list_10k)          # Remove the zero
  
  fit = lmFit(limma_gene_exp, design)
  fit = eBayes(fit)
  results = decideTests(fit)
  
  limma_pvalue = fit$p.value[, 2] %>% as.data.frame()
  adj_limma_pvalue = p.adjust(data.matrix(limma_pvalue, "BH")) %>% as.data.frame()
 
  pval_raw = c(pvalue, limma_pvalue) %>% as.data.frame()
  pvals1 = pval_raw[,1] %>% as.matrix()
  pvals2 = pval_raw[,2] %>% as.matrix()
  
  ##############################################################
  
  calc.cutoff = function(p, B = 20, max=1){
    
    m <- length(p)
    m0 <- m
    bin <- c(-0.1, (1:B)/B*max)
    bin.counts=rep(0,B)
    
    for(i in 1:B){
      bin.counts[i]=sum((p>bin[i])&(p<=bin[i+1]))
    }
    
    tail.means <- rev(cumsum(rev(bin.counts))/(1:B))
    temp <- bin.counts - tail.means
    index <- min((1:B)[temp <= 0])
    cutoff2 <- (index)/B*max
    if(cutoff2 == 1) {cutoff2 <- 1-1/B}
    
    return(cutoff2)
    
  }
  
  cutoff_value1 = calc.cutoff(pvals1, B=20, max=1)
  cutoff_value2 = calc.cutoff(pvals2, B=20, max=1)
  
  cutoff = cbind(c(cutoff_value1), c(cutoff_value2))
  
  colnames(cutoff) = c("cutoff_value1", "cutoff_value2")
  
  ###################################################################
  
  p_vals = pval_raw  %>% filter(pvalue >=cutoff_value1, limma_pvalue>=cutoff_value2)
  
  #################################################################
  
  z_val = as.data.frame(qnorm(as.matrix(p_vals), lower.tail = TRUE))
  colnames(z_val) = c("zvals1", "zvals2")
  
  zvals1 = z_val[,1] %>% as.data.frame()
  zvals2 = z_val[,2] %>% as.data.frame()
  
  ###############################################################
  
  z_val_extremums = as.data.frame(qnorm(as.matrix(cbind(c(cutoff_value1,1),c(cutoff_value2,1))), lower.tail = TRUE))
  
  min_z1 <- z_val_extremums[1,1]
  min_z2 <- z_val_extremums[1,2]
  
  #################################################################
  
  estimate.m0s.pro <- function(p1, p2, B=20){
    
    m <- length(p1)
    
    ##find lambda cutoffs using histogram-based method
    c1 <- calc.cutoff(p1, B=B, max=1)
    c2 <- calc.cutoff(p2, B=B, max=1)
    
    ##estimate m0 for experiment 1
    ind1 <- (p1>=c1)
    m0.1 <- sum(ind1)/(1-c1)
    m0.1 <- min(m0.1, 10000)
    
    ##estimate m0 for experiment 2  
    ind2 <- (p2>=c2)
    m0.2 <- sum(ind2)/(1-c2)
    m0.2 <- min(m0.2, 10000)
    
    ##estimate m00
    ind12 <- ind1 & ind2
    nA <- sum(ind12)
    #pA <- (1-c1)*(1-c2)
    #m00 <- nA/pA
    
    #here, we used converted pvalues to z values to estimnate m00
    # density function for each row of the bivariate z values (x) and 
    # estimated parameters(rho)
    density = function(x, rho)
    {
      sigma = matrix(c(1, rho, rho, 1), 2, 2)
      z = dtmvnorm(x, mean = c(0,0), sigma = sigma, lower = c(min_z1, min_z2))
    }
    
    # log likelihood of the joint densities
    log_likelihood_fn = function(rho){
      
      joint_likelihood = z_val %>% split(.$zvals2) %>% map_dfr(~density(c(.$zvals1,.$zvals2),rho))    
      return(-sum(log(joint_likelihood)))
      
    }
    
    #MLE of the log likelihood function
    optimal_rho = optimize(log_likelihood_fn, lower = -1, upper = 1 ) %>% as.data.frame()
    
    #probability of a random variable greater than cutoff values 
    rho = as.numeric(optimal_rho[1])
    obj_value = optimal_rho[2]
    pA = pmvnorm(lower=c(min_z1, min_z2), upper=c(Inf, Inf), mean=c(0,0), sigma = matrix(c(1, rho, rho, 1), 2, 2))
    
    m00 <- nA/pA
    m00 <- min(m00, 10000)
    
    ##estimate m11
    m11 <- sum(m - m0.1 - m0.2 + m00)
    if (m0.1 == 10000 || m0.2 == 10000 || m00 == 10000){
      m11 = 0
    }
    
    if (m11 < 0){
      m11 = 0
    }
    
    ret <- list()
    ret$ms <- c(obj_value, rho, m, m00, m11)
    names(ret$ms) <- c("obj_value", "optimal_rho", "m", "m00", "m11")
    ret$cutoffs <- c(c1, c2)
    return(ret)
  }
  
  ####################################################################
  
  estimate.m0s <- function(p1, p2, B=20){
    m <- length(p1)
    
    ##find lambda cutoffs using histogram-based method
    c1 <- calc.cutoff(p1, B=B, max=1)
    c2 <- calc.cutoff(p2, B=B, max=1)
    
    ##estimate m0 for experiment 1
    ind1 <- (p1>=c1)
    m0.1 <- sum(ind1)/(1-c1)
    m0.1 <- min(m0.1, 10000)
    
    ##estimate m0 for experiment 2  
    ind2 <- (p2>=c2)
    m0.2 <- sum(ind2)/(1-c2)
    m0.2 <- min(m0.2, 10000)
    
    ##estimate m00
    ind12 <- ind1 & ind2
    nA <- sum(ind12)
    pA <- (1-c1)*(1-c2)
    m00 <- nA/pA
    m00 <- min(m00, 10000)
    
    
    ##estimate m11
    m11 <- sum(m - m0.1 - m0.2 + m00)
    if (m00 == 10000 || m0.1 == 10000 || m0.2 == 10000){
      m11 = 0
    }
    
    if (m11 < 0){
      m11 = 0
    }
    
    ret <- list()
    ret$ms <- c(m00, m11)
    names(ret$ms) <- c("m00.Orr","m11.Orr")
    # ret$cutoffs <- c(c1, c2)
    return(ret)
  }
  
  ###################################################################
  
  limma_ttest <- c(adj_pvalue, adj_limma_pvalue) %>% as.data.frame()
  
  limma_ttest_m00_5 <- limma_ttest[rowSums((limma_ttest[1]>0.05) & (limma_ttest[2]>0.05)), ]
  m00_inter_05 <- nrow(limma_ttest_m00_5)
  
  limma_ttest_m00_1 <- limma_ttest[rowSums((limma_ttest[1]>0.1) & (limma_ttest[2]>0.1)), ]
  m00_inter_1 <- nrow(limma_ttest_m00_1)
  
  limma_ttest_m11_5 <- limma_ttest[rowSums((limma_ttest[1]<=0.05) & (limma_ttest[2]<=0.05)), ]
  m11_inter_05 <- nrow(limma_ttest_m11_5)
  
  limma_ttest_m11_1 <- limma_ttest[rowSums((limma_ttest[1]<=0.1) & (limma_ttest[2]<=0.1)), ]
  m11_inter_1 <- nrow(limma_ttest_m11_1)
  
  ###########################################
  estimator1 = estimate.m0s.pro(pvals1, pvals2, B=20)
  object1 = estimator1$ms %>% as.data.frame()
  
  object11 = estimator1$cutoffs %>% as.data.frame() %>% t()
  colnames(object11) = c("cutoff1", "cutoff2")
  
  estimator2 = estimate.m0s(pvals1, pvals2, B=20)
  object2 = estimator2$ms %>% as.data.frame() %>% t()
  
  object3 = as.data.frame(m00_inter_05)
  object4 = as.data.frame(m00_inter_1)
  
  object5 = as.data.frame(m11_inter_05)
  object6 = as.data.frame(m11_inter_1)
  
  final_result = cbind(object11, object1, object2, object3, object4, object5, object6)
  return(final_result)
}

for (i in 1:50){
  if(i == 1){
    result = simulation(n_samples=4, n_mu=1, n_de=1000)
  }else{
    result = rbind(result, simulation(n_samples=4, n_mu=1, n_de=1000))
  }
} 

write.csv(result, "m11_1k_n4_mu1_result.csv", row.names=FALSE)
