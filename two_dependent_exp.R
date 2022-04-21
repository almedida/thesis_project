simulation3 = function(n_samples, mu_delta, n_de, m11_trt){
  
  m10_genes = sample(selected_genes, n_de) #select 1k genes from the 10k for experiment 1

  #create a data frame using the genes for experiment 1
  m10_exp1 = exp1[m10_genes, ]

  #generate treatment effect for m01 for exp 1
  generate_treatment_effect = function(n_de, genes_df, mu_delta){
    std = apply(genes_df, 1, sd)
    mu = mu_delta*std
    treatment_effect = rnorm(n_de, mu, std)
    
    return(treatment_effect)
  }
  
  #add effects
  m10_exp1_treatment = generate_treatment_effect(n_de, genes_df = m10_exp1, mu_delta)
  treated_m10_exp1 = m10_exp1 + m10_exp1_treatment
  # 
  new_exp1 = exp1 #copy temporary data
  new_exp1[m10_genes, ] = treated_m10_exp1 #replace 1k treated genes
  
  #######################################################################
  rem_selected_genes = setdiff(selected_genes, m10_genes)
  # 
  # #######EXPERIMENT 2 ###########################################################
  # 
  m01_genes = sample(rem_selected_genes, n_de)
  
  m01_exp2 = exp2[m01_genes, ] #create a data frame using genes for experiment 2
   
  # #add effects to m01_Dtrep2
  m01_exp2_treatment = generate_treatment_effect(n_de, genes_df = m01_exp2, mu_delta)
  treated_m01_exp2 = m01_exp2 + m01_exp2_treatment
  # 
  new_exp2 = exp2
  new_exp2[m01_genes, ]=treated_m01_exp2 #replace 1k treated genes
  
  #################EXPERIMENT 1 AND 2#################################
  rem_selected_genes2 = setdiff(rem_selected_genes, m01_genes)
  m11_genes = sample(rem_selected_genes2, m11_trt)
  m11_genes = sample(selected_genes, m11_trt)
  
  m11_exp1 = exp1[m11_genes, ]
  m11_exp2 = exp2[m11_genes, ]
  
  m11_exp1_treatment = generate_treatment_effect(m11_trt, genes_df = m11_exp1, mu_delta)
  treated_m11_exp1 = m11_exp1 + m11_exp1_treatment
  new_exp1[m11_genes, ] = treated_m11_exp1
  
  m11_exp2_treatment = generate_treatment_effect(m11_trt, genes_df = m11_exp2, mu_delta)
  treated_m11_exp2 = m11_exp2 + m11_exp2_treatment
  new_exp2[m11_genes, ]=treated_m11_exp2
  
  #######EXPERIMENT LIMMA ################################
  
  generate_samples = function(treated_exp1_df1, untreated_exp1_df2, treated_exp2_df1, untreated_exp2_df2, n_samples){
    N = ncol(treated_exp1_df1)
    samples  = sample(N, 2*n_samples) #randomly select samples of size 2*n_samples 
    rand_samples1 = samples[1:n_samples]
    rand_samples2 = samples[(n_samples + 1):(2*n_samples)]
    
    treated_exp1 = treated_exp1_df1[rand_samples1]
    untreated_exp1 = untreated_exp1_df2[rand_samples2]
    treated_exp2 = treated_exp2_df1[rand_samples1]
    untreated_exp2 = untreated_exp2_df2[rand_samples2]
    return(list(treated_exp1, untreated_exp1, treated_exp2, untreated_exp2))
  }
  
  limma_sample = generate_samples(treated_exp1_df1 = new_exp1, untreated_exp1_df2 = exp1, 
                                  treated_exp2_df1= new_exp2, untreated_exp2_df2 = exp2, n_samples)
  
  ################## EXPERIMENT 1 LIMMA #########################
  limma1_exp1 = limma_sample[[1]]
  limma2_exp1 = limma_sample[[2]]
  
  ################## EXPERIMENT 2 LIMMA #########################
  limma1_exp2 = limma_sample[[3]]
  limma2_exp2 = limma_sample[[4]]
  
  ###########################################################################
  
  #limma test for experiment 1
  limma_gene_exp1 = c(limma1_exp1, limma2_exp1) %>% as.data.frame()
  
  expr_list_exp1 = factor(
    x = c(rep("exp1_group1", n_samples), rep("exp1_group2", n_samples)),
    levels=c("exp1_group1","exp1_group2")            # Set group 1 to be the first level
  )
  
  exp1_design = model.matrix(~expr_list_exp1)          # Remove the zero
  
  fit = lmFit(limma_gene_exp1, exp1_design)
  fit = eBayes(fit)
  exp1_results = decideTests(fit)
  
  exp1_limma_pval = fit$p.value[, 2] %>% as.data.frame()
  exp1_adj_limma_pval = p.adjust(data.matrix(exp1_limma_pval, "BH")) %>% as.data.frame()
  
  
  ######################################################################################
  
  #######EXPERIMENT 2 ################################
  limma_gene_exp2 = c(limma1_exp2, limma2_exp2) %>% as.data.frame()
  
  expr_list_exp2 = factor(
    x = c(rep("exp2_group1", n_samples), rep("exp2_group2", n_samples)),
    levels=c("exp2_group1","exp2_group2")            # Set group 1 to be the first level
  )
  
  exp2_design = model.matrix(~expr_list_exp2)          # Remove the zero
  
  fit2 = lmFit(limma_gene_exp2, exp2_design)
  fit2 = eBayes(fit2)
  exp2_results = decideTests(fit2)
  
  exp2_limma_pval = fit2$p.value[, 2] %>% as.data.frame()
  exp2_adj_limma_pval = p.adjust(data.matrix(exp2_limma_pval, "BH")) %>% as.data.frame()
  
  ######################################################################################
  
  pval_raw = c(exp1_limma_pval, exp2_limma_pval) %>% as.data.frame()
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
  
  p_vals = pval_raw  %>% filter(exp1_limma_pval >=cutoff_value1, exp2_limma_pval>=cutoff_value2)
  
  #################################################################
  
  z_val = qnorm(as.matrix(p_vals), lower.tail = TRUE) %>% as.data.frame()
  colnames(z_val) = c("zvals1", "zvals2")
  
  zvals1 = z_val[,1] %>% as.data.frame()
  zvals2 = z_val[,2] %>% as.data.frame()
  
  ###############################################################
  
  z_val_extremums = qnorm(as.matrix(cbind(c(cutoff_value1,1),c(cutoff_value2,1))), lower.tail = TRUE) %>% as.data.frame()
  
  min_z1 = z_val_extremums[1,1]
  min_z2 = z_val_extremums[1,2]
  
  #####################################################3
  
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
    
    m01 = m0.1 - m00
    if(m01 < 0){
      m01 = 0
    }
    
    m10 = m0.2 - m00
    if(m10 < 0){
      m10 = 0
    }
    
    
    ret <- list()
    ret$ms <- c(obj_value, rho, m, m00, m11, m01, m10)
    names(ret$ms) <- c("obj_value", "optimal_rho", "m", "m00", "m11", "m01", "m10")
    ret$cutoffs <- c(c1, c2)
    return(ret)
  }
  
  #####################################################################
  
  estimate.m0s.Orr <- function(p1, p2, B=20){
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
    
    m01 = m0.1 - m00
    if(m01 < 0){
      m01 = 0
    }
    m10 = m0.2 - m00
    if(m10 < 0){
      m10 = 0
    }
    
    ret <- list()
    ret$ms <- c(m00, m11, m01, m10)
    names(ret$ms) <- c("m00.Orr","m11.Orr", "m01.Orr", "m10.Orr")
    # ret$cutoffs <- c(c1, c2)
    return(ret)
  }
  
  ############################################################################
  
  m11_intersection = function(df, alpha){
    m11_int = df[rowSums((df[1]<=alpha) & (df[2]<=alpha)), ] 
    m11_int = m11_int %>% nrow()
    return(m11_int)
  }
  
  m00_intersection = function(df, alpha){
    m00_int = df[rowSums((df[1]>alpha) & (df[2]>alpha)), ]
    m00_int =  m00_int %>% nrow()
    return(m00_int)
  }
  
  
  #############################################################################
  limma_ttest = c(exp1_adj_limma_pval, exp2_adj_limma_pval) %>% as.data.frame()
  
  #EE in experiment 1 and 2 controlling FDR at 5% and 10%
  m00_inter_05 = m00_intersection(df=limma_ttest, alpha = 0.05)
  m00_inter_10 = m00_intersection(df=limma_ttest, alpha = 0.1)
  
  #DE in experiment 1 and 2 controlling FDR at 5% and 10%
  m11_inter_05 = m11_intersection(df=limma_ttest, alpha = 0.05)
  m11_inter_10 = m11_intersection(df=limma_ttest, alpha = 0.1)
  
  #####################################################################################
  estimator1 = estimate.m0s.pro(pvals1, pvals2, B=20)
  object1 = estimator1$ms %>% as.data.frame()
  
  object11 = estimator1$cutoffs %>% as.data.frame() %>% t()
  colnames(object11) = c("cutoff1", "cutoff2")
  
  estimator2 = estimate.m0s.Orr(pvals1, pvals2, B=20)
  object2 = estimator2$ms %>% as.data.frame() %>% t()
  
  object3 = data.frame(m00_inter_05)
  object4 = data.frame(m11_inter_05)
  
  #object5 = data.frame(m10_inter_05)
  #object6 = data.frame(m01_inter_05)
  
  object5 = data.frame(m00_inter_10)
  object6 = data.frame(m11_inter_10)
  
  #object9 = data.frame(m10_inter_10)
  #object10 = data.frame(m01_inter_10)
  
  
  
  final_result = cbind(object11, object1, object2, object3, object4, object5, object6)
  return(final_result)
  
}
