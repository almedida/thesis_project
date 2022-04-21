source("dataprep.R")
source("methods.R")

#c1 <- calc.cutoff(adj_pval_liver, B=20, max=1)
#c2 <- calc.cutoff(adj_liver_limma_pval, B=20, max=1)

#cutoff = cbind(c(c1), c(c2))
#colnames(cutoff) = c("cutoff_value1", "cutoff_value2")

###################################################################
#pval_raw = c(pvals1, pvals2) %>% as.data.frame()
#df1 = pvals1 %>% as.data.frame()
#df2 = pvals2 %>% as.data.frame()
#p_vals = pval_raw  %>% filter(df1 >=c1, df2>=c2)


#proposed method and Orr for t test, limma
#ttest
result1 = estimate.m0s.pro(pvals1, pvals2)
result2 = estimate.m0s.Orr(pvals1, pvals2)

object0 = result1$cutoffs %>% as.data.frame() %>% t()
colnames(object0) = c("cutoff1", "cutoff2")
object1 = result1$ms %>% as.data.frame() %>% t()
object2 = result2$ms %>% as.data.frame() %>% t()

result_ttest = cbind(object0, object1, object2)

#limma
result3 = estimate.m0s.pro(kidney_limma_pval, liver_limma_pval)
result4 = estimate.m0s.Orr(kidney_limma_pval, liver_limma_pval)

object3 = result3$cutoffs %>% as.data.frame() %>% t()
colnames(object3) = c("cutoff1", "cutoff2")
object4 = result3$ms %>% as.data.frame() %>% t()
object5 = result4$ms %>% as.data.frame() %>% t()

result_limma = cbind(object3, object4, object5)

#test vs limma
#kidney
result5 = estimate.m0s.pro(pvals1, kidney_limma_pval)
result6 = estimate.m0s.Orr(pvals1, kidney_limma_pval)

object6 = result5$cutoffs %>% as.data.frame() %>% t()
colnames(object6) = c("cutoff1", "cutoff2")
object7 = result5$ms %>% as.data.frame() %>% t()
object8 = result6$ms %>% as.data.frame() %>% t()

result_ttest_limma_kidney = cbind(object6, object7, object8)

#liver
result7 = estimate.m0s.pro(pvals2, liver_limma_pval)
result8 = estimate.m0s.Orr(pvals2, liver_limma_pval)

object9 = result7$cutoffs %>% as.data.frame() %>% t()
colnames(object9) = c("cutoff1", "cutoff2")
object10 = result7$ms %>% as.data.frame() %>% t()
object11 = result8$ms %>% as.data.frame() %>% t()

result_ttest_limma_liver = cbind(object9, object10, object11)

# r = cbind(result_ttest, result_limma, result_ttest_limma_kidney, 
#       result_ttest_limma_liver)
############################################################
#intersection for t test, limma
#ttest
result9 = m11_intersection(adj_pvals, 0.05)
result10 = m11_intersection(adj_pvals, 0.1)
result99 = m00_intersection(adj_pvals, 0.05)
result1010 = m00_intersection(adj_pvals, 0.1)

ttest_05_10 = cbind(result9, result10, result99, result1010)
colnames(ttest_05_10) = c("m11_ttest_05", "m11_ttest_10", "m00_ttest_05", "m00_ttest_10")

#limma
result11 = m11_intersection(adj_limmas, 0.05)
result12 = m11_intersection(adj_limmas, 0.1)
result111 = m00_intersection(adj_limmas, 0.05)
result122 = m00_intersection(adj_limmas, 0.1)

limma_05_10 = cbind(result11, result12, result111, result122)
colnames(limma_05_10) = c("m11_limma_05", "m11_limma_10", "m00_limma_05", "m00_limma_10")

# #t test vs limma
# #kidney
adj_test_limma_kidney = c(adj_pval_kidney, adj_kidney_limma_pval) %>% as.data.frame()
result13 = m11_intersection(adj_test_limma_kidney, 0.05)
result14 = m11_intersection(adj_test_limma_kidney, 0.1)
result133 = m00_intersection(adj_test_limma_kidney, 0.05)
result144 = m00_intersection(adj_test_limma_kidney, 0.1)

limma_vs_ttest_05_10_kidney = cbind(result13, result14, result133, result144)
colnames(limma_vs_ttest_05_10_kidney) = c("m11_limma_vs_ttest_kidney_05", "m11_limma_vs_ttest_kidney_10",
                                   "m00_limma_vs_ttest_kidney_05", "m00_limma_vs_ttest_kidney_10")
 
# #liver
adj_test_limma_liver = c(adj_pval_liver, adj_liver_limma_pval) %>% as.data.frame()
result15 = m11_intersection(adj_test_limma_liver, 0.05)
result16 = m11_intersection(adj_test_limma_liver, 0.1)
result155 = m00_intersection(adj_test_limma_liver, 0.05)
result166 = m00_intersection(adj_test_limma_liver, 0.1)

limma_vs_ttest_05_10_liver = cbind(result15, result16, result155, result166)
colnames(limma_vs_ttest_05_10_liver) = c("m11_limma_vs_ttest_liver_05", "m11_limma_vs_ttest_liver_10",
                                "m00_limma_vs_ttest_liver_05", "m00_limma_vs_ttest_liver_10")
