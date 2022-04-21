rm(list=ls())
library("dplyr")

rmse = function(obs, x){
  result = sqrt(mean((obs - x)^2))
  return(result)
}

#m10 = 1000
#m01 =1000

#############################################

m00 = 8500
m11 = 500

df = read.csv("n4_mu1/5h_n4_mu1_resultDTrep.csv")

rmse1 = rmse(m00, df$m00)
rmse2 = rmse(m11, df$m11)
rmse3 = rmse(m11, df$m01)
rmse4 = rmse(m11, df$m10)
rmse5 = rmse(m00, df$m00.Orr)
rmse6 = rmse(m11, df$m11.Orr)
rmse7 = rmse(m11, df$m01.Orr)
rmse8 = rmse(m11, df$m10.Orr)
rmse9 = rmse(m00, df$m00_inter_05)
rmse10 = rmse(m11, df$m11_inter_05)
rmse11 = rmse(m00, df$m00_inter_10)
rmse12 = rmse(m11, df$m11_inter_10)


col_rho = "Estimated Rho"
col_names = c(col_rho, names(df[6:17])) %>% as.data.frame()

mean_col = c(mean(df$optimal_rho), mean(df$m00), mean(df$m11), mean(df$m01), mean(df$m10), mean(df$m00.Orr), 
             mean(df$m11.Orr), mean(df$m01.Orr), mean(df$m10.Orr), mean(df$m00_inter_05), 
             mean(df$m11_inter_05),mean(df$m00_inter_10), mean(df$m11_inter_10)) %>% as.data.frame()
final_rmse = c(NA, rmse1, rmse2, rmse3, rmse4, rmse5, rmse6, rmse7, rmse8, rmse9,
               rmse10, rmse11, rmse12) %>% as.data.frame()
colnames(final_rmse) = ("RMSE")

estimate_df = c(col_names,mean_col, final_rmse) %>% as.data.frame()
colnames(estimate_df) = c("Methods","Mean Estimate", "RMSE")

write.csv(estimate_df, "5h_n4_mu1_rmse.csv")


###################################################################################

m00 = 7000
m11 = 1000
df = read.csv("n4_mu1/1k_n4_mu1_resultDTrep.csv")

rmse1 = rmse(m00, df$m00)
rmse2 = rmse(m11, df$m11)
rmse3 = rmse(m11, df$m01)
rmse4 = rmse(m11, df$m10)
rmse5 = rmse(m00, df$m00.Orr)
rmse6 = rmse(m11, df$m11.Orr)
rmse7 = rmse(m11, df$m01.Orr)
rmse8 = rmse(m11, df$m10.Orr)
rmse9 = rmse(m00, df$m00_inter_05)
rmse10 = rmse(m11, df$m11_inter_05)
rmse11 = rmse(m00, df$m00_inter_10)
rmse12 = rmse(m11, df$m11_inter_10)


col_rho = "Estimated Rho"
col_names = c(col_rho, names(df[6:17])) %>% as.data.frame()

mean_col = c(mean(df$optimal_rho), mean(df$m00), mean(df$m11), mean(df$m01), mean(df$m10), mean(df$m00.Orr), 
             mean(df$m11.Orr), mean(df$m01.Orr), mean(df$m10.Orr), mean(df$m00_inter_05), 
             mean(df$m11_inter_05),mean(df$m00_inter_10), mean(df$m11_inter_10)) %>% as.data.frame()
final_rmse = c(NA, rmse1, rmse2, rmse3, rmse4, rmse5, rmse6, rmse7, rmse8, rmse9,
               rmse10, rmse11, rmse12) %>% as.data.frame()
colnames(final_rmse) = ("RMSE")

estimate_df = c(col_names,mean_col, final_rmse) %>% as.data.frame()
colnames(estimate_df) = c("Methods","Mean Estimate", "RMSE")

write.csv(estimate_df, "1k_n4_mu1_rmse.csv")

###############################################################

m00 = 1000
m11 = 3000

df = read.csv("n4_mu1/3k_n4_mu1_resultDTrep.csv")

rmse1 = rmse(m00, df$m00)
rmse2 = rmse(m11, df$m11)
rmse3 = rmse(m11, df$m01)
rmse4 = rmse(m11, df$m10)
rmse5 = rmse(m00, df$m00.Orr)
rmse6 = rmse(m11, df$m11.Orr)
rmse7 = rmse(m11, df$m01.Orr)
rmse8 = rmse(m11, df$m10.Orr)
rmse9 = rmse(m00, df$m00_inter_05)
rmse10 = rmse(m11, df$m11_inter_05)
rmse11 = rmse(m00, df$m00_inter_10)
rmse12 = rmse(m11, df$m11_inter_10)


col_rho = "Estimated Rho"
col_names = c(col_rho, names(df[6:17])) %>% as.data.frame()

mean_col = c(mean(df$optimal_rho), mean(df$m00), mean(df$m11), mean(df$m01), mean(df$m10), mean(df$m00.Orr), 
             mean(df$m11.Orr), mean(df$m01.Orr), mean(df$m10.Orr), mean(df$m00_inter_05), 
             mean(df$m11_inter_05),mean(df$m00_inter_10), mean(df$m11_inter_10)) %>% as.data.frame()
final_rmse = c(NA, rmse1, rmse2, rmse3, rmse4, rmse5, rmse6, rmse7, rmse8, rmse9,
               rmse10, rmse11, rmse12) %>% as.data.frame()
colnames(final_rmse) = ("RMSE")

estimate_df = c(col_names,mean_col, final_rmse) %>% as.data.frame()
colnames(estimate_df) = c("Methods","Mean Estimate", "RMSE")

write.csv(estimate_df, "3k_n4_mu1_rmse.csv")


