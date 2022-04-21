rm(list=ls())
source('../../simulation/boxplot.R')
################################################################################

df_5h = read.csv("n4_mu1/5h_n4_mu1_resultDTrep.csv", header=TRUE, sep =",")
df_1k = read.csv("n4_mu1/1k_n4_mu1_resultDTrep.csv", header=TRUE, sep =",")
df_3k = read.csv("n4_mu1/3k_n4_mu1_resultDTrep.csv", header=TRUE, sep =",")
cols = c("m11","m11.Orr", "m11_inter_05", "m11_inter_10")

plot.fn(title = "n = 4, m11 = 500, mu = 1", df_name = df_5h, n_m11 = 500)
plot.fn(title = "n = 4, m11 = 1000, mu = 1", df_name = df_1k, n_m11 = 1000)
plot.fn(title = "n = 4, m11 = 3000, mu = 1", df_name = df_3k, n_m11 = 3000)

#################################################################

df_5h = read.csv("n4_mu2/5h_n4_mu2_resultDTrep.csv", header=TRUE, sep =",")
df_1k = read.csv("n4_mu2/1k_n4_mu2_resultDTrep.csv", header=TRUE, sep =",")
df_3k = read.csv("n4_mu2/3k_n4_mu2_resultDTrep.csv", header=TRUE, sep =",")
cols = c("m11","m11.Orr", "m11_inter_05", "m11_inter_10")

plot.fn(title = "n = 4, m11 = 500, mu = 2", df_name = df_5h, n_m11 = 500)
plot.fn(title = "n = 4, m11 = 1000, mu = 2", df_name = df_1k, n_m11 = 1000)
plot.fn(title = "n = 4, m11 = 3000, mu = 2", df_name = df_3k, n_m11 = 3000)

##########################################################################

df_5h = read.csv("n10_mu1/5h_n10_mu1_resultDTrep.csv", header=TRUE, sep =",")
df_1k = read.csv("n10_mu1/1k_n10_mu1_resultDTrep.csv", header=TRUE, sep =",")
df_3k = read.csv("n10_mu1/3k_n10_mu1_resultDTrep.csv", header=TRUE, sep =",")
cols = c("m11","m11.Orr", "m11_inter_05", "m11_inter_10")

plot.fn(title = "n = 10, m11 = 500, mu = 1", df_name = df_5h, n_m11 = 500)
plot.fn(title = "n = 10, m11 = 1000, mu = 1", df_name = df_1k, n_m11 = 1000)
plot.fn(title = "n = 10, m11 = 3000, mu = 1", df_name = df_3k, n_m11 = 3000)

##########################################################################

df_5h = read.csv("n10_mu2/5h_n10_mu2_resultDTrep.csv", header=TRUE, sep =",")
df_1k = read.csv("n10_mu2/1k_n10_mu2_resultDTrep.csv", header=TRUE, sep =",")
df_3k = read.csv("n10_mu2/3k_n10_mu2_resultDTrep.csv", header=TRUE, sep =",")
cols = c("m11","m11.Orr", "m11_inter_05", "m11_inter_10")

plot.fn(title = "n = 10, m11 = 500, mu = 2", df_name = df_5h, n_m11 = 500)
plot.fn(title = "n = 10, m11 = 1000, mu = 2", df_name = df_1k, n_m11 = 1000)
plot.fn(title = "n = 10, m11 = 3000, mu = 2", df_name = df_3k, n_m11 = 3000)
