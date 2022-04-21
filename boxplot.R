library("pacman")
p_load("ggplot2", "reshape2", "tidyverse")

plot.fn = function(title, df_name, n_m11){
  df_name = df_name[cols]
  #df_name = data.frame(lapply(df_name, function(x) scale(x, center = FALSE, 
                                                 # scale = max(x, na.rm = TRUE)/n_m11)))
  colnames(df_name) = c("New", "Orr", "Int_05", "Int_10")

  df_name  %>% 
    pivot_longer(cols = everything()) %>% 
    ggplot(aes(x = name, y = value)) +
    geom_boxplot() + xlab("Methods") + ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5)) +
    geom_hline(yintercept = n_m11, linetype = "dashed")
}
