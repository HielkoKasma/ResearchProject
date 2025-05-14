## Extract cpg loadings from PC
library(tidyverse)
library(dplyr)
library(tidyr)
library(psych)

#set names of original dataset as rownames
RP_tenper_random_namesfix <- RP_tenper_random_namesfix %>%
  column_to_rownames("X")
#view(RP_tenper_random_namesfix)

#find out if the cpg sites in RP_tenper... are in the same order as the PCA
#This code proves that, since RP_numeric was used for the PCA
identical(RP_tenper_random_namesfix$S1_V2, RP_numeric$S1_V2)

rotation_table <- pca$rotation
nrow(rotation_table) #rotation_table = PC values for each cpg site per PC
rownames(rotation_table) <- rownames(RP_tenper_random_namesfix)
#view(rotation_table) #now, the CpG site corresponding to each PCA value is displayed in rotation_table

#making a dataframe for PC values per CpG site in descending order for PC1
pc1_rotation_table <- rotation_table %>%
  as.data.frame() %>%
  dplyr::select(PC1) %>%
  mutate(PC1 = abs(PC1)) %>%
  arrange(desc(PC1))
any(pc1_rotation_table$PC1 < 0)
#view(pc1_rotation_table)
top100_pc1_cpgs <- rownames(pc1_rotation_table)[1:100]
head(top100_pc1_cpgs)
matching_pc1 <- intersect(top100_pc1_cpgs, rownames(RP_tenper_random_namesfix))
#View(matching_pc1)
beta_pc1 <- RP_tenper_random_namesfix[matching_pc1, ]
#View(beta_pc1)
beta_pc1_t <- t(beta_pc1)
#View(beta_pc1_t)
ids_pc1 <- rownames(beta_pc1_t)
meta_pc <- data.frame(
  Sample = ids_pc1,
  Subject = sub("_V[12]", "", ids_pc1),
  Visit = sub(".*_V", "V", ids_pc1))
#View(meta_pc)
beta_df_pc1 <- cbind(meta_pc, beta_pc1_t)
#View(beta_df_pc1)

library(dplyr)
library(tidyr)
library(irr)

df_ICC1 <- data.frame(cpg_id = top100_pc1_cpgs,
                      ICC = NA,
                      p_value = NA)
for(cpg in top100_pc1_cpgs){
  print(cpg)
  current_df1 <- beta_df_pc1 %>%
    dplyr::select(Subject, Visit, !!sym(cpg))%>%
    pivot_wider(names_from = Visit, values_from = !!sym(cpg))%>% 
    dplyr::select(-Subject)  # Remove the Subject column
  icc_result1 <- irr::icc(current_df1, model = "twoway", type = "consistency", unit = "single")
  icc_value1 <- icc_result1$value
  p_val1 <- icc_result1$p.value
  
  df_ICC1[df_ICC1$cpg_id == cpg, "ICC"] <- icc_value1
  df_ICC1[df_ICC1$cpg_id == cpg, "p_value"] <- p_val1
}
View(df_ICC1)









pc2_rotation_table <- rotation_table %>%
  as.data.frame() %>%
  dplyr::select(PC2) %>%
  mutate(PC2 = abs(PC2)) %>%
  arrange(desc(PC2))
view(pc2_rotation_table)
top100_pc2_cpgs <- rownames(pc2_rotation_table)[1:100]
head(top100_pc2_cpgs)
matching_pc2 <- intersect(top100_pc2_cpgs, rownames(RP_tenper_random_namesfix))
View(matching_pc2)
beta_pc2 <- RP_tenper_random_namesfix[matching_pc2, ]
View(beta_pc2)
beta_pc2_t <- t(beta_pc2)
View(beta_pc2_t)
ids_pc2 <- rownames(beta_pc2_t)
beta_df_pc2 <- cbind(meta_pc, beta_pc2_t)
View(beta_df_pc2)
library(dplyr)
library(tidyr)
library(irr)

df_ICC2 <- data.frame(cpg_id = top100_pc2_cpgs,
                      ICC = NA,
                      p_value = NA)
for(cpg in top100_pc2_cpgs){
  print(cpg)
  current_df2 <- beta_df_pc2 %>%
    dplyr::select(Subject, Visit, !!sym(cpg))%>%
    pivot_wider(names_from = Visit, values_from = !!sym(cpg))%>% 
    dplyr::select(-Subject)  # Remove the Subject column
  icc_result2 <- irr::icc(current_df2, model = "twoway", type = "consistency", unit = "single")
  icc_value2 <- icc_result2$value
  p_val2 <- icc_result2$p.value
  
  df_ICC2[df_ICC2$cpg_id == cpg, "ICC"] <- icc_value2
  df_ICC2[df_ICC2$cpg_id == cpg, "p_value"] <- p_val2
}
View(df_ICC2)








pc3_rotation_table <- rotation_table %>%
  as.data.frame() %>%
  dplyr::select(PC3) %>%
  mutate(PC3 = abs(PC3)) %>%
  arrange(desc(PC3))
view(pc3_rotation_table)
top100_pc3_cpgs <- rownames(pc3_rotation_table)[1:100]
head(top100_pc3_cpgs)
matching_pc3 <- intersect(top100_pc3_cpgs, rownames(RP_tenper_random_namesfix))
View(matching_pc3)
beta_pc3 <- RP_tenper_random_namesfix[matching_pc3, ]
View(beta_pc3)
beta_pc3_t <- t(beta_pc3)
View(beta_pc3_t)
ids_pc3 <- rownames(beta_pc3_t)
beta_df_pc3 <- cbind(meta_pc, beta_pc3_t)
View(beta_df_pc3)
library(dplyr)
library(tidyr)
library(irr)

df_ICC3 <- data.frame(cpg_id = top100_pc3_cpgs,
                      ICC = NA,
                      p_value = NA)
for(cpg in top100_pc3_cpgs){
  print(cpg)
  current_df3 <- beta_df_pc3 %>%
    dplyr::select(Subject, Visit, !!sym(cpg))%>%
    pivot_wider(names_from = Visit, values_from = !!sym(cpg))%>% 
    dplyr::select(-Subject)  # Remove the Subject column
  icc_result3 <- irr::icc(current_df3, model = "twoway", type = "consistency", unit = "single")
  icc_value3 <- icc_result3$value
  p_val3 <- icc_result3$p.value
  
  df_ICC3[df_ICC3$cpg_id == cpg, "ICC"] <- icc_value3
  df_ICC3[df_ICC3$cpg_id == cpg, "p_value"] <- p_val3
}
View(df_ICC3)










pc4_rotation_table <- rotation_table %>%
  as.data.frame() %>%
  dplyr::select(PC4) %>%
  mutate(PC4 = abs(PC4)) %>%
  arrange(desc(PC4))
view(pc4_rotation_table)
top100_pc4_cpgs <- rownames(pc4_rotation_table)[1:100]
head(top100_pc4_cpgs)
matching_pc4 <- intersect(top100_pc4_cpgs, rownames(RP_tenper_random_namesfix))
View(matching_pc4)
beta_pc4 <- RP_tenper_random_namesfix[matching_pc4, ]
View(beta_pc4)
beta_pc4_t <- t(beta_pc4)
View(beta_pc4_t)
ids_pc4 <- rownames(beta_pc4_t)
beta_df_pc4 <- cbind(meta_pc, beta_pc4_t)
View(beta_df_pc4)
library(dplyr)
library(tidyr)
library(irr)

df_ICC4 <- data.frame(cpg_id = top100_pc4_cpgs,
                      ICC = NA,
                      p_value = NA)
for(cpg in top100_pc4_cpgs){
  print(cpg)
  current_df4 <- beta_df_pc4 %>%
    dplyr::select(Subject, Visit, !!sym(cpg))%>%
    pivot_wider(names_from = Visit, values_from = !!sym(cpg))%>% 
    dplyr::select(-Subject)  # Remove the Subject column
  icc_result4 <- irr::icc(current_df4, model = "twoway", type = "consistency", unit = "single")
  icc_value4 <- icc_result4$value
  p_val4 <- icc_result4$p.value
  
  df_ICC4[df_ICC4$cpg_id == cpg, "ICC"] <- icc_value4
  df_ICC4[df_ICC4$cpg_id == cpg, "p_value"] <- p_val4
}
View(df_ICC4)















pc5_rotation_table <- rotation_table %>%
  as.data.frame() %>%
  dplyr::select(PC5) %>%
  mutate(PC5 = abs(PC5)) %>%
  arrange(desc(PC5))
view(pc5_rotation_table)
top100_pc5_cpgs <- rownames(pc5_rotation_table)[1:100]
head(top100_pc5_cpgs)
matching_pc5 <- intersect(top100_pc5_cpgs, rownames(RP_tenper_random_namesfix))
View(matching_pc5)
beta_pc5 <- RP_tenper_random_namesfix[matching_pc5, ]
View(beta_pc5)
beta_pc5_t <- t(beta_pc5)
View(beta_pc5_t)
ids_pc5 <- rownames(beta_pc5_t)
beta_df_pc5 <- cbind(meta_pc, beta_pc5_t)
View(beta_df_pc5)
library(dplyr)
library(tidyr)
library(irr)

df_ICC5 <- data.frame(cpg_id = top100_pc5_cpgs,
                      ICC = NA,
                      p_value = NA)
for(cpg in top100_pc5_cpgs){
  print(cpg)
  current_df5 <- beta_df_pc5 %>%
    dplyr::select(Subject, Visit, !!sym(cpg))%>%
    pivot_wider(names_from = Visit, values_from = !!sym(cpg))%>% 
    dplyr::select(-Subject)  # Remove the Subject column
  icc_result5 <- irr::icc(current_df5, model = "twoway", type = "consistency", unit = "single")
  icc_value5 <- icc_result5$value
  p_val5 <- icc_result5$p.value
  
  df_ICC5[df_ICC5$cpg_id == cpg, "ICC"] <- icc_value5
  df_ICC5[df_ICC5$cpg_id == cpg, "p_value"] <- p_val5
}
View(df_ICC5)










pc6_rotation_table <- rotation_table %>%
  as.data.frame() %>%
  dplyr::select(PC6) %>%
  mutate(PC6 = abs(PC6)) %>%
  arrange(desc(PC6))
view(pc6_rotation_table)
top100_pc6_cpgs <- rownames(pc6_rotation_table)[1:100]
head(top100_pc6_cpgs)
matching_pc6 <- intersect(top100_pc6_cpgs, rownames(RP_tenper_random_namesfix))
View(matching_pc6)
beta_pc6 <- RP_tenper_random_namesfix[matching_pc6, ]
View(beta_pc6)
beta_pc6_t <- t(beta_pc6)
View(beta_pc6_t)
ids_pc6 <- rownames(beta_pc6_t)
beta_df_pc6 <- cbind(meta_pc, beta_pc6_t)
View(beta_df_pc6)
library(dplyr)
library(tidyr)
library(irr)

df_ICC6 <- data.frame(cpg_id = top100_pc6_cpgs,
                      ICC = NA,
                      p_value = NA)
for(cpg in top100_pc6_cpgs){
  print(cpg)
  current_df6 <- beta_df_pc6 %>%
    dplyr::select(Subject, Visit, !!sym(cpg))%>%
    pivot_wider(names_from = Visit, values_from = !!sym(cpg))%>% 
    dplyr::select(-Subject)  # Remove the Subject column
  icc_result6 <- irr::icc(current_df6, model = "twoway", type = "consistency", unit = "single")
  icc_value6 <- icc_result6$value
  p_val6 <- icc_result6$p.value
  
  df_ICC6[df_ICC6$cpg_id == cpg, "ICC"] <- icc_value6
  df_ICC6[df_ICC6$cpg_id == cpg, "p_value"] <- p_val6
}
View(df_ICC6)

results_df

#Plotting ICC vs PC
mean_ICC1 <- mean(df_ICC1$ICC)
mean_ICC2 <- mean(df_ICC2$ICC)
mean_ICC3 <- mean(df_ICC3$ICC)
mean_ICC4 <- mean(df_ICC4$ICC)
mean_ICC5 <- mean(df_ICC5$ICC)
mean_ICC6 <- mean(df_ICC6$ICC)

ICC_means <- c(mean_ICC1,
               mean_ICC2,
               mean_ICC3,
               mean_ICC4,
               mean_ICC5,
               mean_ICC6)
ICC_means
ICC_means_df <- data.frame(t(ICC_means))
ICC_means_df <- data.frame(t(ICC_1000_means))
colnames(ICC_means_df) <- c("ICC1", "ICC2", "ICC3", "ICC4", "ICC5", "ICC6")
rownames(ICC_means_df) <- c("Mean_value")
ICC_means_df

PC_stability_dfform <- data.frame(results_df[5, ])

PC_stability_values <- t(PC_stability_dfform)
colnames(PC_stability_values) <- c("Stability")
ICC_mean_values <- t(ICC_means_df)

PC_vs_ICC <- data.frame(PC_stability_values, ICC_mean_values)

labels_PCAvsICC <- c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6")

library(ggrepel)
library(ggplot2)
ggplot(PC_vs_ICC, aes(x = Stability, y = Mean_value)) +
  geom_point() +
  geom_text_repel(aes(label = labels_PCAvsICC), 
                  size = 3.5, 
                  max.overlaps = Inf) +
  theme_bw() +
  labs(x = "PC Stability Ratio", y = "ICC Mean Value", 
       title = "Stability PCA vs ICC")

#flipping the ratio's, so that ratio of stability is non-self/self
#this results in higher stability = higher ratio, instead of higer stability = lower ratio
stability_ratio_correct <- results_df[3,]/results_df[1,]

#adding new ratios to last row of results_df and removing redundant/incorrect rows
results_df <- rbind(results_df, stability_ratio_correct)
rownames(results_df)[7] <- "Stability ratio (non-self/self)"
results_df <- results_df[-c(5, 6), ]
view(results_df)


ICC100_means<- c(0.1855919,0.9861335,0.4094824,0.4051611,0.565348,0.5740619)

ICC_1000_means<- c(0.1918788,0.8587877,0.391255,0.2320076,0.4235237,0.3792592)


