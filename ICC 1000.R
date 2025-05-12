## Extract cpg loadings from PC
library(tidyverse)
library(dplyr)
library(tidyr)
library(psych)

#set names of original dataset as rownames
RP_tenper_random_namesfix <- RP_tenper_random_namesfix %>%
  column_to_rownames("X")
view(RP_tenper_random_namesfix)

#find out if the cpg sites in RP_tenper... are in the same order as the PCA
#This code proves that, since RP_numeric was used for the PCA
identical(RP_tenper_random_namesfix$S1_V2, RP_numeric$S1_V2)

rotation_table <- pca$rotation
nrow(rotation_table) #rotation_table = PC values for each cpg site per PC
rownames(rotation_table) <- rownames(RP_tenper_random_namesfix)
view(rotation_table) #now, the CpG site corresponding to each PCA value is displayed in rotation_table

#making a dataframe for PC values per CpG site in descending order for PC1
pc1_rotation_table <- rotation_table %>%
  as.data.frame() %>%
  dplyr::select(PC1) %>%
  mutate(PC1 = abs(PC1)) %>%
  arrange(desc(PC1))
any(pc1_rotation_table$PC1 < 0)
view(pc1_rotation_table)
top1000_pc1_cpgs <- rownames(pc1_rotation_table)[1:1000]
head(top1000_pc1_cpgs)
matching_pc1000 <- intersect(top1000_pc1_cpgs, rownames(RP_tenper_random_namesfix))
View(matching_pc1000)
beta_pc1000 <- RP_tenper_random_namesfix[matching_pc1000, ]
View(beta_pc1000)
beta_pc1000_t <- t(beta_pc1000)
View(beta_pc1000_t)
ids_pc1000 <- rownames(beta_pc1000_t)
meta_pc <- data.frame(
  Sample = ids_pc1000,
  Subject = sub("_V[12]", "", ids_pc1000),
  Visit = sub(".*_V", "V", ids_pc1000))
View(meta_pc)
beta_df_pc1000 <- cbind(meta_pc, beta_pc1000_t)
View(beta_df_pc1000)

library(dplyr)
library(tidyr)
library(irr)

df_ICC1000 <- data.frame(cpg_id = top1000_pc1_cpgs,
                      ICC = NA,
                      p_value = NA)
for(cpg in top1000_pc1_cpgs){
  print(cpg)
  current_df1000 <- beta_df_pc1000 %>%
    dplyr::select(Subject, Visit, !!sym(cpg))%>%
    pivot_wider(names_from = Visit, values_from = !!sym(cpg))%>% 
    dplyr::select(-Subject)  # Remove the Subject column
  icc_result1000 <- irr::icc(current_df1000, model = "twoway", type = "consistency", unit = "single")
  icc_value1000 <- icc_result1000$value
  p_val1000 <- icc_result1000$p.value
  
  df_ICC1000[df_ICC1000$cpg_id == cpg, "ICC"] <- icc_value1000
  df_ICC1000[df_ICC1000$cpg_id == cpg, "p_value"] <- p_val1000
}
View(df_ICC1000)
mean(df_ICC1000$ICC) #0.1918788








pc2_rotation_table <- rotation_table %>%
  as.data.frame() %>%
  dplyr::select(PC2) %>%
  mutate(PC2 = abs(PC2)) %>%
  arrange(desc(PC2))
view(pc2_rotation_table)
top1000_pc2_cpgs <- rownames(pc2_rotation_table)[1:1000]
head(top1000_pc2_cpgs)
matching_pc1000 <- intersect(top1000_pc2_cpgs, rownames(RP_tenper_random_namesfix))
View(matching_pc1000)
beta_pc2000 <- RP_tenper_random_namesfix[matching_pc1000, ]
View(beta_pc2000)
beta_pc2000_t <- t(beta_pc2000)
View(beta_pc2000_t)
ids_pc2000 <- rownames(beta_pc2000_t)
beta_df_pc2000 <- cbind(meta_pc, beta_pc2000_t)
View(beta_df_pc2000)
library(dplyr)
library(tidyr)
library(irr)

df_ICC2000 <- data.frame(cpg_id = top1000_pc2_cpgs,
                      ICC = NA,
                      p_value = NA)
for(cpg in top1000_pc2_cpgs){
  print(cpg)
  current_df2000 <- beta_df_pc2000 %>%
    dplyr::select(Subject, Visit, !!sym(cpg))%>%
    pivot_wider(names_from = Visit, values_from = !!sym(cpg))%>% 
    dplyr::select(-Subject)  # Remove the Subject column
  icc_result2000 <- irr::icc(current_df2000, model = "twoway", type = "consistency", unit = "single")
  icc_value2000 <- icc_result2000$value
  p_val2000 <- icc_result2000$p.value
  
  df_ICC2000[df_ICC2000$cpg_id == cpg, "ICC"] <- icc_value2000
  df_ICC2000[df_ICC2000$cpg_id == cpg, "p_value"] <- p_val2000
}
View(df_ICC2000)
mean(df_ICC2000$ICC) #0.8587877







pc3_rotation_table <- rotation_table %>%
  as.data.frame() %>%
  dplyr::select(PC3) %>%
  mutate(PC3 = abs(PC3)) %>%
  arrange(desc(PC3))
view(pc3_rotation_table)
top1000_pc3_cpgs <- rownames(pc3_rotation_table)[1:1000]
head(top1000_pc3_cpgs)
matching_pc3000 <- intersect(top1000_pc3_cpgs, rownames(RP_tenper_random_namesfix))
View(matching_pc3000)
beta_pc3000 <- RP_tenper_random_namesfix[matching_pc3000, ]
View(beta_pc3000)
beta_pc3000_t <- t(beta_pc3000)
View(beta_pc3000_t)
ids_pc3000 <- rownames(beta_pc3000_t)
beta_df_pc3000 <- cbind(meta_pc, beta_pc3000_t)
View(beta_df_pc3000)
library(dplyr)
library(tidyr)
library(irr)

df_ICC3000 <- data.frame(cpg_id = top1000_pc3_cpgs,
                      ICC = NA,
                      p_value = NA)
for(cpg in top1000_pc3_cpgs){
  print(cpg)
  current_df3000 <- beta_df_pc3000 %>%
    dplyr::select(Subject, Visit, !!sym(cpg))%>%
    pivot_wider(names_from = Visit, values_from = !!sym(cpg))%>% 
    dplyr::select(-Subject)  # Remove the Subject column
  icc_result3000 <- irr::icc(current_df3000, model = "twoway", type = "consistency", unit = "single")
  icc_value3000 <- icc_result3000$value
  p_val3000 <- icc_result3000$p.value
  
  df_ICC3000[df_ICC3000$cpg_id == cpg, "ICC"] <- icc_value3000
  df_ICC3000[df_ICC3000$cpg_id == cpg, "p_value"] <- p_val3000
}
View(df_ICC3000)
mean(df_ICC3000$ICC) #0.391255









pc4_rotation_table <- rotation_table %>%
  as.data.frame() %>%
  dplyr::select(PC4) %>%
  mutate(PC4 = abs(PC4)) %>%
  arrange(desc(PC4))
view(pc4_rotation_table)
top1000_pc4_cpgs <- rownames(pc4_rotation_table)[1:1000]
head(top1000_pc4_cpgs)
matching_pc4000 <- intersect(top1000_pc4_cpgs, rownames(RP_tenper_random_namesfix))
View(matching_pc4000)
beta_pc4000 <- RP_tenper_random_namesfix[matching_pc4000, ]
View(beta_pc4000)
beta_pc4000_t <- t(beta_pc4000)
View(beta_pc4000_t)
ids_pc4000 <- rownames(beta_pc4000_t)
beta_df_pc4000 <- cbind(meta_pc, beta_pc4000_t)
View(beta_df_pc4000)
library(dplyr)
library(tidyr)
library(irr)

df_ICC4000 <- data.frame(cpg_id = top1000_pc4_cpgs,
                      ICC = NA,
                      p_value = NA)
for(cpg in top1000_pc4_cpgs){
  print(cpg)
  current_df4000 <- beta_df_pc4000 %>%
    dplyr::select(Subject, Visit, !!sym(cpg))%>%
    pivot_wider(names_from = Visit, values_from = !!sym(cpg))%>% 
    dplyr::select(-Subject)  # Remove the Subject column
  icc_result4000 <- irr::icc(current_df4000, model = "twoway", type = "consistency", unit = "single")
  icc_value4000 <- icc_result4000$value
  p_val4000 <- icc_result4000$p.value
  
  df_ICC4000[df_ICC4000$cpg_id == cpg, "ICC"] <- icc_value4000
  df_ICC4000[df_ICC4000$cpg_id == cpg, "p_value"] <- p_val4000
}
View(df_ICC4000)
mean(df_ICC4000$ICC) #0.2320076














pc5_rotation_table <- rotation_table %>%
  as.data.frame() %>%
  dplyr::select(PC5) %>%
  mutate(PC5 = abs(PC5)) %>%
  arrange(desc(PC5))
view(pc5_rotation_table)
top1000_pc5_cpgs <- rownames(pc5_rotation_table)[1:1000]
head(top1000_pc5_cpgs)
matching_pc5000 <- intersect(top1000_pc5_cpgs, rownames(RP_tenper_random_namesfix))
View(matching_pc5000)
beta_pc5000 <- RP_tenper_random_namesfix[matching_pc5000, ]
View(beta_pc5000)
beta_pc5000_t <- t(beta_pc5000)
View(beta_pc5000_t)
ids_pc5000 <- rownames(beta_pc5000_t)
beta_df_pc5000 <- cbind(meta_pc, beta_pc5000_t)
View(beta_df_pc5000)
library(dplyr)
library(tidyr)
library(irr)

df_ICC5000 <- data.frame(cpg_id = top1000_pc5_cpgs,
                      ICC = NA,
                      p_value = NA)
for(cpg in top1000_pc5_cpgs){
  print(cpg)
  current_df5000 <- beta_df_pc5000 %>%
    dplyr::select(Subject, Visit, !!sym(cpg))%>%
    pivot_wider(names_from = Visit, values_from = !!sym(cpg))%>% 
    dplyr::select(-Subject)  # Remove the Subject column
  icc_result5000 <- irr::icc(current_df5000, model = "twoway", type = "consistency", unit = "single")
  icc_value5000 <- icc_result5000$value
  p_val5000 <- icc_result5000$p.value
  
  df_ICC5000[df_ICC5000$cpg_id == cpg, "ICC"] <- icc_value5000
  df_ICC5000[df_ICC5000$cpg_id == cpg, "p_value"] <- p_val5000
}
View(df_ICC5000)
mean(df_ICC5000$ICC) #0.4235237









pc6_rotation_table <- rotation_table %>%
  as.data.frame() %>%
  dplyr::select(PC6) %>%
  mutate(PC6 = abs(PC6)) %>%
  arrange(desc(PC6))
view(pc6_rotation_table)
top1000_pc6_cpgs <- rownames(pc6_rotation_table)[1:1000]
head(top1000_pc6_cpgs)
matching_pc6000 <- intersect(top1000_pc6_cpgs, rownames(RP_tenper_random_namesfix))
View(matching_pc6000)
beta_pc6000 <- RP_tenper_random_namesfix[matching_pc6000, ]
View(beta_pc6000)
beta_pc6000_t <- t(beta_pc6000)
View(beta_pc6000_t)
ids_pc6000 <- rownames(beta_pc6000_t)
beta_df_pc6000 <- cbind(meta_pc, beta_pc6000_t)
View(beta_df_pc6)
library(dplyr)
library(tidyr)
library(irr)

df_ICC6000 <- data.frame(cpg_id = top1000_pc6_cpgs,
                      ICC = NA,
                      p_value = NA)
for(cpg in top1000_pc6_cpgs){
  print(cpg)
  current_df6000 <- beta_df_pc6000 %>%
    dplyr::select(Subject, Visit, !!sym(cpg))%>%
    pivot_wider(names_from = Visit, values_from = !!sym(cpg))%>% 
    dplyr::select(-Subject)  # Remove the Subject column
  icc_result6000 <- irr::icc(current_df6000, model = "twoway", type = "consistency", unit = "single")
  icc_value6000 <- icc_result6000$value
  p_val6000 <- icc_result6000$p.value
  
  df_ICC6000[df_ICC6000$cpg_id == cpg, "ICC"] <- icc_value6000
  df_ICC6000[df_ICC6000$cpg_id == cpg, "p_value"] <- p_val6000
}
View(df_ICC6000)
mean(df_ICC6000$ICC) #0.3792592


ICC_1000_means<- c(0.1918788,0.8587877,0.391255,0.2320076,0.4235237,0.3792592)
PC_ratios<- c()

