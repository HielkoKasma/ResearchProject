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
top100_pc1_cpgs <- rownames(pc1_rotation_table)[1:100]
head(top100_pc1_cpgs)
matching_pc1 <- intersect(top100_pc1_cpgs, rownames(RP_tenper_random_namesfix))
View(matching_pc1)
beta_pc1 <- RP_tenper_random_namesfix[matching_pc1, ]
View(beta_pc1)
beta_pc1_t <- t(beta_pc1)
View(beta_pc1_t)
ids_pc1 <- rownames(beta_pc1_t)
meta_pc <- data.frame(
  Sample = ids_pc1,
  Subject = sub("_V[12]", "", ids_pc1),
  Visit = sub(".*_V", "V", ids_pc1))
View(meta_pc)
beta_df_pc1 <- cbind(meta_pc, beta_pc1_t)
View(beta_df_pc1)





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
