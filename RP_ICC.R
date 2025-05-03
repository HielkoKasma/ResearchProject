## Extract cpg loadings from PC
library(tidyverse)

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

pc2_rotation_table <- rotation_table %>%
  as.data.frame() %>%
  dplyr::select(PC2) %>%
  mutate(PC2 = abs(PC2)) %>%
  arrange(desc(PC2))
view(pc2_rotation_table)

pc3_rotation_table <- rotation_table %>%
  as.data.frame() %>%
  dplyr::select(PC3) %>%
  mutate(PC3 = abs(PC3)) %>%
  arrange(desc(PC3))
view(pc3_rotation_table)

pc4_rotation_table <- rotation_table %>%
  as.data.frame() %>%
  dplyr::select(PC4) %>%
  mutate(PC4 = abs(PC4)) %>%
  arrange(desc(PC4))
view(pc4_rotation_table)

pc5_rotation_table <- rotation_table %>%
  as.data.frame() %>%
  dplyr::select(PC5) %>%
  mutate(PC5 = abs(PC5)) %>%
  arrange(desc(PC5))
view(pc5_rotation_table)

pc6_rotation_table <- rotation_table %>%
  as.data.frame() %>%
  dplyr::select(PC6) %>%
  mutate(PC6 = abs(PC6)) %>%
  arrange(desc(PC6))
view(pc6_rotation_table)