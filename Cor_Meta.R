meta_full <- metatable_ATL_repeated_subjects_with_dummies
view(meta_full)
meta_clean <- meta_full[,-c(8,10,12,16,17,)]
view(meta_clean)
library(dplyr)
meta_clean <- meta_clean %>%
  mutate(
    summer_autumn = ifelse(new_seasons == "summer_autumn", 1, 0),
    winter_spring = ifelse(new_seasons == "winter_spring", 1, 0)
  )
view(meta_clean)
meta_clean1 <- meta_clean[,-8]
view(meta_clean1)

pca <- prcomp(t(RP_numeric), scale. = F)
pca_data <- data.frame(Sample=rownames(pca$x),
                       X=pca$x[,1], 
                       Y=pca$x[,2]) %>%
  mutate(group_id = str_remove(Sample, "_.*"))  %>%
  left_join(metatable, by = c("Sample" = "pca_ids"))
view(pca_data)

df_PCAvalues <- data.frame(pca_data$X) 
df_PCAvalues$PC2 <- pca_data$Y
df_PCAvalues
colnames(df_PCAvalues) <- c("PC1", "PC2")
df_PCAvalues
rownames(df_PCAvalues) <- pca_data$Sample
df_PCAvalues