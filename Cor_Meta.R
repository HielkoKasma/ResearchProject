meta_full <- metatable_ATL_repeated_subjects_with_dummies
library(tidyverse)
view(meta_full)
meta_clean <- meta_full[,-c(8,10,12,16,17)]
view(meta_clean)
library(dplyr)
meta_clean <- meta_clean %>%
  mutate(
    Season = ifelse(new_seasons == "summer_autumn", 1, 0),
  )
view(meta_clean)
meta_clean1 <- meta_clean[,c(-8, -29,-30)]
view(meta_clean[,c(29, 30)])
view(meta_clean1)
rm(meta_clean1)

#making the table of PCA values per sample for PC1:PC6
pca <- prcomp(t(RP_numeric), scale. = F)
pca_data <- data.frame(Sample=rownames(pca$x),
                       X=pca$x[,1], 
                       Y=pca$x[,2]) %>%
  mutate(group_id = str_remove(Sample, "_.*"))  %>%
  left_join(metatable, by = c("Sample" = "pca_ids"))

df_PCAvalues <- data.frame(pca_data$X) 
df_PCAvalues$PC2 <- pca_data$Y

pca <- prcomp(t(RP_numeric), scale. = F)
pca_data <- data.frame(Sample=rownames(pca$x),
                       X=pca$x[,3], 
                       Y=pca$x[,4]) %>%
  mutate(group_id = str_remove(Sample, "_.*"))  %>%
  left_join(metatable, by = c("Sample" = "pca_ids"))
df_PCAvalues$PC3 <- pca_data$X
df_PCAvalues$PC4 <- pca_data$Y

pca <- prcomp(t(RP_numeric), scale. = F)
pca_data <- data.frame(Sample=rownames(pca$x),
                       X=pca$x[,5], 
                       Y=pca$x[,6]) %>%
  mutate(group_id = str_remove(Sample, "_.*"))  %>%
  left_join(metatable, by = c("Sample" = "pca_ids"))
df_PCAvalues$PC5 <- pca_data$X
df_PCAvalues$PC6 <- pca_data$Y
rownames(df_PCAvalues) <- pca_data$Sample
colnames(df_PCAvalues)[1] <- "PC1"
view(df_PCAvalues)

#Cleaning up metadataframe
meta_clean1 <- meta_clean1 %>%
  mutate(
    Smoking_Status = ifelse(smoking.status == "Ex.smoker", 1, 0),
    Asthma_status = ifelse(asthma.status == "A", 1, 0),
    Sex = ifelse(gender == "female", 1, 0),
    Visit = ifelse(meth_visit == "Visit 1a", 1, 0),
    )
view(meta_clean1)
meta_clean1 <- meta_clean1[,-9]
view(meta_clean1)

meta_clean1_nohos <- meta_clean1[, -c(11:21)]
view(meta_clean1_nohos)

#Code Tatiana - Plotting heatmap

cor.MOFA <- psych::corr.test(df_PCAvalues,
                             
                             meta_clean1_nohos %>%
                               tibble::column_to_rownames('meth_file_id') %>%
                               dplyr::select(where(is.numeric)),
                             method = "spearman",
                             adjust = "BH",
                             minlength = 2)

p_val_adj <- cor.MOFA$p.adj
p_val_reg <- cor.MOFA$p
view(p_val_adj)
cor_coef <- cor.MOFA$r
view(cor_coef)

#plot p val
#install.packages("pheatmap")
library(pheatmap)
symmertic_breaks <- seq(from=0, to = 0.06, length.out = 256)
color=colorRampPalette(c('red4', 'white', 'blue4'))(256)

#heatmap for adjusted p value
desired_order <- c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6")
p_val_adj <- p_val_adj[desired_order, ]
heatmap_p_val <- pheatmap::pheatmap(p_val_adj, main = "adjusted p-value",
                                    cluster_rows = F,
                                    cluster_cols = T,
                                    color = color,
                                    #display_numbers = cor_coef,
                                    breaks = symmertic_breaks,
                                    fontsize_col = 6)

view(df_PCAvalues)

#heatmap for non-adjuested p value
heatmap_p_val <- pheatmap::pheatmap(p_val_reg, main = "non-adjuested p-value",
                                    cluster_rows = F,
                                    cluster_cols = T,
                                    color = color,
                                    #display_numbers = cor_coef,
                                    breaks = symmertic_breaks,
                                    fontsize_col = 6)


