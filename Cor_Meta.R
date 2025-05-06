meta_full <- metatable_ATL_repeated_subjects_with_dummies
library(tidyverse)
view(meta_full)
meta_clean <- meta_full[,-c(8,10,12,16,17)]
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

#The code below is for making the dataframe of PC values per sample
#change the values in teh pca_data code and df_PCAvalues code to get all values for PC1:PC6 and put them in their respective column in the dataframe
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

df_PCAvalues$PC6 <- pca_data$Y
view(df_PCAvalues)

#Code Tatiana

cor.MOFA <- psych::corr.test(df_PCAvalues,
                               
                             meta_clean1 %>%
                               tibble::column_to_rownames('meth_file_id') %>%
                             dplyr::select(where(is.numeric)),
                             method = "spearman",
                             adjust = "BH",
                             minlength = 2)

p_val_adj <- cor.MOFA$p.adj
view(p_val_adj)
cor_coef <- cor.MOFA$r
view(cor_coef)

#plot p val
install.packages("pheatmap")
library(pheatmap)
symmertic_breaks <- seq(from=0, to = 0.06, length.out = 256)
color=colorRampPalette(c('red4', 'white', 'blue4'))(256)
heatmap_p_val <- pheatmap::pheatmap(p_val_adj, main = "adjusted p-value",
                                    cluster_rows = T,
                                    cluster_cols = T,
                                    color = color,
                                    #display_numbers = cor_coef,
                                    breaks = symmertic_breaks,
                                    fontsize_col = 6)
meta_clean1 <- meta_clean1 %>%
  mutate(
    ex_smoker = ifelse(smoking.status == "Ex.smoker", 1, 0),
    non_smoker = ifelse(smoking.status == "Non.smoker", 1, 0),
    asthma = ifelse(asthma.status == "A", 1, 0),
    healthy = ifelse(asthma.status == "H", 1, 0),
    female = ifelse(gender == "female", 1, 0),
    male = ifelse(gender == "male", 1, 0)
  )
view(meta_clean1)
meta_clean1 <- meta_clean1[,-9]

View(meta_clean1)
#aaa
