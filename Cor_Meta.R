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


cor.MOFA <- psych::corr.test(factor_scores %>%
                               tibble::column_to_rownames('sample'),
                             master.table.numeric %>%
                               tibble::column_to_rownames('sample'),
                             method = "spearman",
                             adjust = "BH",
                             minlength = 2)

p_val_adj <- cor.MOFA$p.adj
cor_coef <- cor.MOFA$r

#plot p val
symmertic_breaks <- seq(from=0, to = 0.06, length.out = 256)
color=colorRampPalette(c('red4', 'white', 'blue4'))(256)
heatmap_p_val <- pheatmap::pheatmap(p_val_adj, main = "adjusted p-value",
                                    cluster_rows = T,
                                    cluster_cols = T,
                                    color = color,
                                    #display_numbers = cor_coef,
                                    breaks = symmertic_breaks)