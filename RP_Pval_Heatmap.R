#Had problems with the chisq matrix, so we worked around it
load("C:/Users/hielk/Downloads/chisq_pval_matrix.RData")
view(chisq_pval_matrix)

chisq_pval_matrix
log_chisq_pval_matrix <- -log10(chisq_pval_matrix)
view(log_chisq_pval_matrix)
max(log_chisq_pval_matrix) #to find upper bound for heatmap

symmertic_breaks <- seq(from=0, to = 25, length.out = 256)
color=colorRampPalette(c('blue4', 'white', 'red4'))(256)
heatmap_p_val <- pheatmap::pheatmap(log_chisq_pval_matrix, main = "-log10 P-values Chisq Test Per PC Non Adjusted",
                                    cluster_rows = F,
                                    cluster_cols = F,
                                    color = color,
                                    #display_numbers = cor_coef,
                                    breaks = symmertic_breaks,
                                    fontsize_col = 10,
                                    fontsize_row = 10 )


pvals_vector <- as.vector(chisq_pval_matrix)
pvals_adj <- p.adjust(pvals_vector, method = "BH")
pval_adj_matrix <- matrix(pvals_adj, 
                          nrow = nrow(chisq_pval_matrix), 
                          ncol = ncol(chisq_pval_matrix),
                          dimnames = dimnames(chisq_pval_matrix))
view(pval_adj_matrix)

log_chisq_adjpval_matrix <- -log10(pval_adj_matrix)
view(log_chisq_adjpval_matrix)

heatmap_p_val <- pheatmap::pheatmap(log_chisq_adjpval_matrix, main = "-log10 P-values Chisq Test Per PC Adjusted (BH)",
                                    cluster_rows = F,
                                    cluster_cols = F,
                                    color = color,
                                    #display_numbers = cor_coef,
                                    breaks = symmertic_breaks,
                                    fontsize_col = 10,
                                    fontsize_row = 10 )
