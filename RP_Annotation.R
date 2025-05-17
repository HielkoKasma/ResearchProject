cpg_annotation <- annotation_cpgs_10perc_random
cpg_ann_rel_gene <- cpg_annotation_reletave_to_gene
coding_gene_ann <- coding_gene_annotation_grch38_101
view(cpg_ann_rel_gene)

rownames(cpg_ann_rel_gene) <- cpg_ann_rel_gene$cpg_ID
cpg_ann_rel_gene <- cpg_ann_rel_gene[,-1]
view(cpg_ann_rel_gene)

head(top1000_pc1_cpgs)
length(top1000_pc1_cpgs)
#cpg_ann_rel_gene is the one we need
str(top1000_pc1_cpgs)







#PC1 chisq
top1000_pc1cpg_ann <- cpg_ann_rel_gene[rownames(cpg_ann_rel_gene) %in% top1000_pc1_cpgs, ]
top1000_pc1_cpg_ann <- cpg_ann_rel_gene[top1000_pc1_cpgs, , drop = FALSE]
View(top1000_pc1_cpg_ann)

# how many cpgs are in,up and downstream the gene OVERALL
remaining_cpgs_pc1 <- cpg_ann_rel_gene[!(rownames(cpg_ann_rel_gene) %in% top1000_pc1_cpgs), ]
View(remaining_cpgs_pc1) #87,307 cpg sites now 
total_in_gene_pc1<-sum(grepl("^ENSG", remaining_cpgs_pc1$inside_gene)) #54395
total_up_gene_pc1<-sum(grepl("^ENSG", remaining_cpgs_pc1$upstream_gene)) #20146
total_down_gene_pc1<-sum(grepl("^ENSG", remaining_cpgs_pc1$downstream_gene)) #13233

# how many cpgs are in,up and downstream the gene PC1
in_gene_pc1<-sum(grepl("^ENSG", top1000_pc1_cpg_ann$inside_gene)) #562
up_gene_pc1<-sum(grepl("^ENSG", top1000_pc1_cpg_ann$upstream_gene)) #97
down_gene_pc1<-sum(grepl("^ENSG", top1000_pc1_cpg_ann$downstream_gene)) #65

#making the chisq table
total_cpgs <- nrow(remaining_cpgs_pc1)

# Inside gene
in_table_pc1 <- matrix(c(
  in_gene_pc1,
  1000 - in_gene_pc1,
  total_in_gene_pc1,
  total_cpgs - total_in_gene_pc1
), nrow = 2, byrow = TRUE)
colnames(in_table_pc1) <- c("In_gene", "Not_in_gene")
rownames(in_table_pc1) <- c("PC1", "Overall")
View(in_table_pc1)
chisq.test(in_table_pc1)$p.value #8.652664e-05

# Upstream
up_table_pc1 <- matrix(c(
  up_gene_pc1,
  1000 - up_gene_pc1,
  total_up_gene_pc1,
  total_cpgs - total_up_gene_pc1
), nrow = 2, byrow = TRUE)
colnames(up_table_pc1) <- c("Upstream", "Not_upstream")
rownames(up_table_pc1) <- c("PC1", "Overall")
View(up_table_pc1)
chisq.test(up_table_pc1)$p.value #2.122258e-23

# Downstream
down_table_pc1 <- matrix(c(
  down_gene_pc1,
  1000 - down_gene_pc1,
  total_down_gene_pc1,
  total_cpgs - total_down_gene_pc1
), nrow = 2, byrow = TRUE)
colnames(down_table_pc1) <- c("Downstream", "Not_downstream")
rownames(down_table_pc1) <- c("PC1", "Overall")
View(down_table_pc1)
chisq.test(down_table_pc1)$p.value #3.837709e-14










#PC2
annotated_cpg_top1000_pc2 <- cpg_ann_rel_gene[rownames(cpg_ann_rel_gene) 
                                              %in% top1000_pc2_cpgs, ]
view(annotated_cpg_top1000_pc2)



