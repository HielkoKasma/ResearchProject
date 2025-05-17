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
total_cpgs1 <- nrow(remaining_cpgs_pc1)

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
head(top1000_pc2_cpgs)
top1000_pc2cpg_ann <- cpg_ann_rel_gene[rownames(cpg_ann_rel_gene) %in% top1000_pc2_cpgs, ]
top1000_pc2_cpg_ann <- cpg_ann_rel_gene[top1000_pc2_cpgs, , drop = FALSE]
View(top1000_pc2_cpg_ann)

# how many cpgs are in,up and downstream the gene OVERALL
remaining_cpgs_pc2 <- cpg_ann_rel_gene[!(rownames(cpg_ann_rel_gene) %in% top1000_pc2_cpgs), ]
View(remaining_cpgs_pc2) #87,307 cpg sites now 
total_in_gene_pc2<-sum(grepl("^ENSG", remaining_cpgs_pc2$inside_gene)) #54399
total_up_gene_pc2<-sum(grepl("^ENSG", remaining_cpgs_pc2$upstream_gene)) # 19932
total_down_gene_pc2<-sum(grepl("^ENSG", remaining_cpgs_pc2$downstream_gene)) #13173

# how many cpgs are in,up and downstream the gene PC1
in_gene_pc2<-sum(grepl("^ENSG", top1000_pc2_cpg_ann$inside_gene)) #558
up_gene_pc2<-sum(grepl("^ENSG", top1000_pc2_cpg_ann$upstream_gene)) #311
down_gene_pc2<-sum(grepl("^ENSG", top1000_pc2_cpg_ann$downstream_gene)) #125

#making the chisq table
total_cpgs2 <- nrow(remaining_cpgs_pc2)

# Inside gene
in_table_pc2 <- matrix(c(
  in_gene_pc2,
  1000 - in_gene_pc2,
  total_in_gene_pc2,
  total_cpgs2 - total_in_gene_pc2
), nrow = 2, byrow = TRUE)
colnames(in_table_pc2) <- c("In_gene", "Not_in_gene")
rownames(in_table_pc2) <- c("PC1", "Overall")
View(in_table_pc2)
chisq.test(in_table_pc2)$p.value #2.814641e-05

# Upstream
up_table_pc2 <- matrix(c(
  up_gene_pc2,
  1000 - up_gene_pc2,
  total_up_gene_pc2,
  total_cpgs2 - total_up_gene_pc2
), nrow = 2, byrow = TRUE)
colnames(up_table_pc2) <- c("Upstream", "Not_upstream")
rownames(up_table_pc2) <- c("PC1", "Overall")
View(up_table_pc2)
chisq.test(up_table_pc2)$p.value #7.81515e-10

# Downstream
down_table_pc2 <- matrix(c(
  down_gene_pc2,
  1000 - down_gene_pc2,
  total_down_gene_pc2,
  total_cpgs2 - total_down_gene_pc2
), nrow = 2, byrow = TRUE)
colnames(down_table_pc2) <- c("Downstream", "Not_downstream")
rownames(down_table_pc2) <- c("PC1", "Overall")
View(down_table_pc2)
chisq.test(down_table_pc2)$p.value #0.02568461









#PC3
head(top1000_pc3_cpgs)
top1000_pc3cpg_ann <- cpg_ann_rel_gene[rownames(cpg_ann_rel_gene) %in% top1000_pc3_cpgs, ]
top1000_pc3_cpg_ann <- cpg_ann_rel_gene[top1000_pc3_cpgs, , drop = FALSE]
View(top1000_pc3_cpg_ann)

# how many cpgs are in,up and downstream the gene OVERALL
remaining_cpgs_pc3 <- cpg_ann_rel_gene[!(rownames(cpg_ann_rel_gene) %in% top1000_pc3_cpgs), ]
View(remaining_cpgs_pc3) #87,307 cpg sites now 
total_in_gene_pc3<-sum(grepl("^ENSG", remaining_cpgs_pc3$inside_gene)) #54418
total_up_gene_pc3<-sum(grepl("^ENSG", remaining_cpgs_pc3$upstream_gene)) # 20129
total_down_gene_pc3<-sum(grepl("^ENSG", remaining_cpgs_pc3$downstream_gene)) #13207

# how many cpgs are in,up and downstream the gene PC1
in_gene_pc3<-sum(grepl("^ENSG", top1000_pc3_cpg_ann$inside_gene)) #539
up_gene_pc3<-sum(grepl("^ENSG", top1000_pc3_cpg_ann$upstream_gene)) #114
down_gene_pc3<-sum(grepl("^ENSG", top1000_pc3_cpg_ann$downstream_gene)) #91

#making the chisq table
total_cpgs3 <- nrow(remaining_cpgs_pc3)

# Inside gene
in_table_pc3 <- matrix(c(
  in_gene_pc3,
  1000 - in_gene_pc3,
  total_in_gene_pc3,
  total_cpgs3 - total_in_gene_pc3
), nrow = 2, byrow = TRUE)
colnames(in_table_pc3) <- c("In_gene", "Not_in_gene")
rownames(in_table_pc3) <- c("PC1", "Overall")
View(in_table_pc3)
chisq.test(in_table_pc3)$p.value #5.498688e-08

# Upstream
up_table_pc3 <- matrix(c(
  up_gene_pc3,
  1000 - up_gene_pc3,
  total_up_gene_pc3,
  total_cpgs3 - total_up_gene_pc3
), nrow = 2, byrow = TRUE)
colnames(up_table_pc3) <- c("Upstream", "Not_upstream")
rownames(up_table_pc3) <- c("PC1", "Overall")
View(up_table_pc3)
chisq.test(up_table_pc3)$p.value #3.925944e-18

# Downstream
down_table_pc3 <- matrix(c(
  down_gene_pc3,
  1000 - down_gene_pc3,
  total_down_gene_pc3,
  total_cpgs3 - total_down_gene_pc3
), nrow = 2, byrow = TRUE)
colnames(down_table_pc3) <- c("Downstream", "Not_downstream")
rownames(down_table_pc3) <- c("PC1", "Overall")
View(down_table_pc3)
chisq.test(down_table_pc3)$p.value #1.485411e-07
