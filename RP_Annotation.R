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
#for intergenic
total_inter_gene_pc1 <- nrow(remaining_cpgs_pc1[
  !grepl("^ENSG", remaining_cpgs_pc1$inside_gene) &
    !grepl("^ENSG", remaining_cpgs_pc1$upstream_gene) &
    !grepl("^ENSG", remaining_cpgs_pc1$downstream_gene), 
]) #17853


# how many cpgs are in,up and downstream the gene PC1
in_gene_pc1<-sum(grepl("^ENSG", top1000_pc1_cpg_ann$inside_gene)) #562
up_gene_pc1<-sum(grepl("^ENSG", top1000_pc1_cpg_ann$upstream_gene)) #97
down_gene_pc1<-sum(grepl("^ENSG", top1000_pc1_cpg_ann$downstream_gene)) #65
inter_gene_pc1 <- nrow(top1000_pc1_cpg_ann[
  !grepl("^ENSG", top1000_pc1_cpg_ann$inside_gene) &
    !grepl("^ENSG", top1000_pc1_cpg_ann$upstream_gene) &
    !grepl("^ENSG", top1000_pc1_cpg_ann$downstream_gene),
]) #336


#making the chisq table
total_cpgs1 <- nrow(remaining_cpgs_pc1)

# Inside gene
in_table_pc1 <- matrix(c(
  in_gene_pc1,
  1000 - in_gene_pc1,
  total_in_gene_pc1,
  total_cpgs1 - total_in_gene_pc1
), nrow = 2, byrow = TRUE)
colnames(in_table_pc1) <- c("In_gene", "Not_in_gene")
rownames(in_table_pc1) <- c("PC1", "Overall")
View(in_table_pc1)
chisq_in_gene_pc1<-chisq.test(in_table_pc1)$p.value #8.652664e-05

# Upstream
up_table_pc1 <- matrix(c(
  up_gene_pc1,
  1000 - up_gene_pc1,
  total_up_gene_pc1,
  total_cpgs1 - total_up_gene_pc1
), nrow = 2, byrow = TRUE)
colnames(up_table_pc1) <- c("Upstream", "Not_upstream")
rownames(up_table_pc1) <- c("PC1", "Overall")
View(up_table_pc1)
chisq_up_gene_pc1<-chisq.test(up_table_pc1)$p.value #2.122258e-23

# Downstream
down_table_pc1 <- matrix(c(
  down_gene_pc1,
  1000 - down_gene_pc1,
  total_down_gene_pc1,
  total_cpgs1 - total_down_gene_pc1
), nrow = 2, byrow = TRUE)
colnames(down_table_pc1) <- c("Downstream", "Not_downstream")
rownames(down_table_pc1) <- c("PC1", "Overall")
View(down_table_pc1)
chisq_down_gene_pc1<-chisq.test(down_table_pc1)$p.value #3.837709e-14

# Intergenic
intergene_table_pc1 <- matrix(c(
  inter_gene_pc1,
  1000 - inter_gene_pc1,
  total_inter_gene_pc1,
  total_cpgs1 - total_inter_gene_pc1
), nrow = 2, byrow = TRUE)
colnames(intergene_table_pc1) <- c("Intergenic", "Not_intergenic")
rownames(intergene_table_pc1) <- c("PC1", "Overall")
View(intergene_table_pc1)
chisq_intergene_pc1<-chisq.test(intergene_table_pc1)$p.value #2.289809e-24








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
total_inter_gene_pc2 <- nrow(remaining_cpgs_pc2[
  !grepl("^ENSG", remaining_cpgs_pc2$inside_gene) &
    !grepl("^ENSG", remaining_cpgs_pc2$upstream_gene) &
    !grepl("^ENSG", remaining_cpgs_pc2$downstream_gene), 
]) #18002

# how many cpgs are in,up and downstream the gene PC1
in_gene_pc2<-sum(grepl("^ENSG", top1000_pc2_cpg_ann$inside_gene)) #558
up_gene_pc2<-sum(grepl("^ENSG", top1000_pc2_cpg_ann$upstream_gene)) #311
down_gene_pc2<-sum(grepl("^ENSG", top1000_pc2_cpg_ann$downstream_gene)) #125
inter_gene_pc2 <- nrow(top1000_pc2_cpg_ann[
  !grepl("^ENSG", top1000_pc2_cpg_ann$inside_gene) &
    !grepl("^ENSG", top1000_pc2_cpg_ann$upstream_gene) &
    !grepl("^ENSG", top1000_pc2_cpg_ann$downstream_gene),
]) #187

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
chisq_in_gene_pc2<-chisq.test(in_table_pc2)$p.value #2.814641e-05

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
chisq_up_gene_pc2<-chisq.test(up_table_pc2)$p.value #7.81515e-10

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
chisq_down_gene_pc2<-chisq.test(down_table_pc2)$p.value #0.02568461

# Intergenic
intergene_table_pc2 <- matrix(c(
  inter_gene_pc2,
  1000 - inter_gene_pc2,
  total_inter_gene_pc2,
  total_cpgs2 - total_inter_gene_pc2
), nrow = 2, byrow = TRUE)
colnames(intergene_table_pc2) <- c("Intergenic", "Not_intergenic")
rownames(intergene_table_pc2) <- c("PC1", "Overall")
View(intergene_table_pc2)
chisq_intergene_pc2<-chisq.test(intergene_table_pc2)$p.value #0.1462621








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
total_inter_gene_pc3 <- nrow(remaining_cpgs_pc3[
  !grepl("^ENSG", remaining_cpgs_pc3$inside_gene) &
    !grepl("^ENSG", remaining_cpgs_pc3$upstream_gene) &
    !grepl("^ENSG", remaining_cpgs_pc3$downstream_gene), 
]) #17852

# how many cpgs are in,up and downstream the gene PC1
in_gene_pc3<-sum(grepl("^ENSG", top1000_pc3_cpg_ann$inside_gene)) #539
up_gene_pc3<-sum(grepl("^ENSG", top1000_pc3_cpg_ann$upstream_gene)) #114
down_gene_pc3<-sum(grepl("^ENSG", top1000_pc3_cpg_ann$downstream_gene)) #91
inter_gene_pc3 <- nrow(top1000_pc3_cpg_ann[
  !grepl("^ENSG", top1000_pc3_cpg_ann$inside_gene) &
    !grepl("^ENSG", top1000_pc3_cpg_ann$upstream_gene) &
    !grepl("^ENSG", top1000_pc3_cpg_ann$downstream_gene),
]) #337

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
chisq_in_gene_pc3<-chisq.test(in_table_pc3)$p.value #5.498688e-08

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
chisq_up_gene_pc3<-chisq.test(up_table_pc3)$p.value #3.925944e-18

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
chisq_down_gene_pc3<-chisq.test(down_table_pc3)$p.value #1.485411e-07

# Intergenic
intergene_table_pc3 <- matrix(c(
  inter_gene_pc3,
  1000 - inter_gene_pc3,
  total_inter_gene_pc3,
  total_cpgs3 - total_inter_gene_pc3
), nrow = 2, byrow = TRUE)
colnames(intergene_table_pc3) <- c("Intergenic", "Not_intergenic")
rownames(intergene_table_pc3) <- c("PC1", "Overall")
View(intergene_table_pc3)
chisq_intergene_pc3<-chisq.test(intergene_table_pc3)$p.value #1.01693e-24








#PC4
head(top1000_pc4_cpgs)
top1000_pc4cpg_ann <- cpg_ann_rel_gene[rownames(cpg_ann_rel_gene) %in% top1000_pc4_cpgs, ]
top1000_pc4_cpg_ann <- cpg_ann_rel_gene[top1000_pc4_cpgs, , drop = FALSE]
View(top1000_pc4_cpg_ann)

# how many cpgs are in,up and downstream the gene OVERALL
remaining_cpgs_pc4 <- cpg_ann_rel_gene[!(rownames(cpg_ann_rel_gene) %in% top1000_pc4_cpgs), ]
View(remaining_cpgs_pc4) #87,307 cpg sites now 
total_in_gene_pc4<-sum(grepl("^ENSG", remaining_cpgs_pc4$inside_gene)) #54392
total_up_gene_pc4<-sum(grepl("^ENSG", remaining_cpgs_pc4$upstream_gene)) # 20108
total_down_gene_pc4<-sum(grepl("^ENSG", remaining_cpgs_pc4$downstream_gene)) #13203
total_inter_gene_pc4 <- nrow(remaining_cpgs_pc4[
  !grepl("^ENSG", remaining_cpgs_pc4$inside_gene) &
    !grepl("^ENSG", remaining_cpgs_pc4$upstream_gene) &
    !grepl("^ENSG", remaining_cpgs_pc4$downstream_gene), 
]) #17891

# how many cpgs are in,up and downstream the gene PC1
in_gene_pc4<-sum(grepl("^ENSG", top1000_pc4_cpg_ann$inside_gene)) #565
up_gene_pc4<-sum(grepl("^ENSG", top1000_pc4_cpg_ann$upstream_gene)) #135
down_gene_pc4<-sum(grepl("^ENSG", top1000_pc4_cpg_ann$downstream_gene)) #95
inter_gene_pc4 <- nrow(top1000_pc4_cpg_ann[
  !grepl("^ENSG", top1000_pc4_cpg_ann$inside_gene) &
    !grepl("^ENSG", top1000_pc4_cpg_ann$upstream_gene) &
    !grepl("^ENSG", top1000_pc4_cpg_ann$downstream_gene),
]) #298

#making the chisq table
total_cpgs4 <- nrow(remaining_cpgs_pc4)

# Inside gene
in_table_pc4 <- matrix(c(
  in_gene_pc4,
  1000 - in_gene_pc4,
  total_in_gene_pc4,
  total_cpgs4 - total_in_gene_pc4
), nrow = 2, byrow = TRUE)
colnames(in_table_pc4) <- c("In_gene", "Not_in_gene")
rownames(in_table_pc4) <- c("PC1", "Overall")
View(in_table_pc4)
chisq_in_gene_pc4<-chisq.test(in_table_pc4)$p.value #0.0001924227

# Upstream
up_table_pc4 <- matrix(c(
  up_gene_pc4,
  1000 - up_gene_pc4,
  total_up_gene_pc4,
  total_cpgs4 - total_up_gene_pc4
), nrow = 2, byrow = TRUE)
colnames(up_table_pc4) <- c("Upstream", "Not_upstream")
rownames(up_table_pc4) <- c("PC1", "Overall")
View(up_table_pc4)
chisq_up_gene_pc4<-chisq.test(up_table_pc4)$p.value #1.321539e-12

# Downstream
down_table_pc4 <- matrix(c(
  down_gene_pc4,
  1000 - down_gene_pc4,
  total_down_gene_pc4,
  total_cpgs4 - total_down_gene_pc4
), nrow = 2, byrow = TRUE)
colnames(down_table_pc4) <- c("Downstream", "Not_downstream")
rownames(down_table_pc4) <- c("PC1", "Overall")
View(down_table_pc4)
chisq_down_gene_pc4<-chisq.test(down_table_pc4)$p.value #9.649082e-07

# Intergenic
intergene_table_pc4 <- matrix(c(
  inter_gene_pc4,
  1000 - inter_gene_pc4,
  total_inter_gene_pc4,
  total_cpgs4 - total_inter_gene_pc4
), nrow = 2, byrow = TRUE)
colnames(intergene_table_pc4) <- c("Intergenic", "Not_intergenic")
rownames(intergene_table_pc4) <- c("PC1", "Overall")
View(intergene_table_pc4)
chisq_intergene_pc4<-chisq.test(intergene_table_pc4)$p.value #6.126026e-13










#PC5
head(top1000_pc5_cpgs)
top1000_pc5cpg_ann <- cpg_ann_rel_gene[rownames(cpg_ann_rel_gene) %in% top1000_pc5_cpgs, ]
top1000_pc5_cpg_ann <- cpg_ann_rel_gene[top1000_pc5_cpgs, , drop = FALSE]
View(top1000_pc5_cpg_ann)

# how many cpgs are in,up and downstream the gene OVERALL
remaining_cpgs_pc5 <- cpg_ann_rel_gene[!(rownames(cpg_ann_rel_gene) %in% top1000_pc5_cpgs), ]
View(remaining_cpgs_pc5) #87,307 cpg sites now 
total_in_gene_pc5<-sum(grepl("^ENSG", remaining_cpgs_pc5$inside_gene)) #54323
total_up_gene_pc5<-sum(grepl("^ENSG", remaining_cpgs_pc5$upstream_gene)) # 20074
total_down_gene_pc5<-sum(grepl("^ENSG", remaining_cpgs_pc5$downstream_gene)) #13183
total_inter_gene_pc5 <- nrow(remaining_cpgs_pc5[
  !grepl("^ENSG", remaining_cpgs_pc5$inside_gene) &
    !grepl("^ENSG", remaining_cpgs_pc5$upstream_gene) &
    !grepl("^ENSG", remaining_cpgs_pc5$downstream_gene), 
]) #17966

# how many cpgs are in,up and downstream the gene PC1
in_gene_pc5<-sum(grepl("^ENSG", top1000_pc5_cpg_ann$inside_gene)) #634
up_gene_pc5<-sum(grepl("^ENSG", top1000_pc5_cpg_ann$upstream_gene)) #169
down_gene_pc5<-sum(grepl("^ENSG", top1000_pc5_cpg_ann$downstream_gene)) #115
inter_gene_pc5 <- nrow(top1000_pc5_cpg_ann[
  !grepl("^ENSG", top1000_pc5_cpg_ann$inside_gene) &
    !grepl("^ENSG", top1000_pc5_cpg_ann$upstream_gene) &
    !grepl("^ENSG", top1000_pc5_cpg_ann$downstream_gene),
]) #223

#making the chisq table
total_cpgs5 <- nrow(remaining_cpgs_pc5)

# Inside gene
in_table_pc5 <- matrix(c(
  in_gene_pc5,
  1000 - in_gene_pc5,
  total_in_gene_pc5,
  total_cpgs5 - total_in_gene_pc5
), nrow = 2, byrow = TRUE)
colnames(in_table_pc5) <- c("In_gene", "Not_in_gene")
rownames(in_table_pc5) <- c("PC1", "Overall")
View(in_table_pc5)
chisq_in_gene_pc5<-chisq.test(in_table_pc5)$p.value #0.4641147 insignificant

# Upstream
up_table_pc5 <- matrix(c(
  up_gene_pc5,
  1000 - up_gene_pc5,
  total_up_gene_pc5,
  total_cpgs5 - total_up_gene_pc5
), nrow = 2, byrow = TRUE)
colnames(up_table_pc5) <- c("Upstream", "Not_upstream")
rownames(up_table_pc5) <- c("PC1", "Overall")
View(up_table_pc5)
chisq_up_gene_pc5<-chisq.test(up_table_pc5)$p.value #6.196779e-06

# Downstream
down_table_pc5 <- matrix(c(
  down_gene_pc5,
  1000 - down_gene_pc5,
  total_down_gene_pc5,
  total_cpgs5 - total_down_gene_pc5
), nrow = 2, byrow = TRUE)
colnames(down_table_pc5) <- c("Downstream", "Not_downstream")
rownames(down_table_pc5) <- c("PC1", "Overall")
View(down_table_pc5)
chisq_down_gene_pc5<-chisq.test(down_table_pc5)$p.value #0.001807368

# Intergenic
intergene_table_pc5 <- matrix(c(
  inter_gene_pc5,
  1000 - inter_gene_pc5,
  total_inter_gene_pc5,
  total_cpgs5 - total_inter_gene_pc5
), nrow = 2, byrow = TRUE)
colnames(intergene_table_pc5) <- c("Intergenic", "Not_intergenic")
rownames(intergene_table_pc5) <- c("PC1", "Overall")
View(intergene_table_pc5)
chisq_intergene_pc5<-chisq.test(intergene_table_pc5)$p.value #0.1937474












#PC6
head(top1000_pc6_cpgs)
top1000_pc6cpg_ann <- cpg_ann_rel_gene[rownames(cpg_ann_rel_gene) %in% top1000_pc6_cpgs, ]
top1000_pc6_cpg_ann <- cpg_ann_rel_gene[top1000_pc6_cpgs, , drop = FALSE]
View(top1000_pc6_cpg_ann)

# how many cpgs are in,up and downstream the gene OVERALL
remaining_cpgs_pc6 <- cpg_ann_rel_gene[!(rownames(cpg_ann_rel_gene) %in% top1000_pc6_cpgs), ]
View(remaining_cpgs_pc6) #87,307 cpg sites now 
total_in_gene_pc6<-sum(grepl("^ENSG", remaining_cpgs_pc6$inside_gene)) #54399
total_up_gene_pc6<-sum(grepl("^ENSG", remaining_cpgs_pc6$upstream_gene)) # 20085
total_down_gene_pc6<-sum(grepl("^ENSG", remaining_cpgs_pc6$downstream_gene)) #13182
total_inter_gene_pc6 <- nrow(remaining_cpgs_pc6[
  !grepl("^ENSG", remaining_cpgs_pc6$inside_gene) &
    !grepl("^ENSG", remaining_cpgs_pc6$upstream_gene) &
    !grepl("^ENSG", remaining_cpgs_pc6$downstream_gene), 
]) #17900

# how many cpgs are in,up and downstream the gene PC1
in_gene_pc6<-sum(grepl("^ENSG", top1000_pc6_cpg_ann$inside_gene)) #558
up_gene_pc6<-sum(grepl("^ENSG", top1000_pc6_cpg_ann$upstream_gene)) #158
down_gene_pc6<-sum(grepl("^ENSG", top1000_pc6_cpg_ann$downstream_gene)) #116
inter_gene_pc6 <- nrow(top1000_pc6_cpg_ann[
  !grepl("^ENSG", top1000_pc6_cpg_ann$inside_gene) &
    !grepl("^ENSG", top1000_pc6_cpg_ann$upstream_gene) &
    !grepl("^ENSG", top1000_pc6_cpg_ann$downstream_gene),
]) #289

#making the chisq table
total_cpgs6 <- nrow(remaining_cpgs_pc6)

# Inside gene
in_table_pc6 <- matrix(c(
  in_gene_pc6,
  1000 - in_gene_pc6,
  total_in_gene_pc6,
  total_cpgs6 - total_in_gene_pc6
), nrow = 2, byrow = TRUE)
colnames(in_table_pc6) <- c("In_gene", "Not_in_gene")
rownames(in_table_pc6) <- c("PC1", "Overall")
View(in_table_pc6)
chisq_in_gene_pc6<-chisq.test(in_table_pc6)$p.value #2.814641e-05

# Upstream
up_table_pc6 <- matrix(c(
  up_gene_pc6,
  1000 - up_gene_pc6,
  total_up_gene_pc6,
  total_cpgs6 - total_up_gene_pc6
), nrow = 2, byrow = TRUE)
colnames(up_table_pc6) <- c("Upstream", "Not_upstream")
rownames(up_table_pc6) <- c("PC1", "Overall")
View(up_table_pc6)
chisq_up_gene_pc6<-chisq.test(up_table_pc6)$p.value #8.706796e-08

# Downstream
down_table_pc6 <- matrix(c(
  down_gene_pc6,
  1000 - down_gene_pc6,
  total_down_gene_pc6,
  total_cpgs6 - total_down_gene_pc6
), nrow = 2, byrow = TRUE)
colnames(down_table_pc6) <- c("Downstream", "Not_downstream")
rownames(down_table_pc6) <- c("PC1", "Overall")
View(down_table_pc6)
chisq_down_gene_pc6<-chisq.test(down_table_pc6)$p.value #0.002435336

# Intergenic
intergene_table_pc6 <- matrix(c(
  inter_gene_pc6,
  1000 - inter_gene_pc6,
  total_inter_gene_pc6,
  total_cpgs6 - total_inter_gene_pc6
), nrow = 2, byrow = TRUE)
colnames(intergene_table_pc6) <- c("Intergenic", "Not_intergenic")
rownames(intergene_table_pc6) <- c("PC1", "Overall")
View(intergene_table_pc6)
chisq_intergene_pc6<-chisq.test(intergene_table_pc6)$p.value #8.590559e-11





#making the chisq p-values table: 
chisq_pval_matrix <- matrix(c(
  chisq_in_gene_pc1, chisq_up_gene_pc1, chisq_down_gene_pc1, chisq_intergene_pc1,
  chisq_in_gene_pc2, chisq_up_gene_pc2, chisq_down_gene_pc2, chisq_intergene_pc2,
  chisq_in_gene_pc3, chisq_up_gene_pc3, chisq_down_gene_pc3, chisq_intergene_pc3,
  chisq_in_gene_pc4, chisq_up_gene_pc4, chisq_down_gene_pc4, chisq_intergene_pc4,
  chisq_in_gene_pc5, chisq_up_gene_pc5, chisq_down_gene_pc5, chisq_intergene_pc5,
  chisq_in_gene_pc6, chisq_up_gene_pc6, chisq_down_gene_pc6, chisq_intergene_pc6
), 
nrow = 4, ncol = 6, byrow = FALSE)

# Set row and column names
rownames(chisq_pval_matrix) <- c("In_gene", "Upstream", "Downstream", "Intergenic")
colnames(chisq_pval_matrix) <- c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6")

# View the matrix 
View(chisq_pval_matrix)


