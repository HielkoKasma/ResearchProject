cpg_annotation <- annotation_cpgs_10perc_random
cpg_ann_rel_gene <- cpg_annotation_reletave_to_gene
coding_gene_ann <- coding_gene_annotation_grch38_101

rownames(cpg_ann_rel_gene) <- cpg_ann_rel_gene$cpg_ID
cpg_ann_rel_gene <- cpg_ann_rel_gene[,-1]
view(cpg_ann_rel_gene)

head(top100_pc1_cpgs)
length(top1000_pc1_cpgs)
#cpg_ann_rel_gene is the one we need
str(top1000_pc1_cpgs)

annotated_cpg_top1000_pc1 <- cpg_ann_rel_gene[rownames(cpg_ann_rel_gene) %in% top1000_pc1_cpgs, ]

nrow(annotated_cpg_top1000_pc1)

annotated_cpg_top1000_pc2 <- cpg_ann_rel_gene[rownames(cpg_ann_rel_gene) 
                                              %in% top1000_pc2_cpgs, ]
view(annotated_cpg_top1000_pc2)



