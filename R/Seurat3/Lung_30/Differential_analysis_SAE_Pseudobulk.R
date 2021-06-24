########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################
invisible(lapply(c("Seurat","dplyr","cowplot","pbapply",
                   "magrittr","data.table","future","ggplot2","tidyr"), function(x) {
                           suppressPackageStartupMessages(library(x,character.only = T))
                   }))

# ######################################################################
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat3_functions.R")
save.path = "Yang/Lung_30/DE_analysis/surface_airway_epithelial/"
#======1.2 load  Seurat =========================
# load files
object = readRDS(file = "data/Lung_SCT_30_20200710.rds")
DefaultAssay(object) = "SCT"
Idents(object) = "Doublets"
object <- subset(object, idents = "Singlet")

object$cell_types %<>% gsub("d-S","TASC",.)
SAE <- c("BC1","BC2","BC-p","IC1","IC2","IC3","S","TASC","H","p-C","C1","C2","C3","Ion","NE")
Idents(object) = "cell_types"
object %<>% subset(idents = SAE)

# Option A:
Idents(object) = "orig.ident"
exp <- AverageExpression(object, assays = "SCT",slot = "data")
exp = exp$SCT
group_A = t(data.frame(c("UNC-48-P", "UNC-48-D"),
                        c("UNC-55-P", "UNC-55-D"),
                        c("UNC-66-P", "UNC-66-D"),
                        c("UNC-69-P", "UNC-69-D"),
                        c("UNC-71-P", "UNC-71-D")))
rownames(group_A) = gsub("-P","",group_A[,1])
Rowsum = rowSums(exp[,c(group_A[,1],group_A[,2])])
exp = exp[Rowsum >0,]

FindMarkers_bulk <- function(data, ident.1 = NULL,ident.2 = NULL){
        p_val <- pbsapply(rownames(data),function(gene){
                res = wilcox.test(x = data[gene,ident.1], y = data[gene,ident.2],
                                  paired = TRUE, alternative = "two.sided")
                res$p.value
        })
        ident1_pos = data[,ident.1] > 0
        pct_1 = rowSums(ident1_pos)/ncol(ident1_pos)
        ident2_pos = data[,ident.2] > 0
        pct_2 = rowSums(ident2_pos)/ncol(ident2_pos)
        
        avg_UMI_1 = log2(rowSums(expm1(data[,ident.1]))+1)
        avg_UMI_2 = log2(rowSums(expm1(data[,ident.2]))+1)
        avg_logFC = avg_UMI_1-avg_UMI_2
        p_adj = p.adjust(p_val, method = "BH")
        
        markers = data.frame('p_val'=p_val,
                             'avg_log2FC'=avg_logFC,
                             'pct.1'=pct_1,
                             'pct.2'=pct_2,
                             'p_val_adj'=p_adj,
                             'avg_log2_UMI.1'=avg_UMI_1,
                             'avg_log2_UMI.2'=avg_UMI_2,
                             'gene'=rownames(data))
        markers = markers[order(markers$avg_log2FC,decreasing = T),]
        return(markers)
}

markers_A = FindMarkers_bulk(exp,group_A[,1],group_A[,2])
markers_A$'cluster' = "P vs. D"
write.csv(markers_A,file = paste0(save.path,"Option_A_SAE_bulk_DEGs_P_vs_D.csv"))

# Option B:
object$orig.ident %<>% plyr::mapvalues(from = c("UNC-48-D","UNC-48-T",
                                                "UNC-55-D","UNC-55-T",
                                                "UNC-66-D","UNC-66-T",
                                                "UNC-69-D","UNC-69-T",
                                                "UNC-71-D","UNC-71-T"),
                                       to = c("UNC-48-D+T","UNC-48-D+T",
                                              "UNC-55-D+T","UNC-55-D+T",
                                              "UNC-66-D+T","UNC-66-D+T",
                                              "UNC-69-D+T","UNC-69-D+T",
                                              "UNC-71-D+T","UNC-71-D+T"))
Idents(object) = "orig.ident"
exp <- AverageExpression(object, assays = "SCT",slot = "data")
exp = exp$SCT
group_B = t(data.frame(c("UNC-48-P", "UNC-48-D+T"),
                       c("UNC-55-P", "UNC-55-D+T"),
                       c("UNC-66-P", "UNC-66-D+T"),
                       c("UNC-69-P", "UNC-69-D+T"),
                       c("UNC-71-P", "UNC-71-D+T")))
rownames(group_B) = gsub("-P","",group_B[,1])

Rowsum = rowSums(exp[,c(group_B[,1],group_B[,2])])
exp = exp[Rowsum >0,]

markers_B = FindMarkers_bulk(exp,group_B[,1],group_B[,2])
markers_B$'cluster' = "P vs. D+T"

write.csv(markers_B,file = paste0(save.path,"Option_B_SAE_Pseudobulk_DEGs_P_vs_D+T.csv"))
