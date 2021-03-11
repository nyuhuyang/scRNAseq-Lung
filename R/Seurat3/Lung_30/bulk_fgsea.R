library(dplyr)
library(magrittr)
library(stringr)
library(ggsci)
library(fgsea)
library(openxlsx)
library(progress)
library(BiocParallel)
library(tidyr)
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat3_functions.R")

path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)

TPM = read.csv(file = "Yang/RNA-seq/TPM_log2FC.csv",row.names =1)
#===========================
save.path = paste0("Yang/RNA-seq/fgsea/")
if(!dir.exists(save.path))dir.create(save.path, recursive = T)

TPM$gene = rownames(TPM)
TPM$up_regulated = ""
TPM[(TPM$Mean_log2FC_D_vs_P > 0 & TPM$ttest < 0.05),"up_regulated"] = "green"
TPM[(TPM$Mean_log2FC_D_vs_P < 0 & TPM$ttest < 0.05),"up_regulated"] = "blue"

dataset1 <- TPM[TPM$up_regulated %in% c("green", "blue"),
                c("gene","Mean_log2FC_D_vs_P","up_regulated")] 
colnames(dataset1) = c("gene","avg_logFC","clusters")

sheets = c("BC1","BC2","BC-p","IC1","IC2","IC3","S","d-S",
           "H","p-C","C1","C2","C3","Ion","NE")
read.path = "Yang/Lung_30/DE_analysis/groups/DE_results_Surface Airway Epithelial.xlsx"
DEG_list1 <- lapply(sheets, function(x) {
        tmp = readxl::read_excel(read.path,sheet = x)
        gene = pull(tmp[tmp$p_val_adj < 0.05 & tmp$avg_logFC >= 0.5,"gene"])
        return(gene)
})
DEG_list2 <- lapply(sheets, function(x) {
        tmp = readxl::read_excel(read.path,sheet = x)
        gene = pull(tmp[tmp$p_val_adj < 0.05 & tmp$avg_logFC >= 1,"gene"])
        return(gene)
})
names(DEG_list1) = sheets
names(DEG_list2) = sheets
df_list <- list()
dataset1$avg_logFC %<>% abs()
df_list[[1]] <- FgseaDotPlot(stats = dataset1, pathways = DEG_list1,
             Colv = F, pval = 0.05, padj = 0.25,size = " -log10(pval)",
             font.xtickslab = 10,
             font.ytickslab = 10,
             order.yaxis = sheets,
             cols = c("blue","green","yellow","orange","chocolate1","red"),
             rotate = F,
             rm.na = F,
             save.path = save.path, file.name = "GSEA_strategy_1_option_1",
             width = 2.5, height= 4, do.return = T)

df_list[[2]] <- FgseaDotPlot(stats = dataset1, pathways = DEG_list2,
             Colv = F, pval = 0.05, padj = 0.25,size = " -log10(pval)",
             font.xtickslab = 10,
             font.ytickslab = 10,
             order.yaxis = sheets,
             cols = c("blue","green","yellow","orange","chocolate1","red"),
             rotate = F,
             rm.na = T,
             save.path = save.path, file.name = "GSEA_strategy_1_option_2",
             width = 2.5, height= 3, do.return = T)
dataset2 = TPM[rownames(dataset1),grep("log2FC_D_vs_P_UNC",colnames(TPM))]
colnames(dataset2) = paste0("log2FC",1:6)
dataset2 %<>% tibble::rownames_to_column(var = "gene")
dataset2 %<>% gather(key = "clusters", value = "avg_logFC",-gene)

df_list[[3]] <- FgseaDotPlot(stats = dataset2, pathways = DEG_list1,
             Colv = F, pval = 0.05, padj = 0.25,size = " -log10(pval)",
             font.xtickslab = 10,
             font.ytickslab = 10,
             order.yaxis = sheets,
             cols = c("blue","green","yellow","orange","chocolate1","red"),
             rotate = T,save.path = save.path, file.name = "GSEA_strategy_2_option_1",
             width = 5, height= 4, do.return = T)

df_list[[4]] <- FgseaDotPlot(stats = dataset2, pathways = DEG_list2,
             Colv = F, pval = 0.05, padj = 0.25,size = " -log10(pval)",
             font.xtickslab = 10,
             font.ytickslab = 10,
             order.yaxis = sheets,
             cols = c("blue","green","yellow","orange","chocolate1","red"),
             rotate = T,save.path = save.path, file.name = "GSEA_strategy_2_option_2",
             width = 5, height= 4, do.return = T)
#df_fgseaRes = df_fgseaRes[df_fgseaRes$padj < 0.05, ]

names(df_list) = c("strategy_1_option_1","strategy_1_option_2",
                   "strategy_2_option_1","strategy_2_option_2")
write.xlsx(df_list,file = paste0(save.path,"GSEA_bulk.xlsx"),
           colNames = TRUE, borders = "surrounding",colWidths = c(NA, "auto", "auto"))

write.csv(df_fgseaRes, file = paste0(save.path,"fgsea_bulk_RNAseq.csv"),
          row.names = FALSE)
# scatter plot
data = df_fgseaRes[,c()]
