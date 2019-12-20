########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################

library(dplyr)
library(tidyr)
library(kableExtra)
library(magrittr)
library(gplots)
library(ggplot2)
library(fgsea)
library(tibble)
library(ggpubr)
library(ggsci)
source("../R/Seurat3_functions.R")
path <- paste0("Yang/proximal_distal_terminal/Non-Integration/GSEA")
if(!dir.exists(path)) dir.create(path, recursive = T)
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
DEG_path <- "Yang/proximal_distal_terminal/Non-Integration/DEGs/cell_types_pairwise/SCT/"

# 3.1.1 load plan
df_labels_Pairs <- readxl::read_excel("doc/DEG cell types comparison plan.xlsx")
df_labels_Pairs = df_labels_Pairs[df_labels_Pairs$`GSEA` %in% "GSEA",]

(ident1 <- strsplit(df_labels_Pairs$`Group 1`, "\\+"))

csv <- function(x, y) paste0(DEG_path,"Lung_24_SCT_pairwise_",x,"_",paste(y,collapse = "."),".csv")
csv_files <- mapply(csv, df_labels_Pairs$Comparison, ident1)

CSV_files = list.files(DEG_path)
table(csv_files %in% paste0(DEG_path,CSV_files))

hallmark <- gmtPathways("../seurat_resources/msigdb/h.all.v7.0.symbols.gmt")
c2mark <- gmtPathways("../seurat_resources/msigdb/c2.all.v7.0.symbols.gmt")
tftmark <- gmtPathways("../seurat_resources/msigdb/c3.tft.v7.0.symbols.gmt")

names(hallmark) = gsub("_"," ",names(hallmark))
names(c2mark) = gsub("_"," ",names(c2mark))
names(tftmark) = gsub("_"," ",names(tftmark))

hallmark %>% head() %>% lapply(head)
c2mark %>% head() %>% lapply(head)
tftmark %>% head() %>% lapply(head)

res <- list()
for(i in seq_along(csv_files)){
        res[[i]] = read.csv(file = csv_files[i],row.names = 1, stringsAsFactors=F)
        res[[i]]$gene = rownames(res[[i]])
        res[[i]]$Group1.vs.Group2 = gsub("_.*","",res[[i]]$Group1.vs.Group2)
        res[[i]]$Group1.vs.Group2 = paste0(names(csv_files[i]),"_",res[[i]]$Group1.vs.Group2)
        res[[i]] = res[[i]][order(res[[i]]$avg_logFC, decreasing = T),]
}
res = bind_rows(res)
(clusters <- unique(res$Group1.vs.Group2))
colnames(res)[grep("Group1.vs.Group2",colnames(res))] = "clusters"
# Now, run the fgsea algorithm with 1000 permutations:
fgseaRes = FgseaDotPlot(stats=res, pathways=hallmark, nperm=1000,
                     padj = 0.25,pval = 0.05,
                     order.yaxis.by = c("3_BC","NES"),decreasing = F,
                     order.x = unique(res$clusters),
                     size = " -log10(pval)", fill = "NES",
                     pathway.name = "Hallmark",rotate.x.text = T,
                     title = "multiple lung cells",
                     font.xtickslab=12, font.main=17, font.ytickslab = 8,
                     font.legend = list(size = 15),font.label = list(size = 15),
                     do.return = T,
                     width = 10,height = 7)
write.csv(fgseaRes, file = paste0(path,"Hallmark_FDR0.25_pval0.05.csv"))

fgseaRes = FgseaDotPlot(stats=res, pathways=c2mark, nperm=1000,
                        padj = 0.05,pval = 0.05,
                        order.by = c("3_BC","NES"),decreasing = F,
                        order = unique(res$clusters),
                        size = " -log10(pval)", fill = "NES",
                        pathway.name = "All curated gene sets",rotate.x.text = T,
                        title = "multiple lung cells",
                        font.xtickslab=12, font.main=17, font.ytickslab = 1,
                        font.legend = list(size = 15),font.label = list(size = 15),
                        do.return = T,
                        width = 8,height = 50)
write.csv(fgseaRes, file = paste0(path,"C2mark_FDR0.05_pval0.05.csv"))

fgseaRes = FgseaDotPlot(stats=res, pathways=tftmark, nperm=1000,
                        padj = 0.25,pval = 0.05,
                        order.by = c("3_BC","NES"),decreasing = F,
                        order = unique(res$clusters),
                        size = " -log10(pval)", fill = "NES",
                        pathway.name = "Transcription factor targets",rotate.x.text = T,
                        title = "multiple lung cells",
                        font.xtickslab=12, font.main=17, font.ytickslab = 1,
                        font.legend = list(size = 15),font.label = list(size = 15),
                        do.return = T,
                        width = 8,height = 20)
write.csv(fgseaRes, file = paste0(path,"tft_FDR0.25_pval0.05.csv"))
