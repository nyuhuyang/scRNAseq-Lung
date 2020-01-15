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
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)

# load DEGs files
res =  read.csv("Yang/proximal_distal_terminal/Non-Integration/DEGs/cell_types/Lung_24-FC0_cell_types.csv",
                stringsAsFactors = F)
df_cell_types <- readxl::read_excel("doc/Cell type abbreviation.xlsx")
res$cluster %<>% plyr::mapvalues(
        from = df_cell_types$`Cell types`,
        to = df_cell_types$Abbreviation)
# Load GeneSets
readPathways <- function(folder){
        files <- list.files(folder)
        pathway_list <- vector(mode = "list", length = length(files))
        names(pathway_list) = gsub(".txt","",files,fixed = T) %>% 
                gsub("[0-9]+_","",.)
        for(i in seq_along(files)){
                tmp <- read.csv(paste(folder, files[i],sep = "/"), col.names = "gene")
                pathway_list[[i]] = stringr::str_split(tmp$gene, pattern = " /// |///") %>% 
                        unlist() %>% unique
        }
        return(pathway_list)
}

Development <- readPathways("doc/GeneSets/1-Development")
Physiology <- readPathways("doc/GeneSets/2-Physiology")
Disease <- readPathways("doc/GeneSets/3-Disease")
Cancer <- readPathways("doc/GeneSets/4-Cancer")

Development %>% head() %>% lapply(head)
Physiology %>% head() %>% lapply(head)
Disease %>% head() %>% lapply(head)
Cancer %>% head() %>% lapply(head)

# Now, run the fgsea algorithm with 1000 permutations:
fgseaRes = FgseaDotPlot(stats=res, pathways=Development,
                     padj = 0.25,pval = 0.05,
                     order.yaxis.by = c("BC","NES"),decreasing = F,
                     Rowv = T,Colv = T,
                     size = " -log10(pval)", fill = "NES",
                     pathway.name = "Development",rotate.x.text = T,
                     title = "multiple lung cells",
                     font.xtickslab=7, font.main=17, font.ytickslab = 14,
                     font.legend = list(size = 15),font.label = list(size = 15),
                     do.return = T,
                     width = 12,height = 6)
write.csv(fgseaRes, file = paste0(path,"Development_FDR0.25_pval0.05.csv"))

fgseaRes = FgseaDotPlot(stats=res, pathways=Physiology,
                        padj = 0.25,pval = 0.05,
                        order.yaxis.by = c("BC","NES"),decreasing = F,
                        Rowv = T,Colv = T,
                        size = " -log10(pval)", fill = "NES",
                        pathway.name = "Physiology",rotate.x.text = T,
                        title = "multiple lung cells",
                        font.xtickslab=7, font.main=17, font.ytickslab = 14,
                        font.legend = list(size = 15),font.label = list(size = 15),
                        do.return = T,
                        width = 12,height = 6)
write.csv(fgseaRes, file = paste0(path,"Physiology_FDR0.05_pval0.05.csv"))

fgseaRes = FgseaDotPlot(stats=res, pathways=Disease,
                        padj = 0.25,pval = 0.05,
                        order.yaxis.by = c("BC","NES"),decreasing = F,
                        Rowv = T,Colv = T,
                        size = " -log10(pval)", fill = "NES",
                        pathway.name = "Disease",rotate.x.text = T,
                        title = "multiple lung cells",
                        font.xtickslab=7, font.main=17, font.ytickslab = 14,
                        font.legend = list(size = 15),font.label = list(size = 15),
                        do.return = T,
                        width = 12,height = 4)
write.csv(fgseaRes, file = paste0(path,"Disease_FDR0.25_pval0.05.csv"))

fgseaRes = FgseaDotPlot(stats=res, pathways=Cancer,
                        padj = 0.25,pval = 0.05,
                        order.yaxis.by = c("BC","NES"),decreasing = F,
                        Rowv = T,Colv = T,
                        size = " -log10(pval)", fill = "NES",
                        pathway.name = "Cancer",rotate.x.text = T,
                        title = "multiple lung cells",
                        font.xtickslab=7, font.main=17, font.ytickslab = 12,
                        font.legend = list(size = 15),font.label = list(size = 15),
                        do.return = T,
                        width = 12,height = 6)
write.csv(fgseaRes, file = paste0(path,"Cancer_FDR0.25_pval0.05.csv"))
