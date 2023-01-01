library(Seurat)
library(dplyr)
library(magrittr)
library(readxl)
library(stringr)
library(openxlsx)
library(tidyverse)
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)

#-======================
#after obtain DEGs from human lung HLCA v2
library(arrow)
annotations.path <- "data/human_lung_v2/data/annotations.parquet"
annotations <- read_parquet(annotations.path)
Secretory <- unique(annotations[annotations$ann_level_3 %in% c("Secretory"),"ann_finest_level"])

#-----ven diagram-------
DEGsv2 <- read_excel("output/20221013/2022-10-13-ann_finest_level.xlsx") %>%
    filter(p_val_adj < 0.05)
colnames(DEGsv2) =c("cluster","gene","scores","avg_log2FC","p_val","p_val_adj","pct.1","pct.2")
DEGsv2 =DEGsv2[,c("p_val","avg_log2FC","pct.1","pct.2","p_val_adj","cluster","gene")]
openxlsx::write.xlsx(DEGsv2, file =  "output/20221013/2022-10-13-ann_finest_level.filter.xlsx",
                     colNames = TRUE,rowNames = F,borders = "surrounding")

DEGsv2 <- read_excel("output/20221013/2022-10-13-ann_finest_level.xlsx") %>% 
    filter(p_val_adj < 0.05 & pts > 0.1 & avg_log2FC >1 & group %in% Secretory)  %>%
    group_by(group) %>%
    top_n(100,scores) %>%
    ungroup() %>%
    as.data.frame()
Secretory_pos_genes <- split(DEGsv2,DEGsv2$group) %>% 
    lapply(function(df) df[(df$avg_log2FC > 1),"genes"]) %>% 
    lapply(unique)


DEGs <- read_excel("Yang/Lung_30/hg38/DE_analysis/Lung_30_DEG_Cell.category.xlsx") %>% as.data.frame()
DEGs[DEGs$p_val == 0,"p_val"] = min(DEGs[DEGs$p_val > 0,"p_val"], .Machine$double.xmin)
pos_genes <-  DEGs %>%
    filter(p_val_adj < 0.05 & cluster == "TASC") %>% 
    head(100) %>%
    .[,"gene"] %>% unique
Secretory_pos_genes %<>% c(list("TASC"=pos_genes))

for(cell in Secretory){
    euler_df <- eulerr::euler(Secretory_pos_genes[c("TASC",cell)],shape = "circle")
    
    g <- plot(euler_df, quantities = TRUE, lty = 1:6,
              legend = TRUE)
    jpeg(paste0(path,"/","Venn_",cell,"_TASC_top100.jpeg"), units="in", width=10, height=7,res=600)
    print(g)
    dev.off()
    Progress(which(Secretory %in% cell), length(Secretory))
}

library(ggvenn)
ggvenn(
    Secretory_pos_genes, 
    #fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
    stroke_size = 0.5, set_name_size = 8
)
x <- list(
    A = sample(genes,300), 
    B = sample(genes,525), 
    C = sample(genes,440),
    D = sample(genes,350),
    E = sample(genes,250),
    EF = sample(genes,250)
)
venn.diagram(Secretory_pos_genes[1,], filename = paste0(path,"venn_8.png"))

library(VennDiagram)
openxlsx::write.xlsx(DEGsv2, file =  "output/20221013/2022-10-13-ann_finest_level.Secretory.xlsx",
                     colNames = TRUE,rowNames = F,borders = "surrounding")

DEGs_AT2 <- read_excel("output/20221014/2022-10-14-Transitional Club-AT2 vs AT2 in   .xlsx") %>% as.data.frame() %>%
    filter(p_val_adj < 0.05) %>% 
    head(100) %>%
    .[,"gene"] %>% unique
DEGs_Secretory <- read_excel("output/20221014/2022-10-14-Transitional Club-AT2 vs Club (nasal)_Club (non-nasal)_Goblet (bronchial)_Goblet (nasal)_Goblet (subsegmental) in   .xlsx") %>% as.data.frame() %>%
    filter(p_val_adj < 0.05) %>% 
    head(100) %>%
    .[,"gene"] %>% unique

Secretory_AT2_pos_genes <- list("TASC"=pos_genes,
                                "Transitional Club-AT2 vs AT2" = DEGs_AT2,
                                "Transitional Club-AT2 vs other Secretory" = DEGs_Secretory)
euler_df <- eulerr::euler(Secretory_AT2_pos_genes,shape = "circle")

g <- plot(euler_df, quantities = TRUE, lty = 1:6,
          legend = TRUE)
file.name = "Venn_Secretory_AT2_TASC_top100.jpeg"
jpeg(paste0(path,"/",file.name), units="in", width=10, height=7,res=600)
print(g)
dev.off()



#-----fgsea dotplot-------
# all
DEGsv2 <- read_excel("output/20221013/2022-10-13-ann_finest_level.xlsx") %>% 
    filter(p_val_adj < 0.05 & pts > 0.1 & avg_log2FC >2)  %>%
    group_by(group) %>%
    top_n(250) %>%
    ungroup()
    
pos_genes <- split(DEGsv2,DEGsv2$group) %>% 
    lapply(function(df) df$genes) %>% 
    lapply(unique)
str(pos_genes)
res <- read_excel("Yang/Lung_30/hg38/DE_analysis/Lung_30_DEG_Cell.category.xlsx") %>% 
    filter(pct.1 > 0.1 & avg_log2FC >0) %>% 
    as.data.frame() 
plot <- FgseaDotPlot(stats=res, pathways=pos_genes,Rowv = TRUE,Colv = TRUE,
             cols = c(rev(brewer.pal(10,"RdBu"))), return.plot = TRUE)
plot <- plot + theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
             axis.text.x =    element_text(angle = 90,
                                           hjust = 0,
                                           vjust =  0.5)
             )
jpeg(paste0(path, "/Dotplot_log2FC_2_fdr0.25_p0.05.jpeg"),units="in", width=9, height=9,res=600)
print(plot)
dev.off()


#---
annotations.path <- "data/human_lung_v2/data/annotations.parquet"
annotations <- read_parquet(annotations.path)
Epi <- unique(annotations[annotations$ann_level_1 %in% "Epithelial","ann_finest_level"])
Epi_pos_genes <- pos_genes[Epi]

meta.data = readRDS("output/Lung_30_20210831_metadata_v2.rds")
Lung_Epi <- unique(meta.data[meta.data$Superfamily %in% "Epi","Cell_subtype"]) %>% as.character()

plot1 <- FgseaDotPlot(stats= res[res$cluster %in% Lung_Epi,], pathways=Epi_pos_genes,Rowv = TRUE,Colv = TRUE,cols = rev(brewer.pal(10,"RdBu")), return.plot = TRUE,font.main = 20)

plot1 <- plot1 + theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
                     axis.text.x =    element_text(angle = 90,
                                                   hjust = 0,
                                                   vjust =  0.5)
)
jpeg(paste0(path, "/Dotplot_Epi_log2FC_2_fdr0.25_p0.05.jpeg"),units="in", width=5, height=5 ,res=600)
print(plot1)
dev.off()
