library(Seurat)
library(dplyr)
library(tidyr)
library(kableExtra)
library(magrittr)
library(readxl)
library(cowplot)
library(stringr)
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)
# load data
object = readRDS(file = "data/Lung_29_20200617.rds")
#======== rename ident =================
object %<>% FindClusters(resolution = 3.5)
object %<>% FindClusters(resolution = 2)

df_annotation <- readxl::read_excel("doc/Annotations/20200619_Annotations.xlsx")
# Resolution = 2 
res_2 <- df_annotation[df_annotation$Resolution %in% 2,]
Cluster_2 <- stringr::str_split(string = res_2$Cluster, 
                                pattern = "\\+")
df_res_2 <- data.frame("Cluster" = as.integer(unlist(Cluster_2)), stringsAsFactors = F)
for(m in 1:nrow(df_res_2)) {
        c = as.character(df_res_2$Cluster[m])
        c = paste0("\\b",c,"\\b")
        df_res_2$Celltype[m] = res_2[grep(c, Cluster_2),"Cell type"]
}
# Resolution = 3.5
res_3 <- df_annotation[df_annotation$Resolution %in% 3.5,]
Cluster_3 <- stringr::str_split(string = res_3$Cluster, 
                                pattern = "\\+")
df_res_3 <- data.frame("Cluster" = as.integer(unlist(Cluster_3)), stringsAsFactors = F)
for(m in 1:nrow(df_res_3)) {
        c = as.character(df_res_3$Cluster[m])
        c = paste0("\\b",c,"\\b")
        df_res_3$Celltype[m] = res_3[grep(c, Cluster_3),"Cell type"]
}
#===========================================
meta.data = cbind(object@meta.data,  object@reductions$umap@cell.embeddings)
# Resolution = 2 

meta.data$annotations3 = as.character(meta.data$SCT_snn_res.2)
for(m in 1:nrow(df_res_2)){
        meta.data[meta.data$annotations3 %in% df_res_2$Cluster[m],
                  "annotations3"] = df_res_2$Celltype[m]
}
meta.data[meta.data$annotations3 %in% "S-d" & meta.data$UMAP_1 < 5,"annotations3"] = "S"
# Resolution = 3.5
for(m in 1:nrow(df_res_2)){
        meta.data[meta.data$SCT_snn_res.3.5 %in% df_res_3$Cluster[m],
                  "annotations3"] = df_res_3$Celltype[m]
}
#meta.data[!(meta.data$annotations3 %in% unique(df_annotation$`Cell type`)),"annotations3"] = "Mixed"

exp = FetchData(object, vars = c("KRT5","KRT15","CDH5","CLDN5","CD68","CD163","DCN","ACTA2",
                                 "CD3E", "CD3G","SFTPC","NAPSA","MKI67",
                                 "SCGB3A2","SFTPB","KRT14"))
#Pro <- meta.data$annotations3 %in% "Proliferating"
#meta.data[Pro & (exp$KRT5 > 0 | exp$KRT15 > 0),"annotations3"] = "BC-p"
#meta.data[Pro & (exp$CDH5 > 0 | exp$CLDN5 > 0),"annotations3"] = "En-p"
#meta.data[Pro & (exp$CD68 > 0 | exp$CD163 > 0),"annotations3"] = "M-p"
#meta.data[Pro & (exp$DCN > 0 | exp$ACTA2 > 0),"annotations3"] = "Str-p"
#meta.data[Pro & (exp$CD3E > 0 | exp$CD3G > 0),"annotations3"] = "T-p"
#meta.data[Pro & (exp$SFTPC > 0 | exp$NAPSA > 0),"annotations3"] = "AT2-p"

meta.data$annotations3 %<>% gsub("^S-.*","S",.)
S <- meta.data$annotations3 %in% c("S")
meta.data[S & (exp$SCGB3A2 > 1 | exp$SFTPB > 1),"annotations3"] = "S-d"

SM1 <- meta.data$annotations3 %in% c("SM1")
meta.data[SM1 & (exp$KRT5 > 0 | exp$KRT14 > 0),"annotations3"] = "MEC"
meta.data[SM1 & meta.data$UMAP_2 > 7,"annotations3"] = "F2"
meta.data[meta.data$annotations3 %in% "En-SM" & meta.data$UMAP_1 > 0 &
                  (exp$ACTA2 > 0 | exp$DCN > 0),"annotations3"] = "F3"
meta.data[meta.data$annotations3 %in% "C1" & 
                  (meta.data$UMAP_2 > 3 & meta.data$UMAP_2 < 7),"annotations3"] = "Mixed"

# select by corrdinates
g <- UMAPPlot.1(object, group.by = "annotations3", cols = Singler.colors, 
                label = T, label.repel = T, no.legend = T,do.return = T)+
        rectangle(-8, 1, 7, 13,colour = "black")+ # En-p
        rectangle(-12, -4, -10, -2,colour = "black")+ # T-p
        rectangle(-5, 1, -3, 5,colour = "black")+ # M-p
        rectangle(1.5, 8, 3.7, 13,colour = "black")+ # Str-p
        rectangle(2, 8, -8, -1,colour = "black")+ # BC-p
        rectangle(7, 12.5, -13, -8,colour = "black") # AT2-p

jpeg(paste0(path,"UMAP~.jpeg"), units="in", width=10, height=10,res=600)
print(g)
dev.off()

rectangle_name <- function(df, x_left, x_right, y_bottom, y_top, cell.type){
        tmp = df[,"UMAP_1"] > x_left & df[,"UMAP_1"] < x_right & 
                df[,"UMAP_2"] > y_bottom & df[,"UMAP_2"] < y_top & 
                df$annotations3 %in% "Proliferating"
        df[rownames(df)[tmp], "annotations3"] = cell.type
        return(df)
}

meta.data %<>% rectangle_name(-8, 1, 7, 13, "En-p")
meta.data %<>% rectangle_name(-12, -4, -10, -2, "T-p")
meta.data %<>% rectangle_name(-5, 1, -3, 5, "M-p")
meta.data %<>% rectangle_name(1.5, 8, 3.7, 13, "Str-p")
meta.data %<>% rectangle_name(2, 8, -8, -1, "BC-p")
meta.data %<>% rectangle_name(7, 12.5, -13, -8, "AT2-p")

#meta.data$annotations3 %<>% gsub("Proliferating","Mixed",.)

object[["annotations3"]] = meta.data$annotations3
Idents(object)= "annotations3"

UMAPPlot.1(object,group.by = "annotations3", cols = Singler.colors, 
           label = T, label.repel = T, no.legend = T,do.return = F, do.print = T)

PrepareShiny(sub_object, samples, Rshiny_path, split.by = "annotations3", 
             reduction = "umap",verbose = T)
saveRDS(object, file = paste0("data/Lung_29_20200617.rds"))

pro <- meta.data$annotations3 %in% "Proliferating"

Str <- meta.data[,"UMAP_1"] > 1.5 & meta.data[,"UMAP_1"] < 8 & 
        meta.data[,"UMAP_2"] > 3.7 & meta.data[,"UMAP_2"] < 13
Str <- subset(object, cells = rownames(meta.data)[Str],)
Pro <- subset(object, cells = rownames(meta.data)[pro],)

# - Table: number of cells per cell types (per each sample and total)
df <- table(object$annotations3, object$orig.ident) %>% 
        as.data.frame()
colnames(df) = c("cell.types","samples","Freq")
df %<>% spread("samples","Freq")
rownames(df) = df$cell.types
df = df[order(df$cell.types),]

#write.csv(df, paste0(path,"Lung_24-",con,"_cell.types_by_samples.csv"))
write.csv(df, paste0(path,"Cell_types_by_samples.csv"))