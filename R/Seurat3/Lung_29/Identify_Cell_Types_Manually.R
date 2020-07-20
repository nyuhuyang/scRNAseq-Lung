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
object1 = readRDS(file = "data/Lung_29_20200617.rds")
DefaultAssay(object) = "SCT"
#======== rename ident =================
object %<>% FindClusters(resolution = 4.9)
object %<>% FindClusters(resolution = 4.8)
object %<>% FindClusters(resolution = 3.5)
object %<>% FindClusters(resolution = 2)

df_annotation <- readxl::read_excel("doc/Annotations/20200625_Annotations.xlsx", sheet = "Sheet1")
resolutions = unique(df_annotation$Resolution)
df_res_list <- list()
for(i in seq_along(resolutions)){
        res = resolutions[i]
        temp <- df_annotation[df_annotation$Resolution %in% res,]
        Cluster <- stringr::str_split(string = temp$Cluster, 
                                        pattern = "\\+")
        df_res <- data.frame("Cluster" = as.integer(unlist(Cluster)), stringsAsFactors = F)
        for(m in 1:nrow(df_res)) {
                c = as.character(df_res$Cluster[m])
                c = paste0("\\b",c,"\\b")
                df_res$Celltype[m] = temp[grep(c, Cluster),"Cell type"]
        }
        df_res$Resolution = as.character(res)
        df_res_list[[i]] = df_res
}
df_res = bind_rows(df_res_list)

#===========================================
meta.data = cbind(object@meta.data,  object@reductions$umap@cell.embeddings)

meta.data$annotations3 = "unknown"
for(m in 1:nrow(df_res)){
        resolution = paste0("SCT_snn_res.",df_res$Resolution[m])
        meta.data[meta.data[,resolution] %in% df_res$Cluster[m],
                  "annotations3"] = df_res$Celltype[m]
        Progress(m, nrow(df_res))
}
#meta.data[meta.data$annotations3 %in% "S-d" & meta.data$UMAP_1 < 5,"annotations3"] = "S"
meta.data[meta.data$annotations3 %in% "IC2+" & meta.data$`SCT_snn_res.2` %in% 62,
          "annotations3"] = "IC2"

# label by expression
exp = FetchData(object, vars = c("KRT5","KRT15","CDH5","CLDN5","CD68","CD163","DCN","ACTA2",
                                 "CD3E", "CD3G","SFTPC","NAPSA","MKI67",
                                 "SCGB3A2","SFTPB","KRT14","UBE2C","SERPINB4","SERPINB3","S100A2"))
#meta.data$annotations3 %<>% gsub("^S-.*","S",.)
S <- grepl("^S$",meta.data$annotations3)
meta.data[S & (exp$SCGB3A2 > 1 | exp$SFTPB > 1),"annotations3"] = "S-d"

SM1 <- meta.data$annotations3 %in% c("SM1")
meta.data[SM1 & (exp$KRT5 > 0 | exp$KRT14 > 0),"annotations3"] = "MEC"
meta.data[SM1 & meta.data$UMAP_2 > 7,"annotations3"] = "F2"
meta.data[meta.data$annotations3 %in% "En-SM" & meta.data$UMAP_1 > 0 &
                  (exp$ACTA2 > 0 | exp$DCN > 0),"annotations3"] = "F3"
meta.data[meta.data$annotations3 %in% "C1" & 
                  (meta.data$UMAP_2 > 3 & meta.data$UMAP_2 < 7),"annotations3"] = "Mixed"

# label by corrdinates
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

# BC-p
BC_dn <- (meta.data$annotations3 %in% "BC-p") & (exp$MKI67 == 0) & (exp$UBE2C == 0)
meta.data[BC_dn & (exp$KRT5 > 0 | exp$KRT15 > 0),"annotations3"] = "BC"
meta.data[BC_dn & exp$KRT5 == 0 & exp$KRT15 == 0 & exp$SERPINB4 > 0 ,"annotations3"] = "IC2"
meta.data[BC_dn & exp$KRT5 == 0 & exp$KRT15 == 0 & exp$SERPINB4 == 0 &
                  (exp$SERPINB3 >0 | exp$S100A2 >0),"annotations3"] = "IC1"
# M0
Idents(object) = "annotations3"
M <- subset(object, idents = c("Mon","M0","M1","M2"))
Idents(M) %<>% factor(levels = c("Mon","M2","M1","M0"))
features = c("SPP1","CXCL9","CXCL10","CXCL11","CXCL12","MT1H","MT1G","CCL7","IDO1",
             "CCL18","FABP4","MSR1","APOC1","APOE","FBP1","GRN","C1QB","C1QA")
features = unique(c("IL1B","S100A9","EREG","CXCL3","SERPINB2","CD163","MARCO","CD68","HMOX1",
             "SPP1","MT1G","CXCL10","CXCL9","APOC1","C1QB","C1QA","APOE","CCL18","FABP4",
             "GRN","ACP5","FBP1"))
jpeg(paste0(path,"RidgePlot.jpeg"), units="in", width=10, height=10,res=600)
print(RidgePlot(M, features = features[1:9], ncol = 3))
dev.off()
jpeg(paste0(path,"RidgePlot~.jpeg"), units="in", width=10, height=10,res=600)
print(RidgePlot(M, features = features[10:18], ncol = 3))
dev.off()
jpeg(paste0(path,"RidgePlot~~.jpeg"), units="in", width=10, height=10,res=600)
print(RidgePlot(M, features = features[19:22], ncol = 2))
dev.off()

M1_test <- rowSums((exp[,features[1:9]] > 0)*1L)
M@meta.data[names(M1_test)[M1_test == 0],"annotations3"] = "M1"



table(meta.data$annotations3, meta.data$SCT_snn_res.4.9)

#  unclear label as Mixed
meta.data$annotations3 %<>% gsub("Proliferating|unknown","Mixed",.)

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

T_reg <- subset(object, idents = "T-reg")
exp1 = FetchData(T_reg,vars = c("FOXP3","IL2RA","TIGIT","CTLA4"))
exp1$criteria_1 = exp1$FOXP3 > 0
exp1$criteria_2 = exp1$FOXP3 > 0 | exp1$IL2RA > 0 |exp1$TIGIT > 0 |exp1$CTLA4 > 0
exp1$orig.ident = gsub("_.*","",rownames(exp1))
for(c in c("criteria_1", "criteria_2")){
        df <- table(exp1[,c], exp1$orig.ident) %>% 
                as.data.frame()
        colnames(df) = c(c,"samples","Freq")
        df %<>% spread("samples","Freq")
        write.csv(df, paste0(path,c,"_Treg.csv"))
}


# - Table: number of cells per cell types (per each sample and total)
df <- table(object$annotations3, object$orig.ident) %>% 
        as.data.frame()
colnames(df) = c("cell.types","samples","Freq")
df %<>% spread("samples","Freq")
rownames(df) = df$cell.types
df = df[order(df$cell.types),]

#write.csv(df, paste0(path,"Lung_24-",con,"_cell.types_by_samples.csv"))
write.csv(df, paste0(path,"Cell_types_by_samples~.csv"))