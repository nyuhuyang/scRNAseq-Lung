#conda activate r4.0

library(Seurat)
library(dplyr)
library(tidyr)
library(kableExtra)
library(magrittr)
library(readxl)
library(cowplot)
library(stringr)
library(openxlsx)
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)
# load data
object = readRDS(file = "data/Lung_time_6_20201001.rds")
DefaultAssay(object)  = "SCT"
set.seed(101)
#======== rename ident =================
object %<>% FindClusters(resolution = 0.4)
object %<>% FindClusters(resolution = 2)

df_annotation <- readxl::read_excel("doc/Annotations/20201004_Lung_time_6 - Annotations.xlsx",
                                    sheet = "Sheet1")
resolutions = unique(df_annotation$Resolution) %>% .[complete.cases(.)]
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
                df_res$Celltype[m] = temp[grep(c, Cluster),"Label"]
        }
        df_res$Resolution = as.character(res)
        df_res_list[[i]] = df_res
}
df_res = bind_rows(df_res_list);rm(df_res_list);GC()
df_res = df_res[complete.cases(df_res$Cluster),]
#===========================================
meta.data = object@meta.data

meta.data$annotations3 = "Un"
for(m in 1:nrow(df_res)){
        resolution = paste0("SCT_snn_res.",df_res$Resolution[m])
        meta.data[meta.data[,resolution] %in% df_res$Cluster[m],
                  "annotations3"] = df_res$Celltype[m]
        Progress(m, nrow(df_res))
}

# label by expression
exp = FetchData(object, vars = c("SCGB1A1","SCGB3A2"))
meta.data[meta.data$SCT_snn_res.0.4 %in% c(0,6) & exp$SCGB3A2 > 1,
          "annotations3"] = "S-d"
object@meta.data = meta.data
Idents(object) = "annotations3"
object %<>% AddMetaColor(label= "annotations3", colors = Singler.colors)
UMAPPlot.1(object, group.by= "annotations3", do.print = T, label = T, label.repel = T)

meta.data = cbind(meta.data[,c("orig.ident","annotations3")],
                  meta.data[,-which(grepl(c("orig.ident|annotations3"),colnames(meta.data)))])
meta.data$`integrated_snn_res.0.6` = NULL
meta.data$`cell.type` = NULL
meta.data$`cell.type.colors` = NULL
meta.data = cbind(object@reductions$umap@cell.embeddings,meta.data[,c("orig.ident","annotations3",.)])
write.csv(meta.data, file = paste0(path, "Lung_time_6_coordinates.csv"))

saveRDS(object, file = paste0("data/Lung_time_6_20201001.rds"))


#======== test cluster 4 and 6 ==========
Idents(object) = "SCT_snn_res.0.4"
sub_object = subset(object, idents = c(0,6))
jpeg(paste0(path,"RidgePlot_SCGB3A2.jpeg"), 
     units="in", width=10, height=7,res=300)
RidgePlot(sub_object, features = "SCGB3A2", ncol = 1)
dev.off()

jpeg(paste0(path,"RidgePlot_SCGB3A2_all.jpeg"), 
     units="in", width=10, height=7,res=300)
RidgePlot(object, features = "SCGB3A2", ncol = 1)
dev.off()
exp = FetchData(sub_object, vars = c("SCGB1A1","SCGB3A2"))
sub_object@meta.data %<>% cbind(exp)
sub_object@meta.data$SCT_snn_res.0.4 %<>% as.character
table(sub_object$`SCT_snn_res.0.4`, sub_object$SCGB3A2 >1)
sub_object <- subset(sub_object, subset = SCGB3A2 >=1)
table(sub_object$`SCT_snn_res.0.4`,sub_object$`SCT_snn_res.2`) %>% t %>%
  as.data.frame.matrix() %>% write.table(file = paste0(path,"SCGB3A2.txt"))
