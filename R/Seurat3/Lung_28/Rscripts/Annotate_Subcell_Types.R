########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################
invisible(lapply(c("Seurat","dplyr","kableExtra","cowplot",
                   "magrittr","MAST","future","ggplot2","tidyr"), function(x) {
                           suppressPackageStartupMessages(library(x,character.only = T))
                   }))
source("../R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)

# Need 32 GB
#  ====== read annotation excel =================
annotation <- readxl::read_excel("doc/Harmony annotation Yang.xlsx")
(groups <- unique(annotation$UMAP) %>% sort)
(rds_files <-list.files("data", pattern = paste(groups,collapse = "|")))

# SLURM_ARRAY_TASK_ID
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
if (length(slurm_arrayid)!=1)  stop("Exact one argument must be supplied!")
# coerce the value to an integer
args <- as.numeric(slurm_arrayid)
print(paste0("slurm_arrayid=",args))
(g <- groups[args])
(rds <- rds_files[args])
g_anno = annotation[annotation$UMAP %in% g,]
#======== rename ident =================
object = readRDS(file = paste0("data/",rds))
resolutions = g_anno["res"] %>% pull
(res = unique(resolutions))
for(i in seq_along(res)) object %<>% FindClusters(resolution = res[i])
str2int <- function(char){
        int <- gsub(",","",char) %>% strsplit(" ") %>% unlist %>% as.integer()
        return(int)
}
if(g == "Global"){
        Idents(object) = "SCT_snn_res.0.2"
        sub_object <- subset(object, idents = 6)
        Idents(sub_object) = "conditions"
        sub_object %<>% subset(idents = "proximal")
        
        object[["cell.labels"]] = 0
        meta.data = object@meta.data
        meta.data["barcodes"] = rownames(meta.data)
        meta.data[colnames(sub_object),"cell.labels"] = g_anno$`cell label`[1]
        for(i in 2:length(resolutions)){
                meta.data[meta.data[,paste0("SCT_snn_res.",resolutions[i])] %in% g_anno$cluster[i],
                          "cell.labels"] = g_anno$`cell label`[i]
        }
}
if(g != "Global"){
        object[["cell.labels"]] = 0
        meta.data = object@meta.data
        meta.data["barcodes"] = rownames(meta.data)
        for(i in seq_along(resolutions)){
                meta.data[meta.data[,paste0("SCT_snn_res.",resolutions[i])] %in% str2int(g_anno$cluster[i]),
                          "cell.labels"] = g_anno$`cell label`[i]
        }
}
if(g == "21-ambigous-pca"){
        object %<>% AddModuleScore(features = list(c("ACTA2","ELN","COL1A1","BGN")), name = "F5")
        colnames(object@meta.data)[grep("F51",colnames(object@meta.data))] = "ACTA2+ELN+COL1A1+BGN"
        
        m <- cbind(meta.data, object[["umap"]]@cell.embeddings, object[["ACTA2+ELN+COL1A1+BGN"]])
        m[m$UMAP_1 > -2.5 & m$UMAP_1 < -1.8 & m$UMAP_2 > -8 & m$UMAP_2 < -7.2,
          "cell.labels"] = "Cr"
        m[((m$UMAP_1 > -4.6 & m$UMAP_1 < -4.1 & m$UMAP_2 < -4.7 & m$UMAP_2 > -5.5) |
                   (m$UMAP_1 > -3.8 & m$UMAP_1 < -3.1 & m$UMAP_2 < -6 & m$UMAP_2 > -6.5)) &
                  m$`ACTA2+ELN+COL1A1+BGN` > 0.5,
          "cell.labels"] = "F5"  
        m[m$UMAP_1 > 0.7 & m$UMAP_1 < 1.3 & m$UMAP_2 > -4.6 & m$UMAP_2 < -3.9,
          "cell.labels"] = "DC"
        meta.data = m
}
cell.labels = meta.data[meta.data[,"cell.labels"] != 0,c("barcodes","cell.labels")]
write.csv(cell.labels, file = paste0(path,sub("0206","0219",sub("rds","csv",rds))),row.names = FALSE)
