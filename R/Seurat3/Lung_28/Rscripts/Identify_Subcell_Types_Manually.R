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

#  ====== read annotation excel =================
annotation <- readxl::read_excel("doc/Harmony annotation Yang.xlsx")
(groups <- unique(annotation$UMAP))
(rds_files <-list.files("data", pattern = paste(groups,collapse = "|")))

args <- as.numeric(slurm_arrayid)
print(paste0("slurm_arrayid=",args))
(g <- groups[args])
(rds <- rds_files[args])
g_anno = annotation[annotation$UMAP %in% g,]
#======== rename ident =================
object = readRDS(file = paste0("data/",rds))
resolutions = g_anno["res"] %>% pull %>% sort
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
object[["cell.labels"]] = 0
meta.data = object@meta.data
meta.data["barcodes"] = rownames(meta.data)
for(i in seq_along(resolutions)){
        meta.data[meta.data[,paste0("SCT_snn_res.",resolutions[i])] %in% str2int(g_anno$cluster[i]),
                  "cell.labels"] = g_anno$`cell label`[i]
}
cell.labels = meta.data[meta.data[,"cell.labels"] != 0,c("barcodes","cell.labels")]
write.csv(cell.labels, file = paste0(path,sub("rds","csv",rds)),row.names = FALSE)
