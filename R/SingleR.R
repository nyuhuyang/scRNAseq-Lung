library(SingleR)
library(Seurat)
library(reshape2)
library(pheatmap)
library(kableExtra)
library(dplyr)
library(tidyr)
source("../R/Seurat_functions.R")
source("../R/SingleR_functions.R")
path <- paste0(getwd(),"/output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
#====== 3.1 Create Singler Object  ==========================================
(load(file="data/Lung_Harmony_3_20190109.Rda"))
object@scale.data = NULL
GC();GC();GC();GC();GC();GC();GC();GC();GC();GC();GC();GC();GC();GC();
singler = CreateSinglerObject(as.matrix(object@data), annot = NULL, project.name="EC-SR-5542",
                              min.genes = 500,technology = "10X", species = "Human", citation = "",
                              ref.list = list(),normalize.gene.length = F, variable.genes = "de",
                              fine.tune = F, do.signatures = F, clusters = NULL,
                              numCores = SingleR.numCores/2)
# if singler didn't find all cell labels
print(length(singler$singler[[1]]$SingleR.single$labels) == ncol(object@data))
if(length(singler$singler[[1]]$SingleR.single$labels) != ncol(object@data)){
        all.cell = object@cell.names;length(all.cell)
        know.cell = rownames(singler$singler[[1]]$SingleR.single$labels);length(know.cell)
        object = SubsetData(object, cells.use = know.cell)
}
singler$meta.data$orig.ident = object@meta.data$orig.ident # the original identities, if not supplied in 'annot'
singler$meta.data$xy = object@dr$tsne@cell.embeddings # the tSNE coordinates
singler$meta.data$clusters = object@ident # the Seurat clusters (if 'clusters' not provided)
save(singler,file="./output/singler_Lung_3F_20190109.Rda")
