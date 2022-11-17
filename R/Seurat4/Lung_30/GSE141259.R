invisible(lapply(c("Seurat","dplyr","ggplot2","scater","magrittr","pbapply",
                   "cowplot"), function(x) {
                       suppressPackageStartupMessages(library(x,character.only = T))
                   }))
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_functions.R")

slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
if (length(slurm_arrayid)!=1)  stop("Exact one argument must be supplied!")
# coerce the value to an integer
args <- as.integer(as.character(slurm_arrayid))
print(paste0("slurm_arrayid=",args))

#---------
datasets <- c("WholeLung","HighResolution")[args]
read_path <- paste0("data/annotation_references/GSE141259/",datasets)

#WholeLung
file.rename(from = paste0(read_path,"/GSE141259_",datasets,"_rawcounts.mtx.gz"),
            to  =paste0(read_path,"/matrix.mtx.gz"))
file.rename(from = paste0(read_path,"/GSE141259_",datasets,"_barcodes.txt.gz"),
            to  =paste0(read_path,"/barcodes.tsv.gz"))
file.rename(from = paste0(read_path,"/GSE141259_",datasets,"_genes.txt.gz"),
            to  =paste0(read_path,"/features.tsv.gz"))
"@shell
gunzip matrix.mtx.gz
#remove % at 2nd line
gzip -c matrix.mtx > matrix.mtx.gz
"

object <- Read10X(read_path, gene.column=1)
#line 34 in Read10X data <- readMM(file = matrix.loc); data <- t(data) # for HighResolution
meta.data <- data.table::fread(paste0(read_path,"/GSE141259_",datasets,"_cellinfo.csv.gz")) %>%
    as.data.frame() %>% tibble::column_to_rownames("cell_barcode")

table(colnames(object) == rownames(meta.data))
object <- CreateSeuratObject(counts = object,meta.data = meta.data)

object %<>% NormalizeData(verbose = TRUE) %>%
    FindVariableFeatures(verbose = TRUE) %>%
    ScaleData(verbose = TRUE) %>%
    RunPCA(verbose = TRUE,npcs = 50) %>%
    RunUMAP(verbose = TRUE, dims = 1:50) %>%
    FindNeighbors(reduction = "pca",dims = 1:50) %>%
    FindClusters(resolution = 0.8, algorithm = 1,verbose = F)

object[["umap"]]@cell.embeddings[,"UMAP_1"] = object$umap_1
object[["umap"]]@cell.embeddings[,"UMAP_2"] = object$umap_2
UMAPPlot(object, group.by = "cell_type")
saveRDS(object, file = paste0("data/GSE141259_",datasets,".rds"))
#-------
library(Seurat)
library(ShinyCell)
library(magrittr)
library(SeuratData)
library(SeuratDisk)
library(dplyr)
object$mouse <- object$orig.ident

if(datasets =="WholeLung"){
    meta.to.include =c("cell.type","spline_cluster","metacelltype","grouping",
                       "orig.ident","nCount_RNA","nFeature_RNA","nGene","nUMI",
                       "identifier","res.2","RNA_snn_res.0.8","mouse")
}
if(datasets == "HighResolution"){
    meta.to.include = c("cell_type","meta_celltype","orig.ident","nCount_RNA","nFeature_RNA", 
                        "identifier","sample_id","time_point",
                        "louvain_cluster","percent_mito","n_counts","n_genes",
                        "RNA_snn_res.0.8" ,"mouse")
}
shiny.dir = paste0("shinyApp/GSE141259_",datasets,"/")
table(meta.to.include %in% colnames(object@meta.data))
meta.to.include[!meta.to.include %in% colnames(object@meta.data) ]
scConf = createConfig(object, meta.to.include =meta.to.include,maxLevels = 300)
makeShinyApp(object, scConf, gex.assay = "RNA",gene.mapping = TRUE,
             gex.slot = "data",default.gene1 = "Scgb1a1",default.gene2 = "Krt15",
             default.multigene = c("Krt15","Muc5ac","Scgb1a1","Scgb3a2","Sftpb","Foxa2" ),
             default.dimred = c("UMAP_1","UMAP_2"),shiny.dir = shiny.dir,
             shiny.title = switch(datasets,
                                  "WholeLung"="GSE141259 whole mouse lung",
                                  "HighResolution"="GSE141259 High Resolution mouse lung")
             )

max_exp = qlcMatrix::rowMax(object@assays[["RNA"]]@data)
max_exp_df = data.frame("val"= as.vector(max_exp),row.names = rownames(object))
saveRDS(max_exp_df,paste0(shiny.dir,"sc1maxlvl.rds"))

sc1def = readRDS(paste0(shiny.dir,"sc1def.rds"))
sc1def$grp1 = "cell.type"
sc1def$grp2 = "orig.ident"
saveRDS(sc1def,paste0(shiny.dir,"sc1def.rds"))


format(object.size(object@assays$SCT),unit = "GB")
format(object.size(object@assays$integrated),unit = "GB")
object[["RNA"]]@counts = matrix(0,0,0)
object[["RNA"]]@scale.data = matrix(0,0,0)

Class = sapply(meta.data,class)

for(col in  names(Class)[Class %in% "factor"]){
    meta.data[,col] %<>% as.character()
}

object@meta.data = meta.data[,meta.to.include]


format(object.size(object),unit = "GB")

file.remove(paste0(shiny.dir,"sc1csr_gexpr.h5ad"))
SaveH5Seurat(object, filename = paste0(shiny.dir,"sc1csr_gexpr.h5Seurat"))
Convert(paste0(shiny.dir,"sc1csr_gexpr.h5Seurat", dest = "h5ad"))
file.remove(paste0(shiny.dir,"sc1csr_gexpr.h5Seurat"))
