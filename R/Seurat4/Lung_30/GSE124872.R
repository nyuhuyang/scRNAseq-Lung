invisible(lapply(c("Seurat","dplyr","ggplot2","scater","magrittr","pbapply",
                   "cowplot"), function(x) {
                       suppressPackageStartupMessages(library(x,character.only = T))
                   }))
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_functions.R")

# Need 16GB

#---------
datasets <- c("WholeLung","HighResolution")[args]
read_path <- "data/annotation_references/GSE124872/"

"@shell
gunzip GSE124872_raw_counts_single_cell.RData.gz
"

(load(paste0(read_path,"GSE124872_raw_counts_single_cell.RData")))
#line 34 in Read10X data <- readMM(file = matrix.loc); data <- t(data) # for HighResolution
meta.data <- data.table::fread(paste0(read_path,"GSE124872_Angelidis_2018_metadata.csv.gz")) %>%
    as.data.frame() %>% tibble::column_to_rownames("V1")
meta.data$barcode <- gsub(".*:","",rownames(meta.data))
barcode <- gsub(":.*","",colnames(raw_counts))

barcode_identifider <- data.frame(c("muc3838","old_1"),
                                  c("muc3839","old_2"),
                                  c("muc3840","young_1"),
                                  c("muc3841","young_2"),
                                  c("muc4166","old_3"),
                                  c("muc4167","old_4"),
                                  c("muc4168","old_5"),
                                  c("muc4169","young_3"),
                                  c("muc4170","young_4"),
                                  c("muc4172","young_5"),
                                  c("muc4173","young_6"),
                                  c("muc4174","old_6"),
                                  c("muc4175","old_7"),
                                  c("muc4654","young_7"),
                                  c("muc4657","young_8")) %>% t
meta.data$name <- plyr::mapvalues(meta.data$identifier,
                                  from = barcode_identifider[,1],
                                  to = barcode_identifider[,2])
rownames(meta.data) <- paste0(meta.data$name,":",meta.data$barcode)
table(colnames(raw_counts) == rownames(meta.data))
object <- CreateSeuratObject(counts = raw_counts,meta.data = meta.data)

object %<>% NormalizeData(verbose = TRUE) %>%
    FindVariableFeatures(verbose = TRUE) %>%
    ScaleData(verbose = TRUE) %>%
    RunPCA(verbose = TRUE,npcs = 50) %>%
    RunUMAP(verbose = TRUE, dims = 1:50) %>%
    FindNeighbors(reduction = "pca",dims = 1:50) %>%
    FindClusters(resolution = 0.8, algorithm = 1,verbose = F)


UMAPPlot(object, group.by = "celltype")
saveRDS(object, file = paste0("data/GSE124872_single_cell.rds"))
#===SingleR============================
invisible(lapply(c("SingleR","SingleCellExperiment","magrittr","data.table","Matrix"), function(x) {
    suppressPackageStartupMessages(library(x,character.only = T))
}))
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/SingleR_functions.R")

sce <- SingleCellExperiment(list(logcounts=object[["RNA"]]@data),
                            colData=DataFrame(object@meta.data))

Lung = readRDS(file = "data/Lung_SCT_30_20210831.rds")
Lung$Cell_subtype_orig_ident <- paste0(Lung$Cell_subtype,"_",Lung$orig.ident)
Cell_subtype_orig_ident_exp <- AverageExpression(Lung,assays = "SCT",group.by = "Cell_subtype_orig_ident")
sce_Lung <- SingleCellExperiment(list(logcounts=Cell_subtype_orig_ident_exp$SCT),
                                 colData=DataFrame("Cell_subtype" = gsub("_.*","",colnames(Cell_subtype_orig_ident_exp$SCT)),
                                                   row.names = colnames(Cell_subtype_orig_ident_exp$SCT)))

rownames(sce) %<>% toupper()
common <- Reduce(intersect, list(rownames(sce),
                                 rownames(sce_Lung)
))
length(common)
table(sce_Lung$Cell_subtype)
system.time(trained <- trainSingleR(ref = sce_Lung[common,],
                                    labels=sce_Lung$Cell_subtype))
system.time(pred <- classifySingleR(sce[common,], trained))
# elapsed 4872.846 sec
saveRDS(object = pred, file = "output/GSE124872_scLung_singelR_bulkRNA-Lung30.rds")


object <- readRDS(paste0("data/GSE124872_single_cell.rds"))
pred <- readRDS("output/GSE124872_scLung_singelR_bulkRNA-Lung30.rds")

singlerDF = data.frame("Cell_subtype" = pred$pruned.labels,
                       row.names = rownames(pred))
table(is.na(singlerDF$Cell_subtype))
singlerDF$Cell_subtype[is.na(singlerDF$Cell_subtype)]= "unknown"
table(singlerDF$Cell_subtype)
if(all(colnames(object) == rownames(singlerDF))){
    print("all cellID match!")
    object$Cell_subtype = singlerDF$Cell_subtype
}
#Within the mouse CELL TYPE categories club cells and goblet cells:
#identify those labeled as TASC CELL SUBTYPE (group A) using your method
TASC <- object$celltype %in% c("Club_cells","Goblet_cells") & object$Cell_subtype %in% "TASC"
table(TASC)
object$celltype[TASC] <-"TASC"
#-------
library(Seurat)
library(ShinyCell)
library(magrittr)
library(SeuratData)
library(SeuratDisk)
library(dplyr)
object$orig.ident <- object$identifier
object$mouse <- object$identifier
object$age <- gsub("_[1-9]","",object$name)
meta.data <- object@meta.data
meta.to.include =c("celltype","Cell_subtype","age","orig.ident","nCount_RNA","nFeature_RNA","nGene","nUMI",
                   "identifier","res.2","name","grouping",
                   "batch","cells","cluster",  
                   "RNA_snn_res.0.8","mouse")


shiny.dir = paste0("shinyApp/GSE124872_Lung_single_cells/")
table(meta.to.include %in% colnames(object@meta.data))
meta.to.include[!meta.to.include %in% colnames(object@meta.data) ]
scConf = createConfig(object, meta.to.include =meta.to.include,maxLevels = 300)
makeShinyApp(object, scConf, gex.assay = "RNA",gene.mapping = TRUE,
             gex.slot = "data",default.gene1 = "Scgb1a1",default.gene2 = "Krt15",
             default.multigene = c("Krt15","Muc5ac","Scgb1a1","Scgb3a2","Sftpb","Foxa2" ),
             default.dimred = c("UMAP_1","UMAP_2"),shiny.dir = shiny.dir,
             shiny.title = "GSE124872 Lung single cells"
             )

max_exp = qlcMatrix::rowMax(object@assays[["RNA"]]@data)
max_exp_df = data.frame("val"= as.vector(max_exp),row.names = rownames(object))
saveRDS(max_exp_df,paste0(shiny.dir,"sc1maxlvl.rds"))

sc1def = readRDS(paste0(shiny.dir,"sc1def.rds"))
sc1def$grp1 = "celltype"
sc1def$grp2 = "Cell_subtype"
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
Convert(paste0(shiny.dir,"sc1csr_gexpr.h5Seurat"), dest = "h5ad")
file.remove(paste0(shiny.dir,"sc1csr_gexpr.h5Seurat"))
