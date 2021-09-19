########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
#
# ######################################################################
# conda activate r4.0.3
#devtools::install_github("immunogenomics/harmony", ref= "ee0877a",force = T)
invisible(lapply(c("Seurat","dplyr","ggplot2","cowplot","pbapply","sctransform","harmony","magrittr"), function(x) {
        suppressPackageStartupMessages(library(x,character.only = T))
}))
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)


#########################

# ====== load old Lung data=============
Lung = readRDS(file = "data/Lung_SCT_30_20200710.rds")

# ====== load new Lung data=============
object = readRDS(file = "data/Lung_59_20210814.rds")

Lung$barcode = gsub(".*_","",colnames(Lung)) %>% gsub("-1$","",.)

Lung$orig.ident = gsub("-","_",Lung$orig.ident)
Lung$orig.ident %<>% gsub("VU19_D","VU_19_D",.)
table(Lung$orig.ident %in% object$orig.ident)
Lung$barcode = paste0(Lung$orig.ident,"-",Lung$barcode)
Lung %<>% RenameCells(new.names = Lung$barcode)

table(Lung$barcode %in% colnames(object))

meta_data = Lung@meta.data[,c("barcode","cell_types","cell_types.colors","Doublets")]
rownames(meta_data) = meta_data$barcode
meta_data = meta_data[Lung$barcode %in% colnames(object),]
object$barcode = colnames(object)
table(meta_data$barcode %in% object$barcode)
meta.data = object@meta.data
meta.data$cell_types = NULL
meta.data$cell_types.colors = NULL
meta.data %<>% left_join(meta_data,by = "barcode")
meta.data[is.na(meta.data$cell_types),"cell_types"] = "Unknown"

rownames(meta.data) = meta.data$barcode
table(rownames(meta.data) == rownames(object@meta.data))

object@meta.data = meta.data
object %<>% subset(subset = cell_types != "Unknown")
UMAPPlot.1(object,do.print= T, pt.size = 0.5,raster=FALSE, group.by = "cell_types")
remove(Lung);GC()

object[["new.umap"]] <- CreateDimReducObject(embeddings = object@reductions[["umap"]]@cell.embeddings,
                                             key = "newUMAP_", assay = DefaultAssay(object))

object[["new.tsne"]] <- CreateDimReducObject(embeddings = object@reductions[["tsne"]]@cell.embeddings,
                                             key = "newtSNE_", assay = DefaultAssay(object))

table(colnames(Lung) %in% colnames(object))
table(colnames(object) %in% colnames(Lung))
Lung %<>% subset(cells = colnames(object))

object[["hg19.umap"]] <- CreateDimReducObject(embeddings = Lung@reductions[["umap"]]@cell.embeddings,
                                             key = "hg19UMAP_", assay = DefaultAssay(Lung))

object@reductions$pca = NULL
object@reductions$new.umap = NULL
object@reductions$new.tsne = NULL

#=====================
object %<>% SCTransform(method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = TRUE)
object %<>% FindVariableFeatures(selection.method = "vst",
                                 num.bin = 20,
                                 mean.cutoff = c(0.1, 8), dispersion.cutoff = c(1, Inf))

object %<>% ScaleData
npcs <- 105
object %<>% RunPCA(verbose = T,npcs = npcs, features = VariableFeatures(object))

jpeg(paste0(path,"S1_ElbowPlot.jpeg"), units="in", width=10, height=7,res=600)
ElbowPlot(object, ndims = npcs)
dev.off()

object %<>% FindNeighbors(reduction = "pca",dims = 1:npcs)
object %<>% FindClusters(resolution = 0.8)
object %<>% RunUMAP(reduction = "pca", dims = 1:npcs)

meta.data = object@meta.data
meta.data = meta.data[!duplicated(meta.data$cell_types),]
meta.data = meta.data[order(meta.data$cell_types),]

resolutions = c(seq(0.1,5, by = 0.1),seq(2,5, by = 1))
for(i in 1:length(resolutions)){
        object %<>% FindClusters(resolution = resolutions[i])
        Progress(i,length(resolutions))
}

UMAPPlot.1(object,do.print= T, pt.size = 0.5,raster=FALSE,cols=meta.data$cell_types.colors, group.by = "cell_types")
saveRDS(object, file = "data/Lung_30_20210831.rds")

#=======1.9 save SCT only =======================================
format(object.size(object),unit = "GB")

format(object.size(object@assays$RNA),unit = "GB")
format(object.size(object@assays$integrated),unit = "GB")
object[['RNA']] <- NULL
object[['integrated']] <- NULL
object[["SCT"]]@scale.data = matrix(0,0,0)
format(object.size(object),unit = "GB")
saveRDS(object, file = "data/Lung_SCT_30_20210831.rds")

#======== prepare shiny =============================
library(Seurat)
library(ShinyCell)
library(magrittr)
library(SeuratData)
library(SeuratDisk)
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_functions.R")

object = readRDS(file = "data/Lung_SCT_30_20210831.rds")
#meta.data = object@meta.data
#meta.data = meta.data[, grep("SCT_snn_res",colnames(meta.data), invert = TRUE)]

#object %<>% FindNeighbors(reduction = "umap",dims = 1:2)
#resolutions = c(0.8,1:5)
#for(i in 1:length(resolutions)){
#        object %<>% FindClusters(resolution = resolutions[i], algorithm = 1)
#        Progress(i,length(resolutions))
#}
#colnames(meta.data) %<>% gsub("SCT_snn_res.","0UMAP_res=",.)
#meta_data = readRDS(file = "output/20210901/meta.data_Cell_subtype.rds")
#table(rownames(meta.data) == rownames(meta_data))
#colnames(meta.data) %<>% gsub("cell_types","old_cell_types",.)
#meta.data %<>% cbind(meta_data)
#object@meta.data = meta.data

#test_df = data.frame(min_dist = rep(2:5/10,each = 5),
#                     spread = rep(3:7/5,times = 4))

#for(i in 20) {
#        spread <- test_df[i,"spread"]
#        min.dist <- test_df[i,"min_dist"]
#        file.name = paste0("dist.",min.dist,"_spread.",spread)
#        umap = readRDS(paste0("output/20210901/","umap_",file.name,".rds"))
#        object[[paste0("umap_",file.name)]] <- CreateDimReducObject(embeddings = umap[[paste0("umap_",file.name)]]@cell.embeddings,
#                                                                       key = paste0(i,"UMAP_"), assay = DefaultAssay(object))
#        meta.data = readRDS(paste0("output/20210901/meta.data_",file.name,".rds"))
#        colnames(meta.data) %<>% gsub(paste0(file.name,"."),paste0(i,"UMAP_res="),.)
#        
#        object@meta.data %<>% cbind(meta.data)
#        Progress(i,nrow(test_df))
#}

object$Regions %<>% factor(levels = c("proximal","distal","terminal","COPD"))
meta_data = object@meta.data
meta_data = meta_data[!duplicated(meta_data$orig.ident),]
meta_data %<>% group_by(Regions) %>% arrange(Regions, orig.ident)
object$orig.ident %<>% factor(levels = meta_data$orig.ident)
meta.to.include =c("Cell_subtype","Cell_type","UMAP_land","Family","Superfamily",
                   "old_cell_types",
                   "orig.ident","Regions","Patient",
                   "nCount_SCT","nFeature_SCT","percent.mt",
                   grep("UMAP_res=",colnames(object@meta.data), value = T),
                   "Doublets")
table(meta.to.include %in% colnames(object@meta.data))
scConf = createConfig(object, meta.to.include =meta.to.include,maxLevels = 300)
makeShinyApp(object, scConf, gex.assay = "SCT",gene.mapping = TRUE,
             gex.slot = "data",default.gene1 = "SCGB1A",default.gene2 = "KRT15",
             default.multigene = c( "KRT15","SERPINB3","MUC5AC","SCGB1A1","SCGB3A2","SFTPB","FOXA2"),
             default.dimred = c("UMAP_1","UMAP_2"),shiny.dir = "shinyApp/Lung_30_hg38/",
             shiny.title = "Region-specific Lung scRNA-seq")

max_exp = qlcMatrix::rowMax(object@assays[["SCT"]]@data)
max_exp_df = data.frame("val"= as.vector(max_exp),row.names = rownames(object))
saveRDS(max_exp_df,"shinyApp/Lung_30_hg38/sc1maxlvl.rds")


meta.data = object@meta.data
sc1conf = readRDS("shinyApp/Lung_30_hg38/sc1conf.rds")
meta.data1 = meta.data[!duplicated(meta.data$Cell_subtype),]
meta.data1 = meta.data1[order(meta.data1$Cell_subtype),]
sc1conf$fCL[1] = paste(meta.data1$Cell_subtype.colors,collapse = "|")

meta.data2 = meta.data[!duplicated(meta.data$old_cell_types),]
meta.data2 = meta.data2[order(meta.data2$old_cell_types),]
sc1conf$fCL[6] = paste(meta.data2$old_cell_types.colors,collapse = "|")

sc1conf$fCL[8] = paste(c("#1F78B4","#4ca64c","#E6AB02","#FF0000"),collapse = "|")
sc1conf$fCL[151] = "orange|black"
saveRDS(sc1conf,"shinyApp/Lung_30_hg38/sc1conf.rds")

sc1def = readRDS("shinyApp/Lung_30_hg38/sc1def.rds")
sc1def$grp1 = "Cell_subtype"
sc1def$meta1 = "Cell_subtype"
saveRDS(sc1def,"shinyApp/Lung_30_hg38/sc1def.rds")



object@meta.data = object[["cell_types"]]
object@reductions = NULL
file.remove("shinyApp/Lung_30_hg38/sc1csr_gexpr.h5ad")
SaveH5Seurat(Lung, filename = "shinyApp/Lung_30_hg38/sc1csr_gexpr.h5Seurat")
Convert("shinyApp/Lung_30_hg38/sc1csr_gexpr.h5Seurat", dest = "h5ad")
file.remove("shinyApp/Lung_30_hg38/sc1csr_gexpr.h5Seurat")
