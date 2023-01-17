########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
#
# ######################################################################
# conda activate r4.1.3
#devtools::install_github("immunogenomics/harmony", ref= "ee0877a",force = T)
invisible(lapply(c("Seurat","dplyr","ggplot2","cowplot","pbapply","sctransform","harmony","magrittr"), function(x) {
    suppressPackageStartupMessages(library(x,character.only = T))
}))
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)


########################################################################
#
#  1 Seurat Alignment
#
# ######################################################################
#======1.1 Setup the Seurat objects =========================
# read sample summary list
# read sample summary list
meta.data <- readxl::read_excel("data/RNA-seq/ALI bulk RNA-seq long-term.xlsx",sheet = "metadata") %>%
    as.data.frame()
meta.data$Sex %<>% plyr::mapvalues(from = c("M","F"),
                                   to = c("Male","Female"))
object$Regions %<>% gsub("^Prox$","Proximal",.)

meta.data$Sex %<>% factor(levels = c("Male","Female"))
meta.data$Day %<>% factor(levels = paste0("Day",c(0,3,7,14,18,21,28,56,112,137,224,366)))
table(meta.data$Day)
meta.data$orig.ident <- paste0(meta.data$Donor,"_",meta.data$Regions,"_",meta.data$Day)
table(duplicated(meta.data$orig.ident))
meta.data$orig.ident %<>% make.unique()
rownames(meta.data) <- meta.data$orig.ident
nrow(meta.data)
#======1.2 load  Seurat =========================
FPKM <- readxl::read_excel("data/RNA-seq/ALI bulk RNA-seq long-term.xlsx",sheet = "MergedFPKMs")
FPKM$gene %<>% make.unique()
FPKM %<>% tibble::column_to_rownames("gene")
colnames(FPKM) <- meta.data$orig.ident
FPKM %<>% data.matrix

object <- CreateSeuratObject(FPKM,min.cells = 0,names.delim = "-",min.features = 0,meta.data = meta.data)

table(meta.data$orig.ident %in% object$orig.ident)

#https://github.com/satijalab/seurat/issues/3805
# remove NA from raw data
RowsNA<-names(which(rowSums(is.na(object@assays$RNA@counts))>0))
'%!in%' <- function(x,y)!('%in%'(x,y)) #this is a NOT IN function
RowsKEEP<-rownames(object)[rownames(object) %!in% RowsNA]
object <- subset(object,features=RowsKEEP)

# Determine the ‘dimensionality’ of the dataset  =========
npcs = 50

DefaultAssay(object) <- "RNA"
object %<>% NormalizeData()
object <- FindVariableFeatures(object = object, selection.method = "vst",
                               num.bin = 20, nfeatures = 3000,
                               mean.cutoff = c(0.1, 8), dispersion.cutoff = c(1, Inf))
object %<>% ScaleData(verbose = FALSE)
object %<>% RunPCA(npcs = npcs, verbose = FALSE)


jpeg(paste0(path,"ElbowPlot_RNA.jpeg"), units="in", width=10, height=7,res=600)
print(ElbowPlot(object,ndims = npcs))
dev.off()
object <- JackStraw(object, num.replicate = 20,dims = npcs)
object <- ScoreJackStraw(object, dims = 1:npcs)

for(i in 0:5){
    a = i*10+1; b = (i+1)*10
    jpeg(paste0(path,"JackStrawPlot_",a,"_",b,".jpeg"), units="in", width=10, height=7,res=600)
    print(JackStrawPlot(object, dims = a:b))
    dev.off()
    Progress(i, 5)
}
p.values = object[["pca"]]@jackstraw@overall.p.values
print(npcs <- max(which(p.values[,"Score"] <=0.05)))

npcs = 12

#======1.6 Performing SCTransform  =========================
format(object.size(object),unit = "GB")
options(future.globals.maxSize= object.size(object)*10)
# read and select mitochondial genes
mito <- "^MT-"
message("mito.genes:")

(mito.features <- grep(pattern = mito, x = rownames(object), value = TRUE))
object[["percent.mt"]] <- PercentageFeatureSet(object = object, pattern = mito)

object %<>% SCTransform(method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = TRUE)
object[['RNA']]@scale.data = matrix(0,0,0)
object[['SCT']]@scale.data = matrix(0,0,0)

DefaultAssay(object)
object <- FindVariableFeatures(object = object, selection.method = "vst",
                               num.bin = 20, nfeatures = 3000,
                               mean.cutoff = c(0.1, 8), dispersion.cutoff = c(1, Inf))


s.genes <- cc.genes$s.genes %>% gsub("MLF1IP","CENPU",.)
g2m.genes <- cc.genes$g2m.genes %>% plyr::mapvalues(from = c("FAM64A", "HN1"),
                                                    to = c("PIMREG","JPT1"))
object <- CellCycleScoring(object, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
colnames(object@meta.data) %<>% sub("Phase","cell.cycle.phase",.)


object %<>% ScaleData(verbose = FALSE)
object %<>% RunPCA(verbose = T,npcs = npcs)
jpeg(paste0(path,"S1_ElbowPlot_SCT.jpeg"), units="in", width=10, height=7,res=600)
ElbowPlot(object, ndims = npcs)
dev.off()


#======1.8 UMAP from harmony =========================

object %<>% RunUMAP(reduction = "pca", dims = 1:npcs,return.model = TRUE)
object %<>% FindNeighbors(reduction = "pca",dims = 1:npcs)

object %<>% FindClusters(resolution = 0.8)
#system.time(object %<>% RunTSNE(reduction = "pca", dims = 1:npcs))

saveRDS(object, file = "data/Lung_bulk_time_20230112.rds")

object <- readRDS(file = "data/Lung_bulk_time_20230112.rds")

npcs = 12
object %<>% FindNeighbors(reduction = "pca",dims = 1:npcs)
resolutions = c(seq(0.9,1.5, by = 0.1),seq(2,5, by = 1))
for(i in 1:length(resolutions)){
    object %<>% FindClusters(resolution = resolutions[i], algorithm = 1,verbose = F)
    Progress(i, length(resolutions))
}
object$orig.ident <- colnames(object)
object$Sample <- gsub("Proximal|Prox","P",object$orig.ident) %>%
                gsub("Distal","D",.) %>%
                gsub("UNC","UNC-",.) %>%
                gsub("Day","d",.) %>%
                gsub("_","-",.)
saveRDS(object, file = "data/Lung_bulk_time_20230112.rds")
