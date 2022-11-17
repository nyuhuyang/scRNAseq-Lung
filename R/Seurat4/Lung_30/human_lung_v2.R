library(Seurat)
library(SeuratDisk)
library(arrow)
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)


"@shell
mkdir data/annotation_references/human_lung_v2
cd data/annotation_references/human_lung_v2
wget https://seurat.nygenome.org/hlca_ref_files/counts.rds -P data
wget https://seurat.nygenome.org/hlca_ref_files/annotations.parquet -P data
wget https://seurat.nygenome.org/hlca_ref_files/scanvi.parquet -P data
wget https://seurat.nygenome.org/hlca_ref_files/ts_opt.rds -P reference
"
#---------
counts.path <- "data/annotation_references/human_lung_v2/data/counts.rds"
annotations.path <- "data/annotation_references/human_lung_v2/data/annotations.parquet"
dr.path <- "data/annotation_references/human_lung_v2/data/scanvi.parquet"

mtx <- Read10X(counts.path)
obj <- CreateSeuratObject(counts = mtx)


# load annotations
library(arrow)
annotations <- read_parquet(annotations.path)
rownames(annotations) <- annotations$X
annotations$X <- NULL
obj <- AddMetaData(obj, metadata = annotations)

# load in the scANVI latent space which contains 30 dimensions
latent.space <- read_parquet(dr.path)
rownames(latent.space) <- Cells(obj)
scanvi.dr <- CreateDimReducObject(embeddings = as.matrix(latent.space), key = "SCANVI")
obj[["scanvi"]] <- scanvi.dr

# find neighbors based on scANVI latent space
obj <- FindNeighbors(obj, reduction = "scanvi")
saveRDS(obj, file = "data/annotation_references/human_lung_v2/data/annotation_references/human_lung_v2.rds")
obj <- readRDS(file = "data/annotation_references/human_lung_v2/data/annotation_references/human_lung_v2.rds")

# Run SCTransform on the raw counts
format(object.size(obj)*27,unit = "GB")
options(future.globals.maxSize= object.size(obj)*27)
rm(mtx,annotations,scanvi.dr,latent.space)
GC();GC();GC()


obj <- SCTransform(obj, method = "glmGamPoi") # need ~500GB memory
format(object.size(obj),unit = "GB")
# run sPCA
obj <- RunSPCA(object = obj, assay = "SCT", graph = "RNA_snn")

# Force RunUMAP to run with n.epochs to prevent RunUMAP from running for 0 epochs
# Related: https://scanpy.discourse.group/t/umap-incorrectly-installed/663
obj <- RunUMAP(obj, dims = 1:30, reduction = "scanvi", n.epochs = 200, return.model = TRUE)
obj[["SCT"]]@scale.data = matrix(0,0,0)
# save the full size object to perform marker identification
saveRDS(obj, file = "data/annotation_references/human_lung_v2/data/annotation_references/human_lung_v2.rds")

plot1 <- UMAPPlot(obj, group.by = "ann_level_1",raster=FALSE,label = T)
save.path <- paste0("output/",gsub("-","",Sys.Date()))
if(!dir.exists(save.path)) dir.create(save.path, recursive = T)


jpeg(paste0(save.path,"/ann_level_1.jpeg"),units="in", width=10, height=7,res=600)
print(plot1)
dev.off()



#################################
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_differential_expression.R")
object = readRDS(file ="data/annotation_references/human_lung_v2/data/annotation_references/human_lung_v2.rds")
sub_object <- subset(object, subset = ann_finest_level %in% c("Transitional Club-AT2","AT2"))
marker1 <- FindMarkers_UMI(sub_object,ident.1 = "Transitional Club-AT2",
                           ident.2 = "AT2",
                           group.by = "ann_finest_level",logfc.threshold = 0.1)
marker1$cluster = "Transitional Club-AT2 vs. AT2"
openxlsx::write.xlsx(marker1, file =  paste0(path,"2022-10-13-ann_finest_level.Secretory.xlsx"),
                     colNames = TRUE,rowNames = F,borders = "surrounding")

library(arrow)
annotations.path <- "data/annotation_references/human_lung_v2/data/annotations.parquet"
annotations <- read_parquet(annotations.path)
Secretory <- unique(annotations[annotations$ann_level_3 %in% c("Secretory"),"ann_finest_level"])
sub_object2 <- subset(object, subset = ann_finest_level %in% Secretory)
marker2 <- FindMarkers_UMI(sub_object2,ident.1 = "Transitional Club-AT2",
                           ident.2 = NULL,
                           group.by = "ann_finest_level",logfc.threshold = 0.1)

openxlsx::write.xlsx(marker2, file =  paste0(path,"2022-10-14-Transitional Club-AT2 vs Club (nasal)_Club (non-nasal)_Goblet (bronchial)_Goblet (nasal)_Goblet (subsegmental) in   .xlsx"),
                     colNames = TRUE,rowNames = F,borders = "surrounding")


