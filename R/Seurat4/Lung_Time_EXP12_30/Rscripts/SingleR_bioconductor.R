#====== 3.1 Create Singler Object  ==========================================
# conda activate r4.0.3 linux
invisible(lapply(c("Seurat","SingleR","SingleCellExperiment","magrittr","data.table","Matrix","arrow"), function(x) {
    suppressPackageStartupMessages(library(x,character.only = T))
}))
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_functions.R")
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/SingleR_functions.R")

# Need 64GB ?
set.seed(101)
# SLURM_ARRAY_TASK_ID
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
if (length(slurm_arrayid)!=1)  stop("Exact one argument must be supplied!")
# coerce the value to an integer
args <- as.integer(as.character(slurm_arrayid))
print(paste0("slurm_arrayid=",args))

# ====== load single cell =============
reference = c("pseduBulk_human_lung_v2","pseduBulk_Lung30")[args]
object = readRDS(file = "data/Lung_time15_20220523.rds")
sce <- SingleCellExperiment(list(logcounts=object[["SCT"]]@data),
                                colData=DataFrame(object@meta.data))
rm(object);GC()

# ======= load azimuth Lung data ==============================
if(reference == "pseduBulk_human_lung_v2"){
    #counts.path <- "data/annotation_references/human_lung_v2/data/counts.rds"
    #annotations.path <- "data/annotation_references/human_lung_v2/data/annotations.parquet"
    #dr.path <- "data/annotation_references/human_lung_v2/data/scanvi.parquet"
    
    #mtx <- Read10X(counts.path)
    #obj <- CreateSeuratObject(counts = mtx)
    #annotations <- read_parquet(annotations.path)
    #rownames(annotations) <- annotations$X
    #annotations$X <- NULL
    #obj <- AddMetaData(obj, metadata = annotations)
    #obj %<>% NormalizeData()
    obj <- readRDS("data/annotation_references/human_lung_v2/data/human_lung_v2.rds")
    
    obj$ann_finest_level_sample <- paste0(obj$ann_finest_level,"@",obj$sample)
    obj[['RNA']] <- NULL
    obj[["SCT"]]@counts = matrix(0,0,0)
    obj[["SCT"]]@scale.data = matrix(0,0,0)
    obj[["SCT"]]@SCTModel.list = list()
    format(object.size(obj),unit = "GB")
    options(future.globals.maxSize= object.size(obj)*50)
    df <- table(obj$ann_finest_level_sample) %>% as.data.frame
    df <- df[df$Freq > 10,]
    ann_finest_level <- unique(obj$ann_finest_level)
    sub_ann_finest_level <- unique(gsub("@.*","",df$Var1))
    table(ann_finest_level %in% sub_ann_finest_level)
    obj %<>% subset(subset = ann_finest_level_sample %in% df$Var1)
    exp = AverageExpression(obj,group.by = "ann_finest_level_sample",assays = "SCT")
    
    exp = log1p(exp$SCT)
    meta.data <- obj@meta.data
    meta.data = meta.data[!duplicated(meta.data$ann_finest_level_sample),]
    rownames(meta.data) = meta.data$ann_finest_level_sample
    Lung_v2 <- SingleCellExperiment(list(logcounts=exp),
                                 colData=DataFrame("ann_finest_level"=gsub("@.*","",colnames(exp)),
                                                   row.names = colnames(exp)))
    rm(obj,meta.data,exp);GC()

    
    # ====== conbime data =============
    
    common <- Reduce(intersect, list(rownames(sce),
                                     rownames(Lung_v2)
    ))
    length(common)
    table(Lung_v2$ann_finest_level)
    system.time(trained <- trainSingleR(ref = Lung_v2[common,],
                                        labels=Lung_v2$ann_finest_level))
    system.time(pred <- classifySingleR(sce[common,], trained))
    saveRDS(object = pred, file = "output/Lung_time15_20220523_pseduBulk_human_lung_v2_singleR_pred.rds")
}


if(reference == "pseduBulk_Lung30"){
    object = readRDS(file = "data/Lung_SCT_30_20210831.rds")
    meta.data = readRDS("output/Lung_30_20210831_metadata_v2.rds")
    if(all(colnames(object) == rownames(meta.data))){
        print("all cellID match!")
        object@meta.data = meta.data
    }
    object$Cell_subtype_sample <- paste0(object$Cell_subtype,"@",object$orig.ident)
    exp = AverageExpression(object,group.by = "Cell_subtype_sample",assays = "SCT")
    exp = log1p(exp$SCT)
    meta.data = meta.data[!duplicated(meta.data$Cell_subtype_sample),]
    rownames(meta.data) = meta.data$Cell_subtype_sample
    Lung <- SingleCellExperiment(list(logcounts=exp),
                                    colData=DataFrame("Cell_subtype"=gsub("@.*","",colnames(exp)),
                                                      row.names = colnames(exp)))
    
    rm(object,meta.data,exp);GC()
    
    common <- Reduce(intersect, list(rownames(sce),
                                     rownames(Lung)
    ))
    length(common)
    table(Lung$Cell_subtype)
    system.time(trained <- trainSingleR(ref = Lung[common,],
                                        labels=Lung$Cell_subtype))
    system.time(pred <- classifySingleR(sce[common,], trained))
    saveRDS(object = pred, file = "output/Lung_time15_20220523_pseduBulk_pseduBulk_Lung30_singleR_pred.rds")
    
}
