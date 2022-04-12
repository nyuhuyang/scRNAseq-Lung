#====== 3.1 Create Singler Object  ==========================================
# conda activate r4.0.3 linux
invisible(lapply(c("Seurat","SingleR","SingleCellExperiment","magrittr","data.table","Matrix"), function(x) {
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
reference = c("krasnowLung","Lung30")[args]
dataset = "Lung_SCT_62_20220322.rds"
object = readRDS(file = paste0("data/",dataset))
sce <- SingleCellExperiment(list(logcounts=object[["SCT"]]@data),
                                colData=DataFrame(object@meta.data))
rm(object);GC()

# ======= load azimuth Lung data ==============================
if(reference == "krasnowLung"){
    path = "../seurat_resources/azimuth/human_lung/"
    meta.data = fread(paste0(path,"krasnow_hlca_10x_metadata.csv.gz")) %>% 
        as.data.frame() %>% tibble::column_to_rownames("V1")
    counts = fread(paste0(path,"krasnow_hlca_10x_UMIs.csv.gz")) %>% 
        as.data.frame() %>% tibble::column_to_rownames("V1") %>% 
        as.matrix
    format(object.size(counts),unit = "GB")
    counts %<>% as.sparse
    table(rownames(meta.data) == colnames(counts))
    
    options = c("sc","bulk")[1]
    if(options == "sc"){
        libsizes <- colSums(counts)
        size.factors <- libsizes/mean(libsizes)
        Lung <- SingleCellExperiment(list(logcounts=log1p(t(t(counts)/size.factors))),
                                     colData=DataFrame(meta.data))
        rm(counts,meta.data,libsizes,size.factors);GC()
    }
    
    if(options == "bulk"){
        object = CreateSeuratObject(counts,min.cells = 3,names.delim = "-",min.features = 3,
                                    meta.data = meta.data)
        object %<>% NormalizeData()
        exp = AverageExpression(object,group.by = "free_annotation")
        meta.data = meta.data[!duplicated(meta.data$free_annotation),]
        rownames(meta.data) = meta.data$free_annotation
        Lung <- SingleCellExperiment(list(logcounts=exp$RNA),
                                     colData=DataFrame(meta.data))
        rm(counts,meta.data,object,exp);GC()
    }
    
    # ====== conbime data =============
    
    common <- Reduce(intersect, list(rownames(sce),
                                     rownames(Lung)
    ))
    length(common)
    table(Lung$free_annotation)
    system.time(trained <- trainSingleR(ref = Lung[common,],
                                        labels=Lung$free_annotation))
    system.time(pred <- classifySingleR(sce[common,], trained))
    # elapsed 4872.846 sec
    file_name = sub("20220322","20220322_singleR_krasnowLung",dataset)
    saveRDS(object = pred, file = paste0("output/",file_name))
}


if(reference == "Lung30"){
    object = readRDS(file = "data/Lung_SCT_30_20210831.rds")
    Lung <- SingleCellExperiment(list(logcounts=object[["SCT"]]@data),
                                colData=DataFrame(object@meta.data))
    rm(object);GC()
    
    common <- Reduce(intersect, list(rownames(sce),
                                     rownames(Lung)
    ))
    length(common)
    table(Lung$Cell_subtype)
    system.time(trained <- trainSingleR(ref = Lung[common,],
                                        labels=Lung$Cell_subtype))
    system.time(pred <- classifySingleR(sce[common,], trained))
    # elapsed 4872.846 sec
    file_name = sub("20220322","20220322_singleR_Lung30",dataset)
    saveRDS(object = pred, file = paste0("output/",file_name))
}
