####################################
invisible(lapply(c("Seurat","dplyr","magrittr","tidyr","SCINA"), function(x) {
                           suppressPackageStartupMessages(library(x,character.only = T))
                   }))
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_functions.R")
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/SingleR_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)

# load Seurat object
object = readRDS(file = "data/Lung_SCT_30_20210831.rds")
# load Enrichr markers
signatures <- fgsea::gmtPathways("../seurat_resources/azimuth/Azimuth_Cell_Types_2021.txt")
lapply(signatures,length)
signatures %<>% pbapply::pblapply(function(x) x[x %in% rownames(object)])

# Predict cell types with SCINA
system.time(results <- BigSCINA(exp = object@assays$SCT@data, signatures,
                                N = 5000, max_iter = 100, convergence_n = 10, 
                convergence_rate = 0.999, sensitivity_cutoff = 0.9, rm_overlap = 1, allow_unknown = 1, 
                log_file=paste0(path,'SCINA.log')))
table(results$cell_labels)
saveRDS(results,paste0(path,"SCINA_Lung30_Azimuth_Cell_Types_2021.rds"))

saveRDS(object@meta.data ,"output/20211222/meta.data_SCINA_Lung30_Azimuth_Cell_Types_2021.rds")

SCINA_results <- readRDS("output/20211220/SCINA_Lung30_Azimuth_Cell_Types_2021.rds")
SCINA_results$SCINA_Azimuth_cell_labels = gsub(" CL0.*","",SCINA_results$cell_labels)
SCINA_results$SCINA_Azimuth_cell_labels = gsub(" UBER.*","",SCINA_results$SCINA_Azimuth_cell_labels)
