library(ggplot2)
library(Seurat)
library(dplyr)  
library(magrittr)
library(Matrix)
library(qlcMatrix)
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat3_functions.R")

# SLURM_ARRAY_TASK_ID
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
if (length(slurm_arrayid)!=1)  stop("Exact one argument must be supplied!")
# coerce the value to an integer
args <- as.numeric(slurm_arrayid)
print(paste0("slurm_arrayid=",args))

methods = c("T20","IE", "DS", "3C", "T20_non-modified","IE_non-modified", "DS_non-modified", "3C_non-modified")
(method <- methods[args])
Rshiny_path <- paste0("Rshiny/Lung_28_",method,"_umap/")

samples <- c("All_samples","UNC-44-D","UNC-44-T","CU-12-D",
             "CU-12-D-R","CU-12-T","UNC-48-P","UNC-48-D",
             "UNC-48-T","UNC-51-D","UNC-52-D","UNC-54-D",
             "UNC-55-P","UNC-55-D","UNC-55-T","UNC-57-D",
             "UNC-61-D","UNC-66-P","UNC-66-D","UNC-66-T",
             "UNC-67-D","UNC-69-P","UNC-69-D","UNC-69-T",
             "UNC-70-D","UNC-71-P","UNC-71-D","UNC-71-T",
             "VU19-D")

#Prepare Rshiny from single Seurat objects
PrepareShiny <- function(object, samples, Rshiny_path, split.by = "orig.ident",reduction = "tsne",
                         verbose = F,scale =NULL, assay = NULL){
        if(missing(object) | class(object) != "Seurat") stop("samples is not provided")
        if(missing(samples)) stop("samples is not provided")
        if(missing(Rshiny_path)) stop("Rshiny_path is not provided")
        assay = DefaultAssay(object) %||% assay
        Idents(object) <-  split.by
        avaible_samples <- samples %in% c("All_samples",object@meta.data[,split.by])
        if (!all(avaible_samples))
                stop(paste(paste(samples[!avaible_samples],collapse = " "),
                           "are not exist in the data."))
        max_exp <- list()
        if("All_samples" %in% samples) {
                single_object <- object
        } else single_object <- subset(object, idents = samples)
        max_exp = rowMax(single_object[[assay]]@data) %>% as.vector()
        max_exp = max_exp/log(2)
        names(max_exp) = rownames(single_object)
        
        exp <- list()
        tsne <- list()
        for (i in seq_along(samples)){
                sample <- samples[i]
                if(sample == "All_samples") {
                        single_object <- object
                } else single_object <- subset(object, idents = sample)
                #============== exp csv===============
                data <- GetAssayData(single_object)
                data <- as(data, "sparseMatrix")
                data = data/log(2)
                #bad <- rowMax(data) == 0
                #data = data[!bad,]
                if(!is.null(scale)){
                        range <- rowMax(data) - rowMin(data) # range <- apply(data,1,max) - apply(data,1,min)
                        data = sweep(data, 1, range,"/")*scale
                }
                
                if(verbose) {
                        print(sample)
                        print(format(object.size(data),units="MB"))
                }
                exp[[i]] = data
                #============== tsne csv===============
                tsne[[i]] = Embeddings(single_object, reduction = reduction)
                
                svMisc::progress(i/length(samples)*100)
        }
        names(exp) = samples
        names(tsne) = samples
        shiny_data_path <- paste0(Rshiny_path, "data/")
        if(!dir.exists(shiny_data_path)) dir.create(shiny_data_path, recursive = T)
        save(exp,tsne,max_exp, file = paste0(shiny_data_path,basename(Rshiny_path),".Rda"))
}

#============== expression Rda ===============
object = readRDS(file = paste0("data/Lung_28_",args,"-",method,"-harmony_2020331.rds"))

table(object$orig.ident)
DefaultAssay(object) = "SCT"

PrepareShiny(object, samples, Rshiny_path, split.by = "orig.ident", 
             reduction = "umap",verbose = T)