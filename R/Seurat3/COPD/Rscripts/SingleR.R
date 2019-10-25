invisible(lapply(c("Seurat"), function(x) {
                           suppressPackageStartupMessages(library(x,character.only = T))
                   }))
source("../R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
#====== 3.1 Create Singler Object  ==========================================
(load(file="output/singlerF_Lung_16_distal_20191023.Rda"))

(load(file="data/Lung_16_distal_20191022.Rda"))
DefaultAssay(object) = "RNA"
sce <- as.SingleCellExperiment(object)
remove(object);GC()

# for new SingleR ======================
system.time(singler <- SingleR(test=sce, ref=Lung_sce, labels=Lung_sce$labels))

classifyBigSingleR <- function(test, trained, N=5000, quantile = 0.8, fine.tune = TRUE,
                               tune.thresh = 0.05, sd.thresh = NULL, prune = TRUE,
                               assay.type = "logcounts", check.missing = TRUE,
                               BPPARAM = MulticoreParam()){
        singler_list <- list()
        n = ncol(test)
        s = seq(1,n,by=N)
        for (i in s) {
                print(i)
                A = seq(i,min(i+N-1,n))
                tryCatch({
                        singler_list[[i]] <- classifySingleR(test = test[,A], trained = trained, 
                                                     quantile = quantile, fine.tune = fine.tune,
                                                     tune.thresh = tune.thresh, sd.thresh = sd.thresh, 
                                                     prune = prune,assay.type = assay.type, 
                                                     check.missing = check.missing,
                                                     BPPARAM = BPPARAM)
                }, error=function(e) invisible(e))
                GC()
        }
        return(singler_list)
}

#system.time(singler_list <- classifyBigSingleR(sce, trained, N = 5000,
#                                               BPPARAM = MulticoreParam(5, stop.on.error = TRUE)))
system.time(singler <- classifySingleR(sce, trained))
save(singler,file="output/singlerT_Lung_16_distal_20191024.Rda")
