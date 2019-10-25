invisible(lapply(c("SingleR"), function(x) {
                           suppressPackageStartupMessages(library(x,character.only = T))
                   }))
source("../R/SingleR_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
#====== 3.1 Create Singler Object  ==========================================
(load(file="data/ref_Lung_24_20191024.Rda"))
(load(file="data/Lung_16_distal_20191022.Rda"))
object_data = object@assays$RNA@data
remove(object);GC()
singler <- CreateBigSingleRObject.1(as.matrix(object_data), annot = NULL, project.name="Lung_distal_COPD",
                                    N = 5000, min.genes = 200, technology = "10X",
                                    species = "Human", citation = "", ref.list = list(ref),
                                    normalize.gene.length = F, variable.genes = "de", fine.tune = F,
                                    reduce.file.size = F, do.signatures = F,
                                    temp.dir = getwd(), numCores = SingleR.numCores)
save(singler,file="output/singlerF_Lung_16_distal_20191023.Rda")
