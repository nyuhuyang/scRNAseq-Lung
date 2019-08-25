library(SingleR)
source("../R/SingleR_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
#====== 3.1 Create Singler Object  ==========================================
(load(file = "data/Lung.data_24_20190824.Rda"))
attach(blueprint_encode)
singler <- CreateBigSingleRObject.1(object_data, annot = NULL, project.name="EC-SR-5444~5989",
                                    N = 5000, min.genes = 50, technology = "10X",
                                    species = "Human", citation = "", 
                                    ref.list = list(blueprint_encode),
                                    normalize.gene.length = F, variable.genes = "de", 
                                    fine.tune = F,
                                    reduce.file.size = F, do.signatures = F, do.main.types = T,
                                    temp.dir = getwd(), numCores = SingleR.numCores)

save(singler,file="output/singlerF_Lung_24_20190824.Rda")
