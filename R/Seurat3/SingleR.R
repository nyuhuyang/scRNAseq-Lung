library(SingleR)
source("../R/SingleR_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
#====== 3.1 Create Singler Object  ==========================================
(load(file = "data/Lung.data_harmony_12_20190614.Rda"))
(load(file = "../SingleR/data/Blueprint_encode.RData"))
singler = CreateSinglerSeuratObject(object_data, annot = NULL, project.name="ddSeq_5444_5517",
                                    min.genes = 500, technology = "10X", species = "Human" ,
                                    ref.list = list(Blueprint_encode),fine.tune = T,
                                    normalize.gene.length = F, min.cells = 2, npca = 10,
                                    regress.out = "nUMI", reduce.seurat.object = T)
singler <- CreateBigSingleRObject(object_data, annot = NULL, project.name="ddSeq_5444_5517",
                                    N = 5000, min.genes = 500, technology = "10X",
                                    species = "Human", citation = "", 
                                    ref.list = list(Blueprint_encode),
                                    normalize.gene.length = F, variable.genes = "de", fine.tune = F,
                                    reduce.file.size = F, do.signatures = F, do.main.types = T,
                                    temp.dir = getwd(), numCores = SingleR.numCores)

save(singler,file="output/singlerF_Lung_12_20190614.Rda")
