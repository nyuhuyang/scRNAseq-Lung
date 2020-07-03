#devtools::install_github('dviraran/SingleR')
library(SingleR)
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/SingleR_functions.R")
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
#====== 3.1 Create Singler Object  ==========================================
(load(file="data/Lung_28_Nointeg_20200131.Rda"))
object %<>% subset(subset = nFeature_SCT > 200)
object_data <- object[["SCT"]]@data
remove(object); GC()
ref = readRDS(file='data/ref_blue_encode_GSE107011_20200609.rds')

singler <- CreateBigSingleRObject.1(object_data, annot = NULL, project.name="Lung",
                N = 5000, min.genes = 200, technology = "10X",
                species = "Human", citation = "", ref.list = list(ref),
                normalize.gene.length = F, variable.genes = "de", fine.tune = F,
                reduce.file.size = F, do.signatures = F, do.main.types = T,
                temp.dir = getwd(), numCores = SingleR.numCores)
saveRDS(singler,file="output/singlerT_Lung_20200630.rds")
