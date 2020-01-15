library(Seurat)
library(dplyr)
library(tidyr)
library(kableExtra)
library(magrittr)
library(readxl)
library(cowplot)
source("../R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)
#====== load data ==========================================
(load(file="data/Lung_28_20200103.Rda"))
(load(file = paste0("output/SCINA_20200110.Rda")))
object@meta.data$SCINA = results$cell_labels
df_cell_types <- readxl::read_excel("doc/Cell type abbreviation.xlsx")
object$SCINA %<>% plyr::mapvalues(
        from = df_cell_types$`Cell types`,
        to = df_cell_types$Abbreviation)
Idents(object) = "SCINA"
object %<>% sortIdent()
object %<>% AddMetaColor(label= "SCINA", colors = c(Singler.colors,Singler.colors))

# ==== Plot ======
for(label in c(TRUE, FALSE)){
        Idents(object) = "SCINA"
        PCAPlot.1(object, group.by = "SCINA",cols = ExtractMetaColor(object),label = label,
                  label.repel = T, pt.size = 0.5,label.size = 3, no.legend = T,
                  unique.name = "conditions",
                  do.print = T,do.return = F,
                  title =  "Cell types in 28 samples")
        conditions = c("proximal", "distal","terminal","COPD")
        Idents(object) = "conditions"
        for(i in seq_along(conditions)){
                sub_object <- subset(object, idents = conditions[i])
                Idents(sub_object) = "SCINA"
                
                PCAPlot.1(sub_object, group.by = "SCINA",cols = ExtractMetaColor(sub_object),
                          label = label,
                          label.repel = T, pt.size = 0.5,label.size = 3,no.legend = T,
                          unique.name = "conditions",
                          do.print = T,do.return = F,
                          title = paste("Cell types in",conditions[i]))
                Progress(i,length(conditions))
        }
}

meta.data = cbind(object@reductions$pca@cell.embeddings[,1:2],
                  object[["cell_types"]]) %>%
        cbind(object[["SCINA"]]) %>%
        cbind(object[["conditions"]])
write.csv(meta.data, paste0(path,"PCA_all_coordinates.csv"))