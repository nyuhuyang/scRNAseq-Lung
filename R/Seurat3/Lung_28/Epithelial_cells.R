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
#===================
#Idents(object) = "SCT_snn_res.0.8"
#jpeg(paste0(path,"UMAPPlot_object_SCT_snn_res.0.8.jpeg"), 
#    units="in", width=10, height=7,res=600)
#UMAPPlot.1(object,group.by = "SCT_snn_res.0.8", label = T,do.return = T)+
#        geom_vline(xintercept=-5)+
#        geom_hline(yintercept=-1.5)
#dev.off()
Epi <- subset(object,subset = UMAP_1 > -5 & UMAP_2 < -1.5)
#UMAPPlot.1(Epi,group.by = "SCT_snn_res.0.8", label = T,label.repel = T,
#           do.print = T)
df_cell_types$Epi %<>% as.logical()
Epi_cell_types <- df_cell_types$Abbreviation[df_cell_types$Epi]
Idents(Epi) = "SCINA"
Epi <- subset(Epi, idents = Epi_cell_types)

for(label in c(TRUE, FALSE)){
        Idents(Epi) = "SCINA"
        PCAPlot.1(Epi, group.by = "SCINA",cols = ExtractMetaColor(Epi),label = label,
                  label.repel = T, pt.size = 0.5,label.size = 3, no.legend = T,
                  unique.name = "conditions",
                  do.print = T,do.return = F,
                  title =  "Epithelial cell types in 28 samples")
        conditions = c("proximal", "distal","terminal","COPD")
        Idents(Epi) = "conditions"
        for(i in seq_along(conditions)){
                sub_Epi <- subset(Epi, idents = conditions[i])
                Idents(sub_Epi) = "SCINA"
                
                PCAPlot.1(sub_Epi, group.by = "SCINA",cols = ExtractMetaColor(sub_Epi),
                          label = label,
                          label.repel = T, pt.size = 0.5,label.size = 3,no.legend = T,
                          unique.name = "conditions",
                          do.print = T,do.return = F,
                          title = paste("Epithelial cell types in",conditions[i]))
                Progress(i,length(conditions))
        }
}

meta.data = cbind(Epi@reductions$pca@cell.embeddings[,1:2],
                  Epi[["cell_types"]]) %>%
        cbind(Epi[["SCINA"]]) %>%
        cbind(Epi[["conditions"]])
write.csv(meta.data, paste0(path,"PCA_Epi_coordinates.csv"))

c("AGER","CLDN18","AQP4","EMP2","SCEL")