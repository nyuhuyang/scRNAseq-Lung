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
#====== 2.1 Identify cell types ==========================================
(load(file="data/Lung_28_20200103.Rda"))
DefaultAssay(object)  = "SCT"
(load(file = paste0("output/SCINA_20200110.Rda")))
object@meta.data$SCINA = results$cell_labels
df_cell_types <- readxl::read_excel("doc/Cell type abbreviation.xlsx")
object$SCINA %<>% plyr::mapvalues(
    from = df_cell_types$`Cell types`,
    to = df_cell_types$Abbreviation)
Idents(object) = "SCINA"
object %<>% sortIdent()
object %<>% AddMetaColor(label= "SCINA", colors = c(Singler.colors,Singler.colors))


TSNEPlot.1(object, group.by = "conditions",cols = c("#ffa500","#0000ff","#008000"),label = T,
           label.repel = T, pt.size = 0.5,label.size = 3, repel = T,no.legend = T,
           do.print = T,do.return = F,title = "3 Regions")
UMAPPlot.1(object, group.by = "conditions",cols = c("#ffa500","#0000ff","#008000"),label = T,
           label.repel = T, pt.size = 0.5,label.size = 3, repel = T,no.legend = T,
           do.print = T,do.return = F,title = "3 Regions")

conditions = c("proximal", "distal","terminal","COPD")
Idents(object) = "conditions"
for(i in seq_along(conditions)){
    sub_object <- subset(object, idents = conditions[i])
    Idents(sub_object) = "SCINA"
    
    UMAPPlot.1(sub_object, group.by = "SCINA",cols = ExtractMetaColor(sub_object),label = T,
               label.repel = T, pt.size = 0.5,label.size = 3, repel = T,no.legend = T,
               unique.name = "conditions",
               do.print = T,do.return = F,
               title = paste("Cell types in",conditions[i]))
    TSNEPlot.1(sub_object, group.by = "cell_types",cols = ExtractMetaColor(sub_object),label = T,
               label.repel = T, pt.size = 0.5,label.size = 3, repel = T,no.legend = T,
               unique.name = "conditions",
               do.print = T,do.return = F,
               title = paste("Cell types in",conditions[i]))
    Progress(i,length(conditions))
}

write.csv(table(object$cell_types,object$SCINA),
          paste0(path,"Inherited_vs_SCINA.csv"))