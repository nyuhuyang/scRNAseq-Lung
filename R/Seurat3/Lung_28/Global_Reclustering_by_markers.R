########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################
#devtools::install_github("immunogenomics/harmony", ref= "ee0877a",force = T)
invisible(lapply(c("Seurat","dplyr","cowplot","magrittr","MAST",
                   "future","ggplot2","tidyr","harmony"), function(x) {
                           suppressPackageStartupMessages(library(x,character.only = T))
                   }))
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)

names = c("Singlets-harmony","Singlets-harmony-unlabeled","All-harmony", "All-harmony-unlabeled")
for(i in seq_along(data)){
        save.path <- paste0(path,data[i],"/ReductionsPlots/")
        if(!dir.exists(save.path))dir.create(save.path, recursive = T)
        
        object = readRDS(file = paste0("data/Lung_28_",data[i],"_20200522.rds"))

        Idents(object) = "orig.ident"
        lapply(c(TRUE, FALSE), function(lab) 
                UMAPPlot.1(object, group.by="orig.ident",pt.size = 0.5,label = lab,
                           cols = Singler.colors,
                           label.repel = T,alpha = 0.9,
                           no.legend = T,label.size = 4, repel = T, title = "Harmony Integration",
                           do.print = T, do.return = F,save.path = save.path))
        
        Idents(object) = "cell_types"
        lapply(c(TRUE, FALSE), function(lab)
                UMAPPlot.1(object, group.by="cell_types",pt.size = 0.5,label = lab,
                           label.repel = T,alpha = 0.9,cols = Singler.colors,
                           no.legend = T,label.size = 4, repel = T, title = "First annotation",
                           do.print = T, do.return = F,save.path = save.path))
        Idents(object) = "cell.labels"
        lapply(c(TRUE, FALSE), function(lab)
                UMAPPlot.1(object, group.by="cell.labels",pt.size = 0.5,label = lab,
                           label.repel = T,alpha = 0.9,cols = Singler.colors,
                           no.legend = T,label.size = 4, repel = T, title = "2nd annotation",
                           do.print = T, do.return = F,save.path = save.path))
        Idents(object) = "annotations"
        lapply(c(TRUE, FALSE), function(lab)
                UMAPPlot.1(object, group.by="annotations",pt.size = 0.5,label = lab,
                           cols = ExtractMetaColor(object),
                           label.repel = T, alpha= 0.9,
                           no.legend = T,label.size = 4, repel = T, title = "Last annotation",
                           do.print = T, do.return = F,save.path = save.path))
        FeaturePlot.1(object, features = c("MS4A1","CYTL1","TPSAB1",
                                           "DKK2","CHGA","L1CAM",
                                           "SCGB3A2","MZB1","FOXI1"),
                      ncol = 3,do.print = T, save.path = save.path)
}