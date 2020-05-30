########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################
invisible(lapply(c("Seurat","dplyr","kableExtra","cowplot","harmony",
                     "magrittr","data.table","future","ggplot2","tidyr"), function(x) {
                     suppressPackageStartupMessages(library(x,character.only = T))
                             }))
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)
object = readRDS(file = paste0("data/Lung_28_21-ambigous-pca_20200206.rds"))

# draw rectangle
rectangle <- function(x_left, x_right, y_top, y_bottom,...){
        list(geom_segment(aes(x = x_left, xend = x_right, y = y_top, yend = y_top),...),
             geom_segment(aes(x = x_right, xend = x_right, y = y_top, yend = y_bottom),...),
             geom_segment(aes(x = x_left, xend = x_right, y = y_bottom, yend = y_bottom),...),
             geom_segment(aes(x = x_left, xend = x_left, y = y_top, yend = y_bottom),...))
}

jpeg(paste0(path,"FeaturePlot_Cr_HAPLN1.jpeg"), units="in", width=10, height=7,res=600)
FeaturePlot.1(object, features = "HAPLN1")+
        rectangle(-2.5, -1.8, -7.2, -8)
dev.off()

object %<>% AddModuleScore(features = list(c("IGF1", "CCDC80", "SFRP2", "SERPINF1",
                                             "MFAP5", "FBLN2", "FGF7", "PI16", "OGN",
                                             "PCOLCE2")), name = "F")
colnames(object@meta.data)[grep("F1",colnames(object@meta.data))] = "IGF1+CCDC80+SFRP2+SERPINF1+MFAP5+FBLN2+FGF7+PI16+OGN+PCOLCE2"
jpeg(paste0(path,"FeaturePlot_F3_IGF1.jpeg"), units="in", width=10, height=7,res=600)
FeaturePlot.1(object, features = "IGF1+CCDC80+SFRP2+SERPINF1+MFAP5+FBLN2+FGF7+PI16+OGN+PCOLCE2",
              threshold = 0.5,
              title = "IGF1+CCDC80+SFRP2+SERPINF1+MFAP5+FBLN2+FGF7+PI16+OGN+PCOLCE2") +
        xlim(-10,0)+ ylim(-10,0)+
        rectangle(-3.3, -2.7, -6.5, -7.2)+        
        rectangle(-4.6, -4.1, -4.7, -5.5,colour = "blue")+
        rectangle(-3.8, -3.1, -6.5, -6,colour = "blue")
dev.off()

object %<>% AddModuleScore(features = list(c("ACTA2","ELN","COL1A1","BGN")), name = "F5")
colnames(object@meta.data)[grep("F51",colnames(object@meta.data))] = "ACTA2+ELN+COL1A1+BGN"
jpeg(paste0(path,"FeaturePlot_F5.jpeg"), units="in", width=10, height=7,res=600)
FeaturePlot.1(object, features = "ACTA2+ELN+COL1A1+BGN",
              threshold = 0.5,
              title = "ACTA2+ELN+COL1A1+BGN")+ 
        xlim(-10,0)+ ylim(-10,0)+
        rectangle(-4.6, -4.1, -4.7, -5.5,colour = "blue")+
        rectangle(-3.8, -3.1, -6.5, -6,colour = "blue")
dev.off()

object %<>% AddModuleScore(features = list(c("FCER1A", "CCR7", "CD83", "CD1C")), name = "DC")
colnames(object@meta.data)[grep("DC1",colnames(object@meta.data))] = "FCER1A+CCR7+CD83+CD1C"

jpeg(paste0(path,"FeaturePlot_DC~.jpeg"), units="in", width=10, height=7,res=600)
FeaturePlot.1(object, features = "FCER1A+CCR7+CD83+CD1C")+ #xlim(0,4)+ ylim(-6,-3)+
        rectangle(0.7, 1.3, -4.6, -3.9)
dev.off()

# rename ident
annotation <- readxl::read_excel("doc/Harmony annotation Yang.xlsx")
g_anno = annotation[annotation$UMAP %in% "21-ambigous-pca",]
#======== rename ident =================
object = readRDS(file = paste0("data/",rds))
resolutions = g_anno["res"] %>% pull
(res = unique(resolutions))

for(i in seq_along(res)) object %<>% FindClusters(resolution = res[i])

object[["cell.labels"]] = 0
meta.data = object@meta.data
meta.data["barcodes"] = rownames(meta.data)
for(i in seq_along(resolutions)){
        meta.data[meta.data[,paste0("SCT_snn_res.",resolutions[i])] %in% str2int(g_anno$cluster[i]),
                  "cell.labels"] = g_anno$`cell label`[i]
}

m <- cbind(meta.data, object[["umap"]]@cell.embeddings, object[["ACTA2+ELN+COL1A1+BGN"]])
m[m$UMAP_1 > -2.5 & m$UMAP_1 < -1.8 & m$UMAP_2 > -8 & m$UMAP_2 < -7.2,
  "cell.labels"] = "Cr"
m[((m$UMAP_1 > -4.6 & m$UMAP_1 < -4.1 & m$UMAP_2 < -4.7 & m$UMAP_2 > -5.5) |
           (m$UMAP_1 > -3.8 & m$UMAP_1 < -3.1 & m$UMAP_2 < -6 & m$UMAP_2 > -6.5)) &
          m$`ACTA2+ELN+COL1A1+BGN` > 0.5,
  "cell.labels"] = "F5"  
m[m$UMAP_1 > 0.7 & m$UMAP_1 < 1.3 & m$UMAP_2 > -4.6 & m$UMAP_2 < -3.9,
  "cell.labels"] = "DC"
object@meta.data = m

Idents(object) = "cell.labels"
UMAPPlot.1(object, do.print = T, do.return = F, cols =  Singler.colors,
           label = T, label.repel = T)
