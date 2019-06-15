library(Seurat)
library(dplyr)
library(tidyr)
library(kableExtra)
library(magrittr)
library(readxl)
library(cowplot)
source("../R/Seurat3_functions.R")
path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)
#====== 2.1 pathway analysis ==========================================
(load(file="data/Lung_harmony_12_20190614.Rda"))

df_markers <- readxl::read_excel("doc/Renat.markers.xlsx",sheet = "20190613")

#df_markers <- readxl::read_excel("../seurat_resources/bio-rad-markers.xlsx",sheet = "Human.sub")
colnames(df_markers) = gsub(" ","_",colnames(df_markers))
colnames(df_markers) = gsub(":|\\/","_",colnames(df_markers))
colnames(df_markers) = gsub("\\+","",colnames(df_markers))
#markers = df_markers[,-grep("Alias",colnames(df_markers))]
marker.list <- df2list(df_markers)

marker.list %<>% lapply(function(x) x) %>% 
     lapply(function(x) FilterGenes(object,x)) %>% 
     lapply(function(x) x[!is.na(x)])
#marker.list %>% list2df %>% t %>% kable() %>% kable_styling()
Idents(object) <- "RNA_snn_res.0.6"

for(i in 1:length(marker.list)){
    p <- lapply(marker.list[[i]], function(marker) {
        FeaturePlot(object = object, features = marker,pt.size = 0.5, label=T)+
            NoLegend()+
            ggtitle(paste0(marker,Alias(df = df_markers,gene = marker)))+
            theme(plot.title = element_text(hjust = 0.5,size = 15,face = "plain"))
    })
    jpeg(paste0(path,names(marker.list)[i],".jpeg"),units="in", width=10, height=7,res=600)
    print(do.call(plot_grid, p)+ ggtitle(paste(names(marker.list)[i],"markers"))+
              theme(plot.title = element_text(hjust = 0.5,size = 20)))
    dev.off()
    print(paste0(i,":",length(marker.list)))
}

#======== rename ident =================
object$cell.type <- plyr::mapvalues(object$RNA_snn_res.0.6,
                                            from = 0:19,
                                            to = c("Alveolar macrophages",
                                                   "Alveolar type I&II cells \nDistal secretory cells",
                                                   "Endothelial cells",
                                                   "Macrophages",
                                                   "T cells",
                                                   "Stromal fibroblasts",
                                                   "Monocytes",
                                                   "Ciliated cells",
                                                   "Smooth muscle cells& \nStromal fibroblasts",
                                                   "Secretory cells",
                                                   "Endothelial cells",
                                                   "Alveolar type I&II cells \nDistal secretory cells",
                                                   "B cells",
                                                   "Alveolar macrophages& \nRed blood cells",
                                                   "Basal cells",
                                                   "HSC/progenitor cells",
                                                   "Lymphatic endothelial cells",
                                                   "Endothelial cells",
                                                   "Chondrocytes & Fibroblast",
                                                   "Alveolar type I&II cells \nDistal secretory cells"))
object$cell.type <- as.character(object$cell.type)

# rename NK cells
(load(file="output/singlerF_Lung_12_20190614.Rda"))
singlerDF = data.frame("singler1sub" = singler$singler[[1]]$SingleR.single$labels,
                       row.names = rownames(singler$singler[[1]]$SingleR.single$labels))
head(singlerDF)
NK_cells <-  rownames(singlerDF)[singlerDF$singler1sub %in% "NK_cells"]
# PECAM1_neg
SELE <- WhichCells(object, expression = SELE >0 )


object@meta.data[NK_cells,"cell.type"] = "NK cells"
object@meta.data[SELE,"cell.type"] = "Endothelial PECAM1-neg"

Idents(object) <- "cell.type"
TSNEPlot.1(object, cols = ExtractMetaColor(object),label = T,pt.size = 1,
           label.size = 5, repel = T,no.legend = F,do.print = F,
           title = "Cell types")
##############################
# process color scheme
##############################

table(Idents(object)) %>% kable %>% kable_styling()
singler_colors <- readxl::read_excel("doc/singler.colors.xlsx")
singler_colors1 = as.vector(singler_colors$singler.color1[!is.na(singler_colors$singler.color1)])
singler_colors1[duplicated(singler_colors1)]
length(singler_colors1)
length(unique(object$cell.type))
object <- AddMetaColor(object = object, label= "cell.type", colors = singler_colors1)
Idents(object) <- "cell.type"
object <- sortIdent(object)

TSNEPlot.1(object, cols = ExtractMetaColor(object),label = T,pt.size = 1,
           label.size = 5, repel = T,no.legend = T,do.print = T,
           title = "Cell types")
save(object,file="data/Lung_harmony_12_20190614.Rda")
