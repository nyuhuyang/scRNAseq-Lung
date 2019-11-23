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
(load(file="data/Lung_24_20190918.Rda"))
DefaultAssay(object) <- 'RNA'
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
marker.list <- marker.list[!is.na(marker.list)]
marker.list <- marker.list[sapply(marker.list,length)>0]
marker.list %>% list2df %>% t %>% kable() %>% kable_styling()
Idents(object) <- "RNA_snn_res.0.8"

for(i in 1:length(marker.list)){
    p <- lapply(marker.list[[i]], function(marker) {
        FeaturePlot(object = object, features = marker,pt.size = 0.5, label=T,
                    reduction = "tsne")+
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
Idents(object) = "RNA_snn_res.0.8"
object %<>% RenameIdents("0" = "T cells:Central naive(Tcn)",
                         "1" = "Ciliated cells:1",
                         "2" = "T cells:Resident memory(Trm)",
                         "3" = "Endothelial:HEV",
                         "4" = "Monocytes",
                         "5" = "Macrophages",
                         "6" = "Neutrophils",
                         "7" = "B cells",
                         "8" = "NK/T cytotoxic",
                         "9" = "Smooth muscle cells:1",
                         "10" = "Endothelial cells:Capillary:1",
                         "11" = "Secretory cells:1",
                         "12" = "Fibroblasts:1",
                         "13" = "T cells:7SK.2+",
                         "14" = "Alveolar type 2 cells:A",
                         "15" = "Alveolar type 2 cells:B",
                         "16" = "Ciliated cells:2",
                         "17" = "Intermediate cells:1",
                         "18" = "Endothelial cells:Capillary:2",
                         "19" = "Fibroblasts:2",
                         "20" = "Basal cells",
                         "21" = "Submucosal gland serous cells",
                         "22" = "Mucus-producing cells",
                         "23" = "Macrophages:2",
                         "24" = "Mast cells",
                         "25" = "Hybrid cells",
                         "26" = "Ciliated cells:3",
                         "27" = "Smooth muscle cells:2",
                         "28" = "Secretory cells:2",
                         "29" = "Endothelial cells:Lymphatic",
                         "30" = "Smooth muscle cells:3",
                         "31" = "Chondrocytes",
                         "32" = "Pre-ciliated cells",
                         "33" = "Endothelial cells:Smooth muscle-like",
                         "34" = "Plasma cells",
                         "35" = "Proliferating cells",
                         "36" = "T cells:Interferon-responsive",
                         "37" = "Endothelial cells:Arterial",
                         "38" = "Neuroendocrine cells",
                         "39" = "Alveolar type 1 cells",
                         "40" = "Neurons",
                         "41" = "Dendritic Cells:Plasmacytoid")
object@meta.data$cell.types = Idents(object)
object@meta.data$cell.types %<>% as.character()
meta.data = cbind(object@meta.data, object@reductions$tsne@cell.embeddings)
# Cluster 5 – Macrophages
DC <- meta.data$RNA_snn_res.0.8 %in% 5 & meta.data$tSNE_2 < -21
object@meta.data[DC,"cell.types"] = "Dendritic cells:2"
# Cluster 12
object@meta.data[meta.data$RNA_snn_res.0.8 %in% 12 & meta.data$tSNE_2 < 18,
                 "cell.types"] = "Fibroblasts:6"
meta.data <- cbind(meta.data,FetchData(object, vars = "MFAP5"))
object@meta.data[meta.data$RNA_snn_res.0.8 %in% 12 & meta.data$MFAP5 > 1,
                 "cell.types"] = "Fibroblasts:7"
# Cluster 21 – Submucosal gland serous cells
meta.data <- cbind(meta.data,FetchData(object, vars = "MUC5B"))
object@meta.data[meta.data$RNA_snn_res.0.8 %in% 21 & meta.data$tSNE_1 < 0,
                 "cell.types"] = "T cells:7SK.2+"
object@meta.data[meta.data$RNA_snn_res.0.8 %in% 21 & meta.data$tSNE_1 > 7 ,
                 "cell.types"] = "Submucosal gland mucous cells"
# Rename cluster 22:
object@meta.data[meta.data$RNA_snn_res.0.8 %in% 22 & meta.data$tSNE_1 < 19,
                "cell.types"] = "Intermediate cells:2"
object@meta.data[meta.data$RNA_snn_res.0.8 %in% 22 & meta.data$tSNE_1 >= 19,
                 "cell.types"] = "Mucus-producing cells"
# Rename cluster 31:

object@meta.data[meta.data$RNA_snn_res.0.8 %in% 31 & meta.data$tSNE_2 > 10,
                 "cell.types"] = "Unknwon"
# Cluster 35 – Proliferating cells
meta.data <- cbind(meta.data,FetchData(object, vars = c("MKI67","CDH5","CLDN5",
                                                        "KRT5","DCN","SFTPC",
                                                        "PTPRC","CD3G","CD3E",
                                                        "CD163","CD68")))
object@meta.data[meta.data$RNA_snn_res.0.8 %in% 35 & meta.data$MKI67 >0 & 
                meta.data$CDH5 >0 & meta.data$CLDN5 >0, 
                "cell.types"] = "Endothelial cells:Proliferating"
object@meta.data[meta.data$RNA_snn_res.0.8 %in% 35 & meta.data$MKI67 >0 & 
                     meta.data$KRT5 >0, "cell.types"] = "Basal cells:Proliferating"
object@meta.data[meta.data$RNA_snn_res.0.8 %in% 35 & meta.data$MKI67 >0 & 
                 meta.data$DCN >0,"cell.types"] = "Fibroblasts:Proliferating"
object@meta.data[meta.data$RNA_snn_res.0.8 %in% 35 & meta.data$MKI67 >0 & 
                 meta.data$SFTPC >0,"cell.types"] = "Alveolar type 2 cells:Proliferating"
object@meta.data[meta.data$RNA_snn_res.0.8 %in% 35 & meta.data$MKI67 >0 & 
                     meta.data$PTPRC >0 & meta.data$CD3G >0, 
                 "cell.types"] = "T cells:Proliferating"
object@meta.data[meta.data$RNA_snn_res.0.8 %in% 35 & meta.data$MKI67 >0 & 
                     meta.data$CD163 >0 & meta.data$CD68 >0, 
                 "cell.types"] = "Macrophages:Proliferating"
#object %<>% FindClusters(resolution = 1.5)
object@meta.data[meta.data$RNA_snn_res.1.5 %in% 49,"cell.types"] = "Squamous cells"
object@meta.data[meta.data$RNA_snn_res.1.5 %in% 29 & meta.data$tSNE_1 >20,"cell.types"] = "Mucus-producing cells"
object@meta.data[meta.data$RNA_snn_res.1.5 %in% 28,"cell.types"] = "Fibroblasts:3"
object@meta.data[meta.data$RNA_snn_res.1.5 %in% 28 & 
                     meta.data$tSNE_1 < 2 & meta.data$tSNE_2 < 28,"cell.types"] = "Fibroblasts:4"
object@meta.data[meta.data$RNA_snn_res.1.5 %in% 30 & meta.data$tSNE_1 >20,"cell.types"] = "Hybrid cells:2"

#object %<>% FindClusters(resolution = 1.6)
meta.data <- cbind(meta.data,FetchData(object, vars = "SEPP1"))
object@meta.data[meta.data$RNA_snn_res.1.6 %in% 11 & meta.data$SEPP1 > 1 ,"cell.types"] = "Dendritic Cells:1"

object@meta.data[meta.data$RNA_snn_res.1.6 %in% 28 & meta.data$tSNE_1 >5 ,"cell.types"] = "Fibroblasts:5"
object@meta.data[meta.data$RNA_snn_res.1.6 %in% 33,"cell.types"] = "Secretory cells:Distal"
object@meta.data[meta.data$RNA_snn_res.1.6 %in% 34 & 
                     meta.data$tSNE_1 >5 & meta.data$tSNE_1 < 12,"cell.types"] = "Secretory cells:Distal"
object@meta.data[meta.data$RNA_snn_res.1.6 %in% 36,"cell.types"] = "Endothelial cells:Capillary:3"
object@meta.data[meta.data$RNA_snn_res.1.6 %in% 37 & meta.data$tSNE_2 < -15,"cell.types"] = "Endothelial cells:Proliferating"

object@meta.data[meta.data$RNA_snn_res.1.6 %in% 41,"cell.types"] = "Pericytes"
object@meta.data[meta.data$RNA_snn_res.1.6 %in% 50 & meta.data$tSNE_2 > 10 ,"cell.types"] = "Ionocytes"
object@meta.data[meta.data$RNA_snn_res.1.6 %in% 50 & meta.data$tSNE_1 > 5 ,"cell.types"] = "Unknown"
object@meta.data[meta.data$RNA_snn_res.1.6 %in% 54,"cell.types"] = "Myoepithelial cells"
meta.data <- cbind(meta.data,FetchData(object, vars = "FOXI1"))
object@meta.data[meta.data$RNA_snn_res.0.8 %in% 38 & meta.data$FOXI1 >0,
                 "cell.types"] = "Ionocytes"
inherited <- read.csv(paste0("output/Lung_24_20191105_cell_types.csv"))
inherited$barcode %<>% as.character()
object@meta.data[inherited[inherited$cell.types %in% "Mucous gland cells","barcode"],
                 "cell.types"] = "Submucosal gland mucous cells"

#=======
object$conditions %<>% factor(levels = c("proximal", "distal", "terminal"))
df <- table(object$cell.types,object$conditions) %>% as.data.frame() %>% 
    spread(Var2,Freq)
colnames(df)[1] = "cell.types"
write.csv(df, file = paste0(path,"Cell_types_by_regions.csv"))



#========

Idents(object) = "cell.types"
object %<>% sortIdent()
object <- AddMetaColor(object = object, label= "cell.types", colors = Singler.colors)
for( label in c(T,F)){
    TSNEPlot.1(object, group.by = "cell.types",cols = ExtractMetaColor(object),label = label,
               label.repel = T, pt.size = 0.5,label.size = 3, repel = T,no.legend = T,
               do.print = T,do.return = F,title = "Cell types in all 24 samples")
    UMAPPlot.1(object, group.by = "cell.types",cols = ExtractMetaColor(object),label = label,
               label.repel = T, pt.size = 0.5,label.size = 3, repel = T,no.legend = T,
               do.print = T,do.return = F,title = "Cell types in all 24 samples")
}

TSNEPlot.1(object, group.by = "conditions",cols = c("#ffa500","#0000ff","#008000"),label = T,
           label.repel = T, pt.size = 0.5,label.size = 3, repel = T,no.legend = T,
           do.print = T,do.return = F,title = "3 Regions")
UMAPPlot.1(object, group.by = "conditions",cols = c("#ffa500","#0000ff","#008000"),label = T,
           label.repel = T, pt.size = 0.5,label.size = 3, repel = T,no.legend = T,
           do.print = T,do.return = F,title = "3 Regions")


object@meta.data$conditions %<>% as.character()
Idents(object) = "conditions"
for (con in c("distal","proximal","terminal")){
    subset_object <- subset(object,idents = con)
    Idents(subset_object) = "cell.types"
    subset_object %<>% sortIdent()

    for( label in c(T, F)){
        TSNEPlot.1(subset_object, group.by = "cell.types",
                   cols = ExtractMetaColor(subset_object),label = label,unique.name = "conditions",
                   label.repel = T, pt.size = 0.5,label.size = 3, repel = T,no.legend = T,
                   do.print = T,do.return = F,title = paste("Cell types in",con))
        UMAPPlot.1(subset_object, group.by = "cell.types",
                   cols = ExtractMetaColor(subset_object),label = label,unique.name = "conditions",
                   label.repel = T, pt.size = 0.5,label.size = 3, repel = T,no.legend = T,
                   do.print = T,do.return = F,title = paste("Cell types in",con))
    }
}


save(object,file="data/Lung_24_20190918.Rda")