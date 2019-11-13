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
                         "1" = "Ciliated cells-1",
                         "2" = "T cells:Resident memory(Trm)",
                         "3" = "Endothelial:HEV",
                         "4" = "Monocytes",
                         "5" = "Macrophages",
                         "6" = "Neutrophils",
                         "7" = "B cells",
                         "8" = "NK/T cytotoxic",
                         "9" = "Smooth muscle cells-1",
                         "10" = "Endothelial cells:Capillary-1",
                         "11" = "Secretory cells-1",
                         "12" = "Fibroblasts-1",
                         "13" = "T cells:7SK.2+",
                         "14" = "Alveolar type 2 cells:A",
                         "15" = "Alveolar type 2 cells:B",
                         "16" = "Ciliated cells-2",
                         "17" = "Intermediate cells-1",
                         "18" = "Endothelial cells:Capillary-2",
                         "19" = "Fibroblasts-2",
                         "20" = "Basal cells",
                         "21" = "Submucosal gland serous cells",
                         "22" = "Mucus-producing cells",
                         "23" = "Macrophages:Alveolar",
                         "24" = "Mast cells",
                         "25" = "Hybrid cells",
                         "26" = "Ciliated cells-3",
                         "27" = "Smooth muscle cells-2",
                         "28" = "Secretory cells-2",
                         "29" = "Endothelial cells:Lymphatic",
                         "30" = "Smooth muscle cells-3",
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
meta.data = cbind(object@meta.data,object@reductions$tsne@cell.embeddings)
# Cluster 5 – Macrophages
DC <- meta.data$RNA_snn_res.0.8 %in% 5 & meta.data$tSNE_2 < -21
object@meta.data[DC,"cell.types"] = "Dendritic cells"
# Cluster 21 – Submucosal gland serous cells
meta.data <- cbind(meta.data,FetchData(object, vars = "MUC5B"))
object@meta.data[meta.data$MUC5B >2.5 & meta.data$RNA_snn_res.0.8 %in% 21,
                 "cell.types"] = "Submucosal gland mucous cells"
# Rename cluster 22:
object@meta.data[meta.data$RNA_snn_res.0.8 %in% 22 & meta.data$tSNE_1 < 19,
                "cell.types"] = "Mucus-producing cells-1"
object@meta.data[meta.data$RNA_snn_res.0.8 %in% 22 & meta.data$tSNE_1 >= 19,
                 "cell.types"] = "Mucus-producing cells-2"
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

object %<>% FindClusters(resolution = 1.6)
object@meta.data[meta.data$RNA_snn_res.1.6 %in% 28,"cell.types"] = "Fibroblasts-3"
object@meta.data[meta.data$RNA_snn_res.1.6 %in% 33,"cell.types"] = "Secretory cells:Distal"
object@meta.data[meta.data$RNA_snn_res.1.6 %in% 41,"cell.types"] = "Pericytes"
object@meta.data[meta.data$RNA_snn_res.1.6 %in% 54,"cell.types"] = "Myoepithelial cells (KRT5+ ACTA2+)"
meta.data <- cbind(meta.data,FetchData(object, vars = "FOXI1"))
object@meta.data[meta.data$RNA_snn_res.0.8 %in% 38 & meta.data$FOXI1 >0,
                 "cell.types"] = "Ionocytes"

Idents(object) = "cell.types"
 
object %<>% sortIdent()
object <- AddMetaColor(object = object, label= "cell.types", colors = Singler.colors)
Idents(object) = "conditions"
distal <- subset(object,idents = "distal")
terminal <- subset(object,idents = "terminal")
proximal <- subset(object,idents = "proximal")

Idents(proximal) = "cell.types"
proximal %<>% sortIdent()
TSNEPlot.1(proximal, group.by = "cell.types",cols = ExtractMetaColor(proximal),label = T,
           label.repel = T, pt.size = 0.5,label.size = 3, repel = T,no.legend = T,
           do.print = T,do.return = F,title = "Cell types in proximal")

object$conditions %<>% factor(levels = c("proximal", "distal", "terminal"))
df <- table(object$cell.types,object$conditions) %>% as.data.frame() %>% 
    spread(Var2,Freq)
colnames(df)[1] = "cell.types"
write.csv(df, file = paste0(path,"Cell_types_by_regions.csv"))
save(object,file="data/Lung_24_20190918.Rda")