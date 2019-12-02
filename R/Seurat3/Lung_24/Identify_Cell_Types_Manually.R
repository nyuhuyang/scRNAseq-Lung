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
(load(file="data/Lung_24_20191128.Rda"))
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
object %<>% RenameIdents(#"0" = #"T cells:Central naive(Tcn)",
                         #"1" = #"Ciliated cells:1",
                         "2" = "Ciliated cells:1",
                         "3" = "Endothelial cells:Capillary:1",
                         "4" = "Endothelial cells:HEV",
                         "5" = "Monocytes",
                         #"6" = #"Neutrophils",
                         #"7" = #"B cells",
                         "8" = "Neutrophils",
                         "9" = "B cells",
                         #"10" = #"Endothelial cells:Capillary:1",
                         "11" = "Ciliated cells:2",#"Secretory cells:1",
                         #"12" = #
                         "13" = "Alveolar type 2 cells:A",
                         "14" = "Alveolar type 2 cells:B",
                         "15" = "Dendritic cells:MDC1",
                         #"16" = #"Ciliated cells:2",
                         #"17" = #,
                         #"18" = 
                         "19" = "Macrophages:1",
                         #"20" = #"Basal cells",
                         "21" = "Macrophages:2",
                         "22" = "Submucosal gland:Serous cells",
                         "23" = "Mast cells",#"Macrophages:2",
                         "24" = "Hybrid cells",
                         "25" = "Endothelial cells:Smooth muscle-like",
                         "26" = "Ciliated cells:3",
                         #"27" = #
                         "28" = "Endothelial cells:Lymphatic",#"Secretory cells:2",
                         "29" = "Smooth muscle cells:3",
                         "30" = "Chondrocytes",
                         "31" = "Pre-ciliated cells",
                         "32" = "Plasma cells",
                         "33" = "T cells:IFN-stimulated",
                         #"34" = 
                         "35" = "Neuroendocrine cells",
                         "36" = "Alveolar type 1 cells",
                         "37" = "Neurons",
                         "38" = "Myoepithelial cells")
object@meta.data$cell.types = Idents(object)
object@meta.data$cell.types %<>% as.character()
meta.data = cbind(object@meta.data, object@reductions$tsne@cell.embeddings)

Idents(object) = "RNA_snn_res.0.8"
meta.data <- cbind(meta.data,FetchData(object, vars = c("IL7R","FCN3","CD36","CXCL12","FBLN5","CA4")))
object@meta.data[meta.data$RNA_snn_res.0.8 %in% 4 & meta.data$FCN3 > 1.5 & meta.data$tSNE_1 < -10,
                 "cell.types"] = "Endothelial cells:Capillary:2"
object@meta.data[meta.data$RNA_snn_res.0.8 %in% 4 & meta.data$FBLN5>0.1,
                 "cell.types"] = "Endothelial cells:Artery:2"
object@meta.data[meta.data$RNA_snn_res.0.8 %in% 4 & meta.data$CA4>0.1 & meta.data$FBLN5 <= 0.1,
                 "cell.types"] = "Endothelial cells:Capillary:3"

meta.data <- cbind(meta.data,FetchData(object, vars = c("SEPP1","IL3RA")))
object@meta.data[meta.data$RNA_snn_res.0.8 %in% 15 & meta.data$SEPP1>0,
                 "cell.types"] = "Dendritic cells:MDC2"
object@meta.data[meta.data$RNA_snn_res.0.8 %in% 15 & meta.data$tSNE_1 > -14.5 &
                     meta.data$tSNE_1 < -10 & meta.data$tSNE_2  >3,
                 "cell.types"] = "Dendritic cells:Plasmacytoid"
object@meta.data[meta.data$RNA_snn_res.0.8 %in% 22 & meta.data$tSNE_2 > -20,
                 "cell.types"] = "Submucosal gland:Mucous cells"

#object %<>% FindClusters(resolution = 1.5)
meta.data = cbind(object@meta.data, object@reductions$tsne@cell.embeddings)

object@meta.data[meta.data$RNA_snn_res.1.5 %in% 0,"cell.types"] = "T cells:Central memory,naive(Tcn)"
object@meta.data[meta.data$RNA_snn_res.1.5 %in% 1,"cell.types"] = "T cells:Resident memory(Trm)"
object@meta.data[meta.data$RNA_snn_res.1.5 %in% 7,"cell.types"] = "T cells:7SK.2+"
object@meta.data[meta.data$RNA_snn_res.1.5 %in% 8,"cell.types"] = "Fibroblasts:1"
object@meta.data[meta.data$RNA_snn_res.1.5 %in% 10,"cell.types"] = "NK/T cytotoxic cells"
object@meta.data[meta.data$RNA_snn_res.1.5 %in% 15,"cell.types"] = "Smooth muscle cells:1"

meta.data <- cbind(meta.data,FetchData(object, vars = c("WNT2","ACTA2")))
object@meta.data[meta.data$RNA_snn_res.1.5 %in% 28,"cell.types"] = "Fibroblasts:3"
object@meta.data[meta.data$RNA_snn_res.1.5 %in% 28 & meta.data$WNT2 > 0,"cell.types"] = "Fibroblasts:4"
object@meta.data[meta.data$RNA_snn_res.1.5 %in% 28 & meta.data$tSNE_1 > 26.7,"cell.types"] = "Fibroblasts:5"

object@meta.data[meta.data$RNA_snn_res.1.5 %in% 30 & meta.data$tSNE_1 < 8.5,"cell.types"] = "Pericytes"
object@meta.data[meta.data$RNA_snn_res.1.5 %in% 30 & meta.data$tSNE_1 > 8.5,"cell.types"] = "Smooth muscle cells:2"

object@meta.data[meta.data$RNA_snn_res.1.5 %in% 31,"cell.types"] = "Fibroblasts:2"
object@meta.data[meta.data$RNA_snn_res.1.5 %in% 34,"cell.types"] = "Smooth muscle cells:2"
object@meta.data[meta.data$RNA_snn_res.1.5 %in% 43,"cell.types"] = "Endothelial cells:Capillary:5"
object@meta.data[meta.data$RNA_snn_res.1.5 %in% 45, "cell.types"] = "Squamous"
object@meta.data[meta.data$RNA_snn_res.1.5 %in% 46 & meta.data$tSNE_1 > 15,"cell.types"] = "Ionocytes"
object@meta.data[meta.data$RNA_snn_res.1.5 %in% 46 & meta.data$tSNE_2 > -10,"cell.types"] = "Unknown"
object@meta.data[meta.data$RNA_snn_res.1.5 %in% 46 & meta.data$tSNE_1 < 0,"cell.types"] = "Unknown"

 
#object %<>% FindClusters(resolution = 1.6)
meta.data = cbind(object@meta.data, object@reductions$tsne@cell.embeddings)
object@meta.data[meta.data$RNA_snn_res.1.6 %in% 15,"cell.types"] = "Basal cells"
object@meta.data[meta.data$RNA_snn_res.1.6 %in% 17,"cell.types"] = "Secretory cells:1"
object@meta.data[meta.data$RNA_snn_res.1.6 %in% 24,"cell.types"] = "Intermediate cells:1"

meta.data <- cbind(meta.data,FetchData(object, vars = c("MUC5AC")))
object@meta.data[meta.data$RNA_snn_res.1.6 %in% 21 & meta.data$MUC5AC > 0 ,"cell.types"] = "Mucus-producing cells"
object@meta.data[meta.data$RNA_snn_res.1.6 %in% 21 & meta.data$MUC5AC == 0 ,"cell.types"] = "Intermediate cells-2"
meta.data <- cbind(meta.data,FetchData(object, vars = c("SFTPB")))
object@meta.data[meta.data$RNA_snn_res.1.6 %in% 27 & meta.data$SFTPB > 0 ,"cell.types"] = "Secretory cells:Distal"
object@meta.data[meta.data$RNA_snn_res.1.6 %in% 27 & meta.data$SFTPB == 0 ,"cell.types"] = "Secretory cells:3"
object@meta.data[meta.data$RNA_snn_res.1.6 %in% 28,"cell.types"] = "Endothelial cells:Capillary:4"
object@meta.data[meta.data$RNA_snn_res.1.6 %in% 30 & meta.data$tSNE_2 > 15 & meta.data$tSNE_2 < 20 ,
                 "cell.types"] = "Endothelial cells:Artery:1"
object@meta.data[meta.data$RNA_snn_res.1.6 %in% 30 & meta.data$tSNE_2 > 20 ,
                 "cell.types"] = "Endothelial cells:Capillary:5"
object@meta.data[meta.data$RNA_snn_res.1.6 %in% 37,"cell.types"] = "Secretory cells:2"
object@meta.data[meta.data$RNA_snn_res.1.6 %in% 44,"cell.types"] = "Alveolar type 2 cells:C"
object@meta.data[meta.data$RNA_snn_res.1.6 %in% 44 & meta.data$RNA_snn_res.0.8 %in% 4,"cell.types"] = "Alveolar type 2 cells:c"
object@meta.data[meta.data$RNA_snn_res.1.6 %in% 45 & meta.data$tSNE_1 > -8 &
                 meta.data$tSNE_1 < 0 & meta.data$tSNE_2  >1 & meta.data$tSNE_2 < 3,"cell.types"] = "Endothelial cells:Proliferating"
object@meta.data[meta.data$RNA_snn_res.1.6 %in% 45 & meta.data$tSNE_1 > -8 &
                     meta.data$tSNE_1 < 0 & meta.data$tSNE_2  >1 & meta.data$tSNE_2 < 3,"cell.types"] = "Endothelial cells:Proliferating"
object@meta.data[meta.data$RNA_snn_res.1.6 %in% 45 & meta.data$tSNE_1 > -15 &
                     meta.data$tSNE_1 < -10 & meta.data$tSNE_2  > -3 & meta.data$tSNE_2 < 1,"cell.types"] = "T cells:Proliferating"
object@meta.data[meta.data$RNA_snn_res.1.6 %in% 45 & meta.data$tSNE_1 > -20 &
                     meta.data$tSNE_1 < -12 & meta.data$tSNE_2  > 1 & meta.data$tSNE_2 < 10,"cell.types"] = "Macrophages:Proliferating"
object@meta.data[meta.data$RNA_snn_res.1.6 %in% 47, "cell.types"] = "Basal cells:Proliferating"

#== clear up === 
object@meta.data[object$cell.types %in% 7, "cell.types"] = "Stromal cells:Proliferating"
object@meta.data[object$cell.types %in% 12, "cell.types"] = "Fibroblasts:6"
object@meta.data[object$cell.types %in% 10, "cell.types"] = "Secretory cells:Distal:2"

object@meta.data[object$cell.types %in% 17:18, "cell.types"] = "Intermediate cells:1"
object@meta.data[object$cell.types %in% c(0,1,6,20,34), "cell.types"] = "Unknown"

#=======
object$conditions %<>% factor(levels = c("proximal", "distal", "terminal"))
df <- table(object$cell.types,object$conditions) %>% as.data.frame() %>% 
    spread(Var2,Freq)
colnames(df)[1] = "cell.types"
write.csv(df, file = paste0(path,"Cell_types_by_regions.csv"))

df <- table(object$cell.types,object$orig.ident) %>% as.data.frame() %>% 
    spread(Var2,Freq)
colnames(df)[1] = "cell.types"
write.csv(df, file = paste0(path,"Cell_types_by_samples.csv"))

#========

Idents(object) = "cell.types"
object %<>% sortIdent()
object <- AddMetaColor(object = object, label= "cell.types", colors = c(Singler.colors,Singler.colors))
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

save(object,file="data/Lung_24_20191128.Rda")
