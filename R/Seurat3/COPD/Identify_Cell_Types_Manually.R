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
(load(file = "data/Lung_16_distal_20191017.Rda"))
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
Idents(object) <- "integrated_snn_res.0.8"

for(i in 1:length(marker.list)){
    p <- lapply(marker.list[[i]], function(marker) {
        FeaturePlot(object = object, features = marker,pt.size = 0.5, label=T,
                    reduction = "umap")+
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
object %<>% RenameIdents("0" = "NK cells",
                         "1" = "T cells",
                         "2" = "Ciliated cells",
                         "3" = "Macrophages",
                         "4" = "Stromal/fibroblasts",
                         "5" = "Endothelial cells",
                         "6" = "Endothelial cells",
                         "7" = "Alveolar cells/Distal secretory cells",
                         "8" = "Secretory cells",
                         "9" = "Smooth muscle cells",
                         "10" = "Basal cells",
                         "11" = "B cells",
                         "12" = "Monocytes",
                         "13" = "Secretory cells",
                         "14" = "Monocytes",
                         "15" = "Monocytes",
                         "16" = "Macrophages",
                         "17" = "Lymphatic endothelial cells",
                         "18" = "Chondrocytes",
                         "19" = "Unknown",
                         "20" = "Secretory cells")
object@meta.data$cell.type = Idents(object)
object@meta.data$cell.type %<>% as.character()
object@meta.data$singler1main %<>% as.character()
T_HSC_DC <- grep("T-cells|HSC|DC",object@meta.data[,"singler1main"])
object@meta.data[T_HSC_DC,"cell.type"] = object@meta.data[T_HSC_DC,"singler1main"]
object@meta.data$cell.type %<>% gsub("T cells","CD8+ T-cells",.)
Idents(object) = "cell.type"

object %<>% sortIdent()
object <- AddMetaColor(object = object, label= "cell.type", colors = Singler.colors)
UMAPPlot.1(object, group.by = "cell.type",cols = ExtractMetaColor(object),label = T,
           label.repel = T, pt.size = 0.5,label.size = 4, repel = T,no.legend = T,
           do.print = T,do.return = F,title = "Cell types")

UMAPPlot.1(object, group.by = "cell.type",split.by = "conditions",cols = ExtractMetaColor(object),label = F,label.repel = F, pt.size = 0.5,label.size = 4, repel = T,no.legend = T,do.print = T,do.return = F,border = T,,title = "Cell types")


save(object,file="data/Lung_24_20190824.Rda")

#========== P+D+T samples combined (a separate plot with colors based on type of sample - P, D, T)
# All cell types
Idents(object) = "group"
object %<>% subset(idents ="UNC-44", invert = T)
Idents(object) = "cell.type"
for(con in c("distal","proximal","terminal")){
    cellUse = object$conditions %in% con
    sub_obj <- object[,cellUse]
    p <- UMAPPlot.1(object = sub_obj, label = F,label.repel = F, group.by = "cell.type", 
                    cols = ExtractMetaColor(object),
                    do.return = T, no.legend = F, title = paste("UMAP plot for all cell types in",con),
                    pt.size = 0.2,alpha = 1, label.size = 5, do.print = F,unique.name = T)
    jpeg(paste0(path,"UMAPPlot_cell.type_",con,".jpeg"), 
         units='in', width=10, height=7,res=600)
    print(p)
    dev.off()
}

Idents(object) = "integrated_snn_res.0.6"
for(con in c("distal","proximal","terminal")){
    cellUse = object$conditions %in% con
    sub_obj <- object[,cellUse]
    p <- UMAPPlot.1(object = sub_obj, label = T,label.repel = T, group.by = "integrated_snn_res.0.6", 
                    cols = ExtractMetaColor(object),
                    do.return = T, no.legend = T, title = paste("UMAP plot for all clusters in",con),
                    pt.size = 0.2,alpha = 1, label.size = 5, do.print = F,unique.name = T)
    jpeg(paste0(path,"UMAPPlot_clusters_",con,".jpeg"), 
         units='in', width=10, height=7,res=600)
    print(p)
    dev.off()
}
table(object$integrated_snn_res.0.6,object$conditions)
table(object$cell.type,object$conditions)

#================
samples = c("proximal","distal","terminal")
# proximal
args =1
(con <- samples[args])
(load(file = paste0("data/Lung_24",con,"_20190918.Rda")))
table(Idents(object))
Idents(object) = "integrated_snn_res.0.8"
UMAPPlot.1(object = object, label = T,label.repel = T, group.by = "integrated_snn_res.0.8", 
           cols = ExtractMetaColor(object),
           do.return = T, no.legend = T, title = paste("UMAP plot for all clusters in",con),
           pt.size = 0.2,alpha = 1, label.size = 5, do.print = F,unique.name = T)
object %<>% RenameIdents("0" = "T cells:CD4+",
                         "1" = "T cells:TRM",
                         "2" = "Intermediate mucous-like",
                         "3" = "Mucous gland cells",
                         "4" = "Goblet",
                         "5" = "Macrophages",
                         "6" = "B cells",
                         "7" = "T cells:TRM",
                         "8" = "Endothelial cells",
                         "9" = "Ciliated cells",
                         "10" = "Vascular Smooth muscle cells",
                         "11" = "T cells:7SK.2+",
                         "12" = "Endothelial cells",
                         "13" = "Basal cells",
                         "14" = "Fibroblasts",
                         "15" = "Neutrophil",
                         "16" = "Cartilage",
                         "17" = "Mucus-producing goblet cells",
                         "18" = "Mast cells",
                         "19" = "Myo-epithelial cells",
                         "20" = "Pericytes",
                         "21" = "NK cells",
                         "22" = "Proliferating basal cells",
                         "23" = "T cells:CD4+",
                         "24" = "Ionocytes/NEC",
                         "25" = "Neurons",
                         "26" = "Lymphatic endothelial cells",
                         "27" = "NK-like cells")
object@meta.data$manual = as.character(Idents(object))
object@meta.data = cbind(object@meta.data,object@reductions$umap@cell.embeddings)
plasma <- object@meta.data$UMAP_1 >7 & object@meta.data$UMAP_2 < -10
object@meta.data[plasma,"manual"] = "B cells:Plasma"
pre_ciliated <- object$integrated_snn_res.0.8==9 & object@meta.data$UMAP_2 < 14.25 & object@meta.data$UMAP_1 > 1.5
object@meta.data[pre_ciliated,"manual"] = "Pre-ciliated cells"
basal_like <- object$integrated_snn_res.0.8==11 &  object@meta.data$UMAP_1 < 2.2 & object@meta.data$UMAP_2 > -2.2
object@meta.data[basal_like,"manual"] = "Basal-like cells"
ASM <- object$integrated_snn_res.0.8==19 &  object@meta.data$UMAP_1 < -7.8
object@meta.data[ASM,"manual"] = "Airway smooth muscle"

(load(file = "data/Epi_23proximal_20190904.Rda"))
serous <- colnames(Epi)[Epi$integrated_snn_res.0.8 == 1]
object@meta.data[serous,"manual"] = "Serous gland cells"
ESML <- object$integrated_snn_res.0.8==10 &  object@meta.data$UMAP_2 < -5
object@meta.data[ESML,"manual"] = "Endothelial smooth muscle-like cells"
prlif_mf <- object$integrated_snn_res.0.8==5 &  object@meta.data$UMAP_1 < 7.2
object@meta.data[prlif_mf,"manual"] = "Prolifereating macrophage"
monocytes <- object$integrated_snn_res.0.8==5 &  object@meta.data$UMAP_2 < 2
object@meta.data[monocytes,"manual"] = "Monocytes"
DC <- object$integrated_snn_res.0.8==5 &  object@meta.data$UMAP_1 >10
object@meta.data[DC,"manual"] = "Dendritic cells"

UMAPPlot.1(object = object, label = F,label.repel = F, group.by = "manual", 
           cols = Singler.colors, legend.size = 10,
           do.return = T, no.legend = F, title = "UMAP plot for all cell types in proximal",
           pt.size = 0.2,alpha = 0.85, label.size = 3, do.print = T,unique.name = "conditions")
table(object$manual,object$orig.ident) %>% kable() %>% kable_styling()
# Calculate average nUMI in each cell
table(object$manual,object$) %>% kable() %>% kable_styling()


Idents(object) = "orig.ident"
# distal=============
args =2
(con <- samples[args])
(load(file = paste0("data/Lung_24",con,"_20190918.Rda")))
table(Idents(object))

UMAPPlot.1(object = object, label = T,label.repel = T, group.by = "integrated_snn_res.0.8", 
           cols = ExtractMetaColor(object),
           do.return = T, no.legend = T, title = paste("UMAP plot for all clusters in",con),
           pt.size = 0.2,alpha = 1, label.size = 5, do.print = F,unique.name = T)
object %<>% RenameIdents("0" = "fibroblasts",
                         "1" = "T cells:CD8+",
                         "2" = "Ciliated cells",
                         "3" = "Endothelial cells",
                         "4" = "Endothelial cells",
                         "5" = "Distal secretory",
                         "6" = "Ciliated cells",
                         "7" = "Alveolar type 2",
                         "8" = "T cells:CD8+",
                         "9" = "Alveolar macrophages",
                         "10" = "Smooth muscle",
                         "11" = "T cells:CD4+",
                         "12" = "Alveolar macrophages",
                         "13" = "Basal cells",
                         "14" = "Alveolar macrophages",
                         "15" = "B cells",
                         "16" = "Mast cells",
                         "17" = "Smooth muscle",
                         "18" = "B cells:Plasma",
                         "19" = "Alveolar macrophages",
                         "20" = "Lymphatic Endothelial cells",
                         "21" = "Pre-ciliated cells",
                         "22" = "Endothelial cells",
                         "23" = "Alveolar type 1",
                         "24" = "Alveolar macrophages",
                         "25" = "fibroblasts")
object@meta.data$manual = Idents(object)
UMAPPlot.1(object = object, label = T,label.repel = T, group.by = "manual", 
           cols = Singler.colors,
           do.return = T, no.legend = T, title = "UMAP plot for all cell types in proximal",
           pt.size = 0.2,alpha = 0.85, label.size = 4, do.print = T,unique.name = "conditions")

# terminal=============
args =3
(con <- samples[args])
(load(file = paste0("data/Lung_24",con,"_20190918.Rda")))
table(Idents(object))

UMAPPlot.1(object = object, label = T,label.repel = T, group.by = "integrated_snn_res.0.8", 
           cols = ExtractMetaColor(object),
           do.return = T, no.legend = T, title = paste("UMAP plot for all clusters in",con),
           pt.size = 0.2,alpha = 1, label.size = 5, do.print = F,unique.name = T)
object %<>% RenameIdents("0" = "T cells:CD8+",
                         "1" = "Alveolar macrophages",
                         "2" = "Alveolar type 2",
                         "3" = "Endothelial cells",
                         "4" = "Alveolar macrophages",
                         "5" = "NK cells",
                         "6" = "Alveolar macrophages",
                         "7" = "Ciliated cells",
                         "8" = "Distal secretory SCGB3A2+ SFTPB+",
                         "9" = "Smooth muscle",
                         "10" = "B cells", 
                         "11" = "Fibroblasts",
                         "12" = "Alveolar macrophages",
                         "13" = "Endothelial cells",
                         "14" = "Basal cells",
                         "15" = "T cells:CD8+",
                         "16" = "Alveolar macrophages",
                         "17" = "Mast cells",
                         "18" = "Endothelial cells",
                         "19" = "Endothelial cells",
                         "20" = "Lymphatic Endothelial cells",
                         "21" = "B cells:Plasma")
object@meta.data$manual = as.character(Idents(object))
object@meta.data = cbind(object@meta.data,object@reductions$umap@cell.embeddings)
basal <- object@meta.data$manual %in% "Basal cells"
hist(object@meta.data[basal,"UMAP_1"],breaks = 30)
type1 <- object@meta.data$UMAP_1 < -1 & object@meta.data$UMAP_2 >6 
object@meta.data[(basal & type1),"manual"] = "Alveolar type 1"

UMAPPlot.1(object = object, label = T,label.repel = T, group.by = "manual", 
           cols = Singler.colors,
           do.return = T, no.legend = T, title = "UMAP plot for all cell types in proximal",
           pt.size = 0.2,alpha = 0.85, label.size = 4, do.print = T,unique.name = "conditions")

# samples

#number of expressed genes in each cluster (per sample)
samples = c("proximal")
for(con in samples){
    (load(file = paste0("data/Lung_24",con,"_20190918.Rda")))
    table(object$manual, object$orig.ident) %>% as.data.frame() %>% spread(Var2,Freq)
    df <- table(object$manual, object$orig.ident) %>% as.data.frame()
    colnames(df) =c("cell.types","samples","nGene")
    
    df$cell.types <- as.character(df$cell.types)
    df$samples <- as.character(df$samples)
    
    for(i in seq_len(nrow(df))) {
        cells <- object$manual %in% df[i,"cell.types"] & object$orig.ident %in% df[i,"samples"]
        df[i,"nGene"] = as.integer(mean(object@meta.data[cells,"nFeature_RNA"]))
        svMisc::progress(i/nrow(df)*100)
    }
    df %<>% spread("samples","nGene")
    write.csv(df,paste0(path,"Lung_24-",con,"_nGene_by_samples.csv"))
}

samples = c("distal","terminal")
for(con in samples){
    (load(file = paste0("data/Lung_24",con,"_20190918.Rda")))
    df <- table(object$integrated_snn_res.0.8, object$orig.ident) %>% as.data.frame()
    colnames(df) =c("cluster_res.0.8","samples","nGene")
    
    df$cluster_res.0.8 <- as.character(df$cluster_res.0.8)
    df$samples <- as.character(df$samples)
    
    for(i in seq_len(nrow(df))) {
        cells <- object$integrated_snn_res.0.8 %in% df[i,"cluster_res.0.8"] & object$orig.ident %in% df[i,"samples"]
        df[i,"nGene"] = as.integer(mean(object@meta.data[cells,"nFeature_RNA"]))
        svMisc::progress(i/nrow(df)*100)
    }
    df %<>% spread("samples","nGene")
    df[is.na(df)] = 0
    write.csv(df,paste0(path,"Lung_24-",con,"_nGene_by_samples.csv"))
}
