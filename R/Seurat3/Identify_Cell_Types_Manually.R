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
(load(file="data/Lung_24_20190824.Rda"))
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
DefaultAssay(object) <- 'SCT'
Idents(object) <- "integrated_snn_res.0.6"

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

# Epi cell types
#====== 2.2 Identify Epi cell types ==========================================
(load(file="data/Lung_24_20190824.Rda"))

df_markers <- readxl::read_excel("doc/Renat.markers.xlsx",sheet = "20190613")

#df_markers <- readxl::read_excel("../seurat_resources/bio-rad-markers.xlsx",sheet = "Human.sub")
colnames(df_markers) = gsub(" ","_",colnames(df_markers))
colnames(df_markers) = gsub(":|\\/","_",colnames(df_markers))
colnames(df_markers) = gsub("\\+","",colnames(df_markers))
#markers = df_markers[,-grep("Alias",colnames(df_markers))]
marker.list <- df2list(df_markers)

marker.list %<>% lapply(function(x) x) %>% 
    lapply(function(x) FilterGenes(Epi,x)) %>% 
    lapply(function(x) x[!is.na(x)])
#marker.list %>% list2df %>% t %>% kable() %>% kable_styling()
Idents(Epi) <- "RNA_snn_res.0.3"

for(i in 1:length(marker.list)){
    p <- lapply(marker.list[[i]], function(marker) {
        FeaturePlot(object = Epi, features = marker,pt.size = 0.5, label=T)+
            NoLegend()+
            ggtitle(paste0(marker,Alias(df = df_markers,gene = marker)))+
            theme(plot.title = element_text(hjust = 0.5,size = 15,face = "plain"))
    })
    jpeg(paste0(path,names(marker.list)[i],".jpeg"),units="in", width=10, height=7,res=600)
    print(do.call(cowplot::plot_grid, p)+ ggtitle(paste(names(marker.list)[i],"markers"))+
              theme(plot.title = element_text(hjust = 0.5,size = 20)))
    dev.off()
    print(paste0(i,":",length(marker.list)))
}

#======== rename ident =================
Epi$cell.type <- plyr::mapvalues(Epi$RNA_snn_res.0.8,
                                    from = 1:9,
                                    to = c("Alveolar type II cells",
                                           "M-ciliated cells",
                                           "Airway secretory cells",
                                           "D-ciliated cells",
                                           "Alveolar type I cells",
                                           "Distal airway secretory cells",
                                           "Airway basal stem cells",
                                           "Unknown Alveolar cells",
                                           "Alveolar type II cells"))
Idents(Epi) = "cell.type"
Epi %<>% sortIdent()
table(Idents(Epi)) %>% kable() %>% kable_styling()
# process color scheme
singler_colors <- readxl::read_excel("doc/singler.colors.xlsx")
singler_colors1 = as.vector(singler_colors$singler.color1[!is.na(singler_colors$singler.color1)])
singler_colors1[duplicated(singler_colors1)]
length(singler_colors1)
Epi <- AddMetaColor(object = Epi, label= "cell.type", colors = singler_colors1)
Epi %<>% sortIdent()

TSNEPlot.1(Epi, cols = ExtractMetaColor(Epi),label = F,pt.size = 1,
           label.size = 5, repel = T,no.legend = F,do.print = T,
           title = "Epithelial cell types")

p5 <- UMAPPlot(Epi, group.by="cell.type",pt.size = 1,label = F,
               cols = ExtractMetaColor(Epi),
               label.size = 4, repel = T)+ggtitle("Epithelial cell types")+
    theme(plot.title = element_text(hjust = 0.5,size=15,face = "plain"))

jpeg(paste0(path,"Epi_umap_cell.type.jpeg"), units="in", width=10, height=7,res=600)
print(p5)
dev.off()
save(Epi, file = "data/Lung_24_20190824.Rda")
