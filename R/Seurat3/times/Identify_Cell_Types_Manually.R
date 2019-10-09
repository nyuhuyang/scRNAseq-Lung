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
(load(file="data/Lung_8_20190808.Rda"))
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
Idents(object) <- "integrated_snn_res.1.2"

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
"Basal cells",



#======== rename ident =================
object %<>% RenameIdents("0" = "Secretory cells",
                         "1" = "Distal secretory cells",
                         "2" = "Ciliated cells",
                         "3" = "Basal cells",
                         "4" = "Secretory cells",
                         "5" = "Basal cells",
                         "6" = "Ciliated cells",
                         "7" = "Ciliated cells",
                         "8" = "Basal cells",
                         "9" = "Distal secretory cells",
                         "10" = "Secretory cells",
                         "11" = "Secretory cells",
                         "12" = "Distal secretory cells",
                         "13" = "Secretory cells",
                         "14" = "Distal secretory cells",
                         "15" = "Ciliated cells",
                         "16" = "Ciliated cells",
                         "17" = "Ciliated cells",
                         "18" = "Ciliated cells",
                         "19" = "Basal cells",
                         "20" = "Mast cells",
                         "21" = "Ciliated cells",
                         "22" = "Distal secretory cells")


object[['cell.type']] = Idents(object)

TSNEPlot.1(object, cols = ExtractMetaColor(object),label = T,pt.size = 1,
           label.size = 5, repel = T,no.legend = F,do.print = F,
           title = "Cell types")
##############################
# process color scheme
##############################
object <- sortIdent(object)
table(Idents(object)) %>% kable %>% kable_styling()
singler_colors <- readxl::read_excel("doc/singler.colors.xlsx")
singler_colors1 = as.vector(singler_colors$singler.color1[!is.na(singler_colors$singler.color1)])
singler_colors1[duplicated(singler_colors1)]
length(singler_colors1)
length(unique(object$cell.type))
object <- AddMetaColor(object = object, label= "cell.type", colors = SingleR::singler.colors)
Idents(object) <- "cell.type"

object <- sortIdent(object)

TSNEPlot.1(object, cols = ExtractMetaColor(object),label = T,pt.size = 1,label.repel = T,
           label.size = 5, repel = T,no.legend = F,do.print = T,
           title = "Cell types")
UMAPPlot.1(object, cols = ExtractMetaColor(object),label = T,pt.size = 1,label.repel = T,
           label.size = 5, repel = T,no.legend = F,do.print = T,
           title = "Cell types")
save(object,file="data/Lung_8_20190808.Rda")

Idents(object) = "orig.ident"
(samples = c("All","Day-0","Day-3","Day-7","Day-14","Day-21","Day-28",
             "Day-56","Day-122"))
for(sample in samples[-1]){
        sub_object <- subset(object, idents = (if(sample == "All") samples[-1] else sample))
        meta.data = cbind.data.frame(sub_object@meta.data,
                                     sub_object@reductions$tsne@cell.embeddings,
                                     sub_object@reductions$umap@cell.embeddings)
        meta.data = meta.data[,c("tSNE_1","tSNE_2","UMAP_1","UMAP_2","integrated_snn_res.1.2")]
        meta.data$integrated_snn_res.1.2 = as.numeric(as.character(meta.data$integrated_snn_res.1.2))
        
        meta.data = meta.data[order(meta.data$integrated_snn_res.1.2),]
        print(colnames(meta.data))
        write.csv(meta.data, paste0(path,sample,"_tSNE_UMAP_coordinates.csv"))
        
        data = as.matrix(DelayedArray::t(sub_object@assays$SCT@data))
        tsne_data = cbind(meta.data[,3:5], data[match(rownames(meta.data),
                                                      rownames(data)), ])
        
        tsne_data = tsne_data[,-1:2]
        write.csv(DelayedArray::t(tsne_data), paste0(path,sample,"_Expression_data.csv"))
}



#====== 2.2 Identify Epi cell types ==========================================
(load(file="data/Epi_harmony_12_20190627.Rda"))

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