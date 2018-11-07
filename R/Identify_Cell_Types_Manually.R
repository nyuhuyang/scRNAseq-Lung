library(Seurat)
library(SingleR)
library(dplyr)
library(tidyr)
library(kableExtra)
source("../R/Seurat_functions.R")
source("../R/SingleR_functions.R")
path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)
#====== 2.1 pathway analysis ==========================================
lnames = load(file = "data/Lung_20181101.Rda")
list_files = list.files("doc/GeneSets")
gene_set_list <- lapply(paste0("doc/GeneSets/",list_files), function(x){
    read.delim(x, stringsAsFactors = F)[,1]})
names(gene_set_list) <- gsub(".txt","",list_files)
gene_set_list <- lapply(gene_set_list, function(x) HumanGenes(Lung,x))
gene_set_df <- list2df(gene_set_list)
str(gene_set_list)

for(i in 1:length(gene_set_list)){
    Lung <- .AddModuleScore(Lung, genes.list = gene_set_list[i],
                                   ctrl.size = 5,enrich.name = names(gene_set_list[i]))
}

Colmin <- apply(Lung@meta.data[,names(gene_set_list)], 2, min)
Lung@meta.data[,names(gene_set_list)] <- sweep(Lung@meta.data[,names(gene_set_list)], 2, Colmin,"-")

for(j in 1:length(gene_set_list)){
    p <- SingleFeaturePlot.1(object = Lung,
                              feature = names(gene_set_list[j]), threshold = 0.3)
    jpeg(paste0(path,names(gene_set_list[j]),".jpeg"), units="in", width=10, height=7,res=600)
    print(p)
    dev.off()
}

#====== 2.2 marker gene analysis ==========================================
# combine marker gene table
immgen_main = read.csv("../SingleR/output/immgen_main.csv",row.names =1,header = T,
                      stringsAsFactors = F)
Renat.markers <- readxl::read_excel("doc/Renat.markers.xlsx")
colnames(Renat.markers) <- gsub("\\s","_",colnames(Renat.markers))

immgen_Renat <- merge(t(Renat.markers), t(immgen_main),by="row.names",all=TRUE)
rownames(immgen_Renat) <- immgen_Renat$Row.names
immgen_Renat <- immgen_Renat[,-which(colnames(immgen_Renat) == "Row.names")]
immgen_Renat.list <- df2list(t(immgen_Renat))


marker.list <- lapply(immgen_Renat.list, function(x) {
    HumanGenes(object = Lung, marker.genes= x, unique = T)
})
marker.df <- list2df(marker.list)
marker.list <- df2list(marker.df[1:9,])
str(marker.list)

FeaturePlot.1 <- function(object = Lung, x){
    p <- FeaturePlot(object = object, 
                     reduction.use = "tsne",
                     features.plot = x, min.cutoff = NA, do.return =T,
                     cols.use = c("lightgrey","blue"), pt.size = 0.5)
    return(p)
}

dev.off() # to accelerate 
for(i in 1:length(marker.list)){
    p <- FeaturePlot.1(object = Lung, x = marker.list[[i]])
    p1 <- do.call(plot_grid, p)
    p1 <- p1 + ggtitle(paste(names(marker.list)[i],"markers"))+
        theme(plot.title = element_text(hjust = 0.5,size = 20, face = "bold"))
    jpeg(paste0(path,names(marker.list)[i],".jpeg"),
         units="in", width=10, height=7,res=600)
    print(p1)
    print(paste0(i,":",length(marker.list)))
    dev.off()
}

CD14_Monocytes <-  HumanGenes(Lung,c("CD14","LYZ","S100A9","CCL2","CCR2"))
CD16_Monocytes <- HumanGenes(Lung,c("FCGR3A","MS4A7","VMO1","CCR2"))

marker.list <- list("Monocytes" = c(CD14_Monocytes,CD16_Monocytes))
#====== 2.3 ==========================================
Lung <- SetAllIdent(Lung, id = "res.0.8")
table(Lung@ident)
idents <- as.data.frame(table(Lung@ident))
old.ident.ids <- idents$Var1
new.cluster.ids <- c("Macrophages",
                     "Endothelial cells",
                     "Alveolar type II cells",
                     "Stromal/fibroblast",
                     "Neutrophils",
                     "T/B/NK cells",
                     "Macrophages",
                     "Smooth muscle cells",
                     "Mucus-producing cells\nEpithelial cells",
                     "Basal cells",
                     "Ciliated cells",
                     "Red blood cells",
                     "Distal secretory cells",
                     "basophils/Mast cells",
                     "Macrophages",
                     "Lymphatic endothelial cells",
                     "Macrophages")
Lung@ident <- plyr::mapvalues(x = Lung@ident,
                              from = old.ident.ids,
                              to = new.cluster.ids)
jpeg(paste0(path,"TSNEplot.jpeg"),
     units="in", width=10, height=7,res=600)
TSNEPlot.1(object = Lung,do.label = T, group.by = "ident", 
           do.return = TRUE, no.legend = T, colors.use = singler.colors,
           pt.size = 2,label.size = 5,label.repel = T)+
    ggtitle("Manually lable all cell types")+
    theme(text = element_text(size=20),							
          plot.title = element_text(hjust = 0.5, face = "bold"))
dev.off()

Lung <- StashIdent(Lung, save.name = "manual")

#' @param df Nx2 data.frame
#' @param label group label
#' @export majority mejority element in each group
FindMajority <- function(df, label) {
    levels = sort(unique(df[,label]))
    majority <- lapply(levels, function(x) {
        names(which.max(table(df[(df[,label] %in% x),
                                 -which(colnames(df) == label)])))
        })
    majority <- unlist(majority)
    names(majority) <- levels
    
    return(majority)
}

colors <- FindMajority(Lung@meta.data[,c("singler2sub.colors","manual")], "manual")
colors[c('Alveolar type II cells','Basal cells','Ciliated cells','Distal secretory cells')]
Lung@ident <- factor(Lung@ident,level = sort(unique(new.cluster.ids)))

singler_colors <- readxl::read_excel("./doc/singler.colors.xlsx")
singler_colors = as.vector(singler_colors$singler.color1[!is.na(singler_colors$singler.color1)])
length(singler_colors)

Lung <- AddMetaColor(object = Lung, colors = singler_colors)
df_colors <- ExtractMetaColor(object = Lung)
df_colors

p3 <- TSNEPlot.1(object = Lung, dim.1 = 1, dim.2 = 2, 
                 group.by = "ident", do.return = TRUE,pt.size = 2,
                 colors.use = ExtractMetaColor(Lung),no.legend = F,
                 do.label =F,label.size=4, label.repel = T,force=1)+
    ggtitle("Manually label cell types")+
    theme(text = element_text(size=15),
          plot.title = element_text(hjust = 0.5))
jpeg(paste0(path,"PlotTsne_manual.jpeg"), units="in", width=10, height=7,
     res=600)
p3
dev.off()


unknow.cell  <- FeaturePlot(Lung, features.plot = "CLDN5",cols.use = c("lightgrey","blue"),
            do.identify = T)
New <- SubsetData(Lung, cells.use = unknow.cell)
TSNEPlot.1(New)
