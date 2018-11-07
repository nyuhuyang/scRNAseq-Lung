library(SingleR)
library(Seurat)
library(dplyr)
library(reshape2)
library(pheatmap)
library(kableExtra)

source("../R/Seurat_functions.R")
source("../R/SingleR_functions.R")
#====== 3.1 Create Singler Object  ==========================================
lname1 = load(file = "./data/Lung_20181101.Rda")

singler = CreateSinglerObject(as.matrix(Lung@data), annot = NULL, project.name=Lung@project.name,
                              min.genes = 200,technology = "10X", species = "Human", citation = "",
                              normalize.gene.length = F, variable.genes = "de",
                              fine.tune = F, do.signatures = F, clusters = NULL)
# Did singler find all cell labels?
length(singler$singler[[1]]$SingleR.single$labels) == ncol(Lung@data)
singler$meta.data$orig.ident = Lung@meta.data$orig.ident # the original identities, if not supplied in 'annot'
singler$meta.data$xy = Lung@dr$tsne@cell.embeddings # the tSNE coordinates
singler$meta.data$clusters = Lung@ident # the Seurat clusters (if 'clusters' not provided)
save(singler,file="./data/singler_Lung20181101.RData")
#====== 3.2 SingleR specifications ==========================================
# Step 1: Spearman coefficient
lnames = load(file = "./data/singler_Lung20181101.RData")
lnames
SingleR.DrawScatter(sc_data = singler$seurat@data,cell_id = 10, 
                    ref = immgen, sample_id = 232)

# Step 2: Multiple correlation coefficients per cell types are aggregated 
# to provide a single value per cell type per single-cell. 
# In the examples below we use the 80% percentile of correlation values.
# for visualization purposes we only present a subset of cell types (defined in labels.use)
out = SingleR.DrawBoxPlot(sc_data = singler$seurat@data,cell_id = 10, 
                          ref = immgen,main_types = T,
                          labels.use=c('B cells','T cells','DC','Macrophages','Monocytes','NK cells',
                                       'Mast cells','Neutrophils','Fibroblasts','Endothelial cells'))
print(out$plot)

##############################
# Human Primary Cell Atlas (HPCA)
###############################
SingleR.DrawHeatmap(singler$singler[[1]]$SingleR.single.main, top.n = Inf,
                    clusters = singler$meta.data$orig.ident)
#Or by all cell types (showing the top 50 cell types):
jpeg(paste0(path,"/DrawHeatmap.jpeg"), units="in", width=10, height=7,
     res=600)
print(SingleR.DrawHeatmap(singler$singler[[2]]$SingleR.single, top.n = 50))
dev.off()
jpeg(paste0(path,"/DrawHeatmap_normF.jpeg"), units="in", width=10, height=7,
     res=600)
print(SingleR.DrawHeatmap(singler$singler[[2]]$SingleR.single,top.n = 50,normalize = F))
dev.off()
#Next, we can use the fine-tuned labels to color the t-SNE plot:
       
out = SingleR.PlotTsne.1(singler$singler[[2]]$SingleR.single,
                       singler$meta.data$xy,do.label=T,
                       do.letters = F,labels = singler$singler[[2]]$SingleR.single$labels,
                       label.size = 4, label.repel = T,dot.size = 3,do.legend = F,alpha = 1,
                       force=2)
jpeg(paste0(path,"/PlotTsne_sub2.jpeg"), units="in", width=10, height=7,
     res=600)
out+  ggtitle("Supervised sub-cell type labeling by Blueprint+ Encode")+
        theme(text = element_text(size=20),
              plot.title = element_text(hjust = 0.5,size = 18, face = "bold"))
dev.off()
# main types-------
out = SingleR.PlotTsne.1(singler$singler[[1]]$SingleR.single.main,
                         singler$meta.data$xy,do.label=T,
                         do.letters = F,labels = singler$singler[[1]]$SingleR.single.main$labels,
                         label.size = 5, dot.size = 3,do.legend = F,alpha = 1,
                         label.repel = T,force=2)
jpeg(paste0(path,"/PlotTsne_main1.jpeg"), units="in", width=10, height=7,
     res=600)
out+  ggtitle("Supervised main-cell type labeling by HPCA")+
        theme(text = element_text(size=20),
              plot.title = element_text(hjust = 0.5,size = 18, face = "bold"))
dev.off()
#Finally, we can also view the labeling as a table compared to the original identities:

kable(table(singler$singler[[1]]$SingleR.single$labels,
            singler$meta.data$orig.ident)) %>%
        kable_styling()
kable(table(singler$meta.data$orig.ident,
            singler$singler[[1]]$SingleR.single.main$labels)) %>%
        kable_styling()
kable(table(singler$meta.data$orig.ident,singler$seurat@ident)) %>%
        kable_styling()

##############################
# add singleR label to Seurat
###############################
singlerDF = data.frame("singler1sub"=singler$singler[[1]]$SingleR.single$labels,
                       "singler1main"=singler$singler[[1]]$SingleR.single.main$labels,
                       "singler2sub"=singler$singler[[2]]$SingleR.single$labels,
                       "singler2main"=singler$singler[[2]]$SingleR.single.main$labels,
                       "cell.names" = rownames(singler$singler[[1]]$SingleR.single$labels))
knowDF = data.frame("cell.names"= Lung@cell.names)
ident.DF = full_join(singlerDF,knowDF, by="cell.names")
ident.DF<- apply(ident.DF,2,as.character)
rownames(ident.DF) = ident.DF[,"cell.names"]
ident.DF = ident.DF[,-which(colnames(ident.DF) == "cell.names")]
apply(ident.DF, 2, function(x) length(unique(x))) # check number of labels
#ident.DF[is.na(ident.DF)] <- "unknown"
Lung <- AddMetaData(object = Lung,
                   metadata = as.data.frame(ident.DF))
Lung <- SetAllIdent(object = Lung, id = "singler2sub")
##############################
# process singler.color
##############################
singler_colors <- readxl::read_excel("./doc/singler.colors.xlsx")
singler_colors = as.vector(singler_colors$singler.color1[!is.na(singler_colors$singler.color1)])
length(singler_colors)
length(unique(singler$singler[[1]]$SingleR.single.main$labels))
sort(unique(singler$singler[[1]]$SingleR.single$labels)) %>%
        kable() %>%
        kable_styling()

Lung <- SetAllIdent(Lung, id= "singler1main")
Lung <- AddMetaColor(object = Lung, colors = singler_colors[3:27])

p3 <- TSNEPlot.1(object = Lung, dim.1 = 1, dim.2 = 2, 
                group.by = "ident", do.return = TRUE,pt.size = 2,
                colors.use = ExtractMetaColor(Lung),no.legend = T,
                do.label =T,label.size=4, label.repel = T,force=1)+
        ggtitle("Supervised sub cell type labeling by HPCA")+
        theme(text = element_text(size=15,face ="bold"),
              plot.title = element_text(hjust = 0.5))
jpeg(paste0(path,"PlotTsne_main1.jpeg"), units="in", width=10, height=7,
     res=600)
p3# +xlim(-50, 40)
dev.off()

save(Lung,file="./output/Lung_alignment20181031.Rda")
##############################
# subset Seurat
###############################
head(Lung@meta.data)
table(Lung@meta.data$orig.ident)
table(Lung@ident)
g <- SplitTSNEPlot(Lung,group.by = "ident",split.by = "projects",
                   no.legend = T,do.label =T,label.size=3,pt.size = 2,
                   return.plots =T, label.repel = T,force=2)

jpeg(paste0(path,"SplitTSNEPlot.jpeg"), units="in", width=10, height=7,
     res=600)
print(do.call(plot_grid, g))
dev.off()
#---------------
g1 <- SplitTSNEPlot(Lung,group.by = "ident",split.by = "orig.ident",
                   no.legend = T,do.label =T,label.size=3,pt.size = 2,
                   return.plots =T, label.repel = T,force=2)
splited.Lung <- SplitSeurat(Lung,split.by = "orig.ident")
(samples <- splited.Lung[[length(splited.Lung)]])
for(i in 1:length(g1)){
        jpeg(paste0(path,"TSNEPlot_",samples[i],".jpeg"), units="in", width=10, height=7,
             res=600)
        print(g1[[i]])
        print(paste0(i,":",length(g1)))
        dev.off()
}
#---------------
g2 <- SplitTSNEPlot(Lung,group.by = "ident",split.by = "orig.ident",
                    no.legend = T,do.label =F,label.size=3,pt.size = 2,
                    size =15,return.plots =T, label.repel = T,force=2)
jpeg(paste0(path,"splited_TSNEPlot.jpeg"), units="in", width=10, height=7,
     res=600)
print(do.call(plot_grid,g2))
dev.off()

save(Lung, file = "./data/Lung_20181101.Rda")
