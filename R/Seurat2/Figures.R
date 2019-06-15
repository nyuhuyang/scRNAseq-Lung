library(SingleR)
library(Seurat)
library(reshape2)
library(pheatmap)
library(kableExtra)
library(dplyr)
library(tidyr)
library(magrittr)
library(readxl)
library(gplots)
library(ggrepel)
source("../R/Seurat_functions.R")
source("../R/SingleR_functions.R")
source("R/util.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)

#=== load data ======
(load(file="data/Lung_MNN_9_20181101.Rda"))
save(object,file = "data/Lung_MNN_9_20190405.Rda")
##############################
# Dot plot
##############################
object %<>% SetAllIdent("res.0.8")
object = RenameIdent(object, old.ident.name = 16, new.ident.name = 6)
object = StashIdent(object, "res.0.8")

tsne_data <- data.frame("tSNE_1" =object@dr$tsne@cell.embeddings[,"tSNE_1"],
                        "tSNE_2" =object@dr$tsne@cell.embeddings[,"tSNE_2"],
                        "Cell_type" = object@meta.data$manual)
write.csv(tsne_data,paste0(path,"tsne_data.csv"))
# process color scheme======
singler_colors <- readxl::read_excel("doc/singler.colors.xlsx")
singler_colors1 = as.vector(singler_colors$singler.color1[!is.na(singler_colors$singler.color1)])
singler_colors2 = as.vector(singler_colors$singler.color2[!is.na(singler_colors$singler.color2)])
singler_colors1[duplicated(singler_colors1)]
singler_colors2[duplicated(singler_colors2)]
length(singler_colors1);length(singler_colors2)
apply(object@meta.data[,c("manual","res.0.8")],
      2,function(x) length(unique(x)))
object <- AddMetaColor(object = object, label= "res.0.8", colors = singler_colors2)
object <- SetAllIdent(object = object, id = "res.0.8")

# TSEN plot====
TSNEPlot.1(object,label.repel = F,do.label = T,do.print = F,label.size = 5,
           colors.use = ExtractMetaColor(object),no.axes=T)
TSNEPlot.1(object,do.print = T,label.size = 5,no.legend = F, do.label = F)

# expression data

meta.data = t(object@meta.data[,c("orig.ident","res.0.8","projects","manual")])
data = object@data; format(object.size(data),unit="GB")
data = as.data.frame(as.matrix(data));format(object.size(data),unit="GB")
data = sapply(data,as.character)
rownames(data) = rownames(object@data)
table(colnames(meta.data) == colnames(data))
exp = rbind.data.frame(meta.data,data)
format(object.size(exp),unit="GB")

data[1:10,1:2]
exp[1:10,1:2]
exp = exp[,order(exp["res.0.8",])]
rownames(exp)[1:4] = c("samples","cluster","sequencing method", "cell type")
write.csv(exp,paste0(path,"lung_exp.csv"))
# rename ident, Subset epithelial ========
table(object@ident)
idents <- as.data.frame(table(object@ident))
(old.ident.ids <- idents$Var1)
(new.ident.ids = as.vector(singler_colors$`20181101_res.0.8_manual`))
object@ident <- plyr::mapvalues(x = object@ident,
                                from = old.ident.ids,
                                to = new.ident.ids)
epithelial_ident <- c("9 Basal cells",
                      "10 Ciliated cells",
                      "8 Mucus-producing cells\nEpithelial cells",
                      "12 Distal secretory cells",
                      "2 Alveolar type II cells")
epithelial <- SubsetData(object, ident.use = epithelial_ident)
epithelial@ident = factor(epithelial@ident, levels = epithelial_ident)
dotplot_markers <- read_excel("doc/Renat.markers.xlsx",sheet = "Dotplot_epithelial")
markers.to.plot <- FilterGenes(object,dotplot_markers$Gene)

jpeg(paste0(path,"Dotplot_epithelial.jpeg"), units="in", width=10, height=7,res=600)
DotPlot.1(epithelial, genes.plot = markers.to.plot,
        cols.use = c("lightgrey","blue"), x.lab.rot = T, plot.legend = T,
        dot.scale = 8, size = 8,do.return = T)
dev.off()

# explore cell type markers ---(skip)---
Hpca_Blueprint_encode_main <- read.csv("../SingleR/output/Hpca_Blueprint_encode_sub.csv",header = 1)
markers_All <- read.csv("output/20190311/markers.All.csv")
markers <- FilterGenes(object, markers_All[markers_All$cluster == 13,"gene"][1:30])
(df <- SearchAllMarkers(Hpca_Blueprint_encode_main,markers))
df[grep("Endothelial_cells",colnames(Hpca_Blueprint_encode_main)),]

####################################
# dotplot all
####################################
(object_ident <- unique(object@ident) %>% sort(decreasing = T))
object@ident = factor(object@ident, levels = object_ident)

dotplot_markers <- read_excel("doc/Renat.markers.xlsx",sheet = "Dotplot_all")
markers.to.plot <- FilterGenes(object,dotplot_markers$Gene)

# a. A subset of most significant differentially expressed markers in the list A: p val adj <0.1^(-200)
markers = read.csv(paste0("output/20190311/markers.All.csv"))

markers.to.plot = unique(markers[(markers$p_val_adj<0.1^300 & markers$avg_logFC >4),"gene"])
length(markers.to.plot)

# b. 19 Lineage TF-list above (see item 2-c2)
TF = c("HMGA1", "TP63", "SOX15", "SOX7", "ELF3", "SPDEF", "HES4", "RARRES1",
       "XBP1", "SOX4", "FOXJ1", "RFX2", "TP73", "MYB", "SOX2", "HOPX", "NKX2-1", "ETV5", "FOXA2")
markers.to.plot <- FilterGenes(object, TF)


# c. FGF-FGFR family genes (see in the attachment) - this is for the second paper
FGF_FGF <- readxl::read_excel("doc/FGF FGFR pathway genes.xlsx")
markers.to.plot <- FilterGenes(object, toupper(FGF_FGFR$gene))

# d. EGF-EGFR family genes (see in the attachment)
EGF_EGFR <- readxl::read_excel("doc/EGF-EGFR family genes.xlsx") 
markers.to.plot <- FilterGenes(object, EGF_EGFR$gene)


jpeg(paste0(path,"d-Dotplot_all.jpeg"), units="in", width=10, height=7,res=600)
DotPlot.1(object, genes.plot = rev(markers.to.plot),
        cols.use = c("lightgrey","blue"), x.lab.rot = T, plot.legend = T,
        dot.scale = 6, size = 8,do.return = T)
dev.off()

####################################
# dotplot epi
####################################
# get DE genes from Seurat
epithelial@meta.data$cluster_manual = paste0(epithelial@meta.data$res.0.8,"_",
                                             epithelial@meta.data$manual)
epithelial %<>% SetAllIdent("cluster_manual")
(unique(epithelial@ident) %>% sort(decreasing = T))
epithelial@ident = factor(epithelial@ident, levels = c("9_Basal cells",
                                                       "8_Mucus-producing cells\nEpithelial cells",
                                                       "12_Distal secretory cells",
                                                       "10_Ciliated cells",
                                                       "2_Alveolar type II cells"))

# a. A subset of most significant differentially expressed markers in the list A: p val adj <0.1^(-200)
markers.Epi = read.csv(paste0("output/20190311/markers.Epi.csv"))
markers.to.plot = unique(markers.Epi[markers.Epi$p_val_adj<0.1^200,"gene"])
length(markers.to.plot)

# b. 19 Lineage TF-list above (see item 2-c2)
TF = c("HMGA1", "TP63", "SOX15", "SOX7", "ELF3", "SPDEF", "HES4", "RARRES1",
       "XBP1", "SOX4", "FOXJ1", "RFX2", "TP73", "MYB", "SOX2", "HOPX", "NKX2-1", "ETV5", "FOXA2")
markers.to.plot <- FilterGenes(epithelial, TF)

# c. FGF-FGFR family genes (see in the attachment) - this is for the second paper
FGF_FGF <- readxl::read_excel("doc/FGF FGFR pathway genes.xlsx")
markers.to.plot <- FilterGenes(epithelial, toupper(FGF_FGFR$gene))

# d. EGF-EGFR family genes (see in the attachment)
EGF_EGFR <- readxl::read_excel("doc/EGF-EGFR family genes.xlsx") 
markers.to.plot <- FilterGenes(epithelial, EGF_EGFR$gene)


jpeg(paste0(path,"d-Dotplot_epi.jpeg"), units="in", width=10, height=7,res=600)
DotPlot.1(epithelial, genes.plot = rev(markers.to.plot),
          cols.use = c("lightgrey","blue"), x.lab.rot = T, plot.legend = T,
          dot.scale = 6, size = 8,do.return = T)
dev.off()

#=== epithelial tsne =========
remove<- FeaturePlot(epithelial,"KRT19",do.identify = T)
write.csv(remove, "output/remove_epi.csv")
epithelial <- SubsetData(epithelial, cells.use = epithelial@cell.names[!epithelial@cell.names %in% remove])
epithelial <- SetAllIdent(object = epithelial, id = "res.0.8")

TSNEPlot.1(epithelial,label.repel = F,do.label = T,do.print = T,label.size = 5,
           colors.use = ExtractMetaColor(epithelial))
TSNEPlot.1(epithelial,do.print = T,label.size = 5,no.legend = F,
           colors.use = ExtractMetaColor(epithelial))
save(epithelial, file = "data/epithelial_Harmony_5_20190123.Rda")


##############################
#  DoHeatmap- all ============
##############################
object <- SubsetData(object, ident.remove = "16 Macrophages")

object@meta.data$heatmap = paste(object@meta.data$res.0.8, object@meta.data$manual)
object %<>% SetAllIdent(id = "heatmap")
table(object@ident)
group.order = c( "9 Basal cells","10 Ciliated cells","8 Mucus-producing cells\nEpithelial cells",
                 "12 Distal secretory cells","2 Alveolar type II cells",
                 "3 Stromal/fibroblast","7 Smooth muscle cells",
                 "1 Endothelial cells","15 Lymphatic endothelial cells","0 Macrophages",
                 "6 Macrophages","14 Macrophages","4 Neutrophils","5 T/B/NK cells",
                 "13 basophils/Mast cells","11 Red blood cells")

# 0. A subset of most significant differentially expressed markers in the list A: p val adj <0.1^(-200)
markers = read.csv(paste0("output/20190311/markers.All.csv"))
markers$order <- plyr::mapvalues(x = markers$cluster,
                                from = c(9,10,8,12,2,3,7,1,15,0,6,14,4,5,13,11),
                                to =   0:15)
markers = markers[order(markers$order),]
markers = markers[markers$cluster != 16,]
markers = markers[(markers$p_val_adj<0.1^50),]
dim(markers)
top <-  markers %>% 
        group_by(cluster) %>% 
        top_n(5000, avg_logFC) %>% 
        as.data.frame
markers.to.plot = top$gene
length(markers.to.plot)

# 2. FGF-FGFR family genes (see in the attachment) - this is for the second paper
FGF_FGF <- readxl::read_excel("doc/FGF FGFR pathway genes.xlsx")
markers.to.plot <- FilterGenes(object, toupper(FGF_FGF$gene))

# 3. EGF-EGFR family genes (see in the attachment)
EGF_EGFR <- readxl::read_excel("doc/EGF-EGFR family genes.xlsx") 
markers.to.plot <- FilterGenes(object, EGF_EGFR$gene)

g <- DoHeatmap(object, genes.use = markers.to.plot,use.scaled = T,group.order = group.order, 
               title= "Top DE genes (adjusted p value < 10^(-50)) in all cell types",
               #title = "FGF-FGFR family genes in all cell types",
               #title = "EGF-EGFR family genes in all cell types",               
               slim.col.label = TRUE,group.cex = 6,
                 group.label.rot = T,cex.row = 0,remove.key =F)
jpeg(paste0(path,"0-DoHeatmap-all.jpeg"), units="in", width=10, height=7,
     res=600)
print(g)
dev.off()

exp <- DoHeatmap_exp(object, genes.use = markers.to.plot,use.scaled = T,group.order = group.order)
write.csv(exp, paste0(path,"0-DoHeatmap-all~.csv"))
#---DoHeatmap vertical bar----
g1 <- MakeCorlorBar(object, markers,do.print = F, do.return = T)
jpeg(paste0(path,"0-DoHeatmap-all_legend.jpeg"),units="in", width=10, height=7,res=600)
print(g1)
dev.off()


##############################
#  DoHeatmap- epi ============
##############################
table(object@ident)
epithelial_ident <- c("9 Basal cells",
                      "10 Ciliated cells",
                      "8 Mucus-producing cells\nEpithelial cells",
                      "12 Distal secretory cells",
                      "2 Alveolar type II cells")
epithelial <- SubsetData(object, ident.use = epithelial_ident)
epithelial@ident = factor(epithelial@ident, levels = epithelial_ident)


# 1. 19 Lineage TF-list above (see item 2-c2)
TF = c("HMGA1", "TP63", "SOX15", "SOX7", "ELF3", "SPDEF", "HES4", "RARRES1",
       "XBP1", "SOX4", "FOXJ1", "RFX2", "TP73", "MYB", "SOX2", "HOPX", "NKX2-1", "ETV5", "FOXA2")
markers.to.plot <- FilterGenes(epithelial, TF)

#4. Heatmap for epithelial clusters only: CL9, CL10, CL8, CL12, CL2 (in this order, if possible) using all markers for epithelial clusters with adjusted p value <10(-50)
markers = read.csv(paste0("output/20190311/markers.Epi.csv"))
markers$order <- plyr::mapvalues(x = markers$cluster,
                                 from = c(9,10,8,12,2),
                                 to =   unique(markers$cluster))
markers = markers[order(markers$order),]
markers = markers[(markers$p_val_adj<0.1^50),]
dim(markers)
top <-  markers %>% 
        group_by(cluster) %>% 
        top_n(500, avg_logFC) %>% 
        as.data.frame
markers.to.plot <- top$gene
length(markers.to.plot)
g <- DoHeatmap(epithelial, genes.use = markers.to.plot,use.scaled = T, 
               #title= paste("19 Lineage TF in each cluster"),
               title= paste("Top DE genes (adjusted p value < 10^(-50)) in all epithelial cell types"),
               slim.col.label = TRUE,group.cex = 15,
               group.label.rot = T,cex.row = 0,remove.key =F)
jpeg(paste0(path,"4-DoHeatmap-epi.jpeg"), units="in", width=10, height=7,
     res=600)
print(g)
dev.off()

exp <- DoHeatmap_exp(epithelial, genes.use = markers.to.plot,use.scaled = T,group.order = epithelial_ident)
write.csv(exp, paste0(path,"4-DoHeatmap-epi.csv"))
#---DoHeatmap vertical bar----
g1 <- MakeCorlorBar(epithelial, markers,do.print = F, do.return = T)
jpeg(paste0(path,"4-DoHeatmap-epi_legend.jpeg"),units="in", width=10, height=7,res=600)
print(g1)
dev.off()
##############################
# heatmap.2
##############################
heatmap_markers <- read_excel("doc/Renat.markers.xlsx",sheet = "Dotplot_epithelial")
heatmap_markers <- FilterGenes(epithelial, heatmap_markers$Gene)
y = epithelial@scale.data[heatmap_markers,]

x = epithelial@meta.data[colnames(y),c("res.0.8","manual")]
z = merge(x,t(y), by = "row.names")
rownames(z) = z$Row.names
z = z[,-1]
colnames(z)[1:2] = c("tSNE_cluster","cell_type")
write.csv(t(z),file = paste0(path,"epithelial_exp.csv"))
## Column clustering (adjust here distance/linkage methods to what you need!)
hc <- hclust(as.dist(1-cor(as.matrix(y), method="spearman")), method="complete")
label_colors <- epithelial@meta.data[hc$labels,c("manual","manual.colors")] %>%
        apply(2,as.character)
head(label_colors)

jpeg(paste0(path,"/Heatmap2_epithelial_cells.jpeg"), units="in", width=10, height=7,res=600)
heatmap.2(as.matrix(y),
          Colv = as.dendrogram(hc), Rowv= TRUE,
          ColSideColors = label_colors[,2], trace ="none",labCol = FALSE,dendrogram = "column",
          key.xlab = "scale log nUMI",
          cexRow = 0.8,
          margins = c(2,5),
          #scale = "row",
          breaks = seq(-3,3,length.out = 101),
          col = bluered,
          main = "Hierarchical clustering of lung epithelial cells")
par(lend = 1)           # square line ends for the color legend
legend(0, 0.8,       # location of the legend on the heatmap plot
       legend = label_colors[!duplicated(label_colors[,1]),1], # category labels
       col = label_colors[!duplicated(label_colors[,1]),2],  # color key
       lty= 1,             # line style
       lwd = 10,            # line width
       cex=0.8
)
dev.off()

#####################
#Volcano plot - DEG analysis
####################
# FindPairMarkers
epithelial %<>% SetAllIdent("res.0.8")
ident.1 <- list(9,9,12,12,c(8,9,10,12))
ident.2 <- list(c(2,8,10,12),c(8,12,10),8,2,2)
subfolder <- paste0(path,"DEG/")
gde.pair <- FindPairMarkers(epithelial, ident.1 = ident.1, ident.2 = ident.2,
                            logfc.threshold = 0.25, min.cells.group =3,
                            return.thresh = 0.01, only.pos = FALSE, save.path = subfolder)
write.csv(gde.pair, paste0(subfolder,"pairwise_comparision.csv"))
head(gde.pair,10) %>% kable %>% kable_styling

titles <- c("Basal cells (CL9) vs. all other epithelial cells (CL8, CL12, CL10, CL2)",
            "Basal cells (CL9) vs. all other airway epithelial cells (CL8, CL12, CL10)",
            "Distal (CL12) vs. Secretory (CL8)",
            "Distal (CL12) vs. Alveolar (CL2)",
            "Airway epithelial (CL9, CL8, CL10, CL12) vs Alveolar (CL2)")
# Volcano plot=========
(clusters <- unique(gde.pair$cluster1.vs.cluster2))
for(i in 1:length(clusters)){
        df <- gde.pair[gde.pair$cluster1.vs.cluster2 %in% clusters[i],]
        g <- ggplot(df,aes(avg_logFC,-log10(p_val_adj))) + 
                geom_point() + 
                ggtitle(titles[i]) + 
                geom_text_repel(aes(label = gene), 
                                data=df[(df$avg_logFC > 1 | df$avg_logFC < -1)& 
                                         df$p_val_adj < head(sort(df$p_val_adj),40) %>% tail(1) ,]) +
                geom_point(color = ifelse((df$avg_logFC > 1  & df$p_val_adj < 0.05) , "red",
                                   ifelse((df$avg_logFC < -1 & df$p_val_adj < 0.05), "blue","gray")))
        
        jpeg(paste0(path,"Volcano_plot",clusters[i],".jpeg"), units="in", width=10, height=7,res=600)
        print(g)
        dev.off()
}

