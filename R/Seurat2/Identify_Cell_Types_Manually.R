library(Seurat)
library(dplyr)
library(tidyr)
library(kableExtra)
library(gplots)
source("../R/Seurat_functions.R")
source("../R/SingleR_functions.R")
path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
#====== 3.1 identify phenotype for each cluster  ==========================================
(load(file="data/MCL_Harmony_20_20181231.Rda"))

df_markers <- readxl::read_excel("../SingleR/output/bio-rad-markers.xlsx")
df_markers <- readxl::read_excel("./doc/Renat.markers.xlsx")
colnames(df_markers) = gsub(" ","_",colnames(df_markers))
colnames(df_markers) = gsub(":|\\/","_",colnames(df_markers))
colnames(df_markers) = gsub("\\+","",colnames(df_markers))
col_index <- grep("Alias",colnames(df_markers))
if(!identical(col_index,integer(0))) {
        markers = df_markers[,-grep("Alias",colnames(df_markers))]
} else markers = df_markers
marker.list <- df2list(markers)

marker.list %<>% lapply(function(x) HumanGenes(object,x))  %>%
        lapply(function(x) x[1:12]) %>% lapply(function(x) x[!is.na(x)])

marker.list %>% list2df %>% t %>% kable() %>% kable_styling()

for(i in 1:length(marker.list)){
        p <- lapply(marker.list[[i]], function(marker) {
                SingleFeaturePlot.1(object = object, feature = marker,pt.size = 0.5,
                                    gradient.use = c("lightblue", "blue3"),threshold=0.1)+
                        ggtitle(paste0(marker,Alias(df = df_markers,gene = marker)))+ #
                        theme(plot.title = element_text(hjust = 0.5,size = 12, face = "bold"))
                })

        jpeg(paste0(path,names(marker.list)[i],".jpeg"),
             units="in", width=10, height=7,res=600)
        print(do.call(plot_grid, p)+ ggtitle(paste(names(marker.list)[i],"markers"))+
                      theme(plot.title = element_text(hjust = 0.5,size = 20, face = "bold")))
        dev.off()
        print(paste0(i,":",length(marker.list)))
}
#=======3.2 dotplot================
markers.to.plot <- unlist(marker.list)
marker_order <- c(3,13,11,16,5,6,15,10,14,0,12,4,7,8,2,1,9)
object@ident <- factor(x = object@ident, levels = marker_order) # Relevel object@ident
jpeg(paste0(path,"S2B_dotplot.jpeg"), units="in", width=10, height=7,
     res=600)
DotPlot(object, genes.plot = rev(markers.to.plot),
        cols.use = c("blue","red"), x.lab.rot = T, plot.legend = T,
        dot.scale = 8, do.return = T)
dev.off()

#====== 3.4 Rename ident =========================
(load(file="data/MCL_Harmony_20_20181231.Rda"))
object <- SetAllIdent(object, id = "res.0.6")
table(object@ident)

idents <- as.data.frame(table(object@ident))
old.ident.ids <- idents$Var1
new.cluster.ids <- c("Stromal_fibroblasts",
                     "Ciliated_cells",
                     "Mucus-producing_cells/Secretory_club_cells",
                     "Macrophages/Monocytes",
                     "Distal_secretory_cells",
                     "Vascular_endothelial_cells",
                     "Vascular_endothelial_cells",
                     "Secretory_club_cells",
                     "Squamous_cells",
                     "Epithelial_cells",
                     "Smooth_muscle_cells",
                     "T/NK_cells",
                     "Stromal_fibroblasts?",
                     "B_cells",
                     "Smooth_muscle_cells",
                     "Vascular_endothelial_cells",
                     "Lymphatic_endothelial_cells")

object@ident <- plyr::mapvalues(x = object@ident,
                              from = old.ident.ids,
                              to = new.cluster.ids)
TSNEPlot.1(object = object,do.label = T, group.by = "ident", 
           do.return = TRUE, no.legend = T, colors.use = singler.colors,
           pt.size = 1,label.size = 4,label.repel = T)+
        ggtitle("Discovery cell types by marker genes")+
        theme(text = element_text(size=20),							
              plot.title = element_text(hjust = 0.5,size=18, face = "bold"))
object <- StashIdent(object = object, save.name = "manual")
object@meta.data$manual %>% table() %>% kable() %>% kable_styling()
##############################
# process color scheme
##############################
singler_colors <- readxl::read_excel("./doc/singler.colors.xlsx")
singler_colors2 = as.vector(singler_colors$singler.color2[!is.na(singler_colors$singler.color2)])
singler_colors2[duplicated(singler_colors2)]
length(singler_colors2)
apply(object@meta.data[,c("singler2sub","manual")],
      2,function(x) length(unique(x)))
object@meta.data[,c("manual")] %>% table() %>% kable() %>% kable_styling()
object <- AddMetaColor(object = object, label= "manual", colors = singler_colors2[1:14])
object <- SetAllIdent(object = object, id = "manual")
TSNEPlot.1(object, colors.use = ExtractMetaColor(object),no.legend = F)
##############################
# draw tsne plot
##############################
p4 <- TSNEPlot.1(object = object, do.label = F, group.by = "ident",
                 do.return = TRUE, no.legend = F,
                 colors.use = ExtractMetaColor(object),
                 pt.size = 1,label.size = 3,force = 2)+
        ggtitle("Identify cell types by marker genes")+
        theme(text = element_text(size=10),							
              plot.title = element_text(hjust = 0.5,size = 18, face = "bold")) 

jpeg(paste0(path,"PlotTsne_manual.jpeg"), units="in", width=10, height=7,
     res=600)
print(p4)
dev.off()

save(object,file="data/Lung_Harmony_3_20190109.Rda")

g <- SplitTSNEPlot(object,group.by = "ident",split.by = "orig.ident",
                   #select.plots = c(1,3,2),#c(6:8,1:5)
                   no.legend = T,do.label =T,label.size=3,size=20,
                   return.plots =T, label.repel = T,force=2)
jpeg(paste0(path,"SplitTSNEPlot3.jpeg"), units="in", width=10, height=7,
     res=600)
g[[3]]
dev.off()
table(object.all@meta.data$manual, object@meta.data$orig.ident) %>%
        prop.table(margin = 2) %>% kable %>%
        kable_styling()

#====== 2.2 focus on MCL heterogeneity, and T cell subsets and NK cells in the tSNE plots =====
df_markers <- readxl::read_excel("doc/MCL-markers.xlsx")
(markers <- df_markers[,1] %>% as.matrix %>% as.character %>% HumanGenes(object,marker.genes = .))

table(object@meta.data$orig.ident)
table(object@ident)
object@meta.data$orig.ident = gsub("BH|DJ|MD|NZ","Normal",object@meta.data$orig.ident)

df_samples <- readxl::read_excel("doc/181227_scRNAseq_info.xlsx")
colnames(df_samples) <- tolower(colnames(df_samples))
tests <- paste0("test",c(3))
for(test in tests){
        sample_n = which(df_samples$tests %in% c("control",test))
        samples <- unique(df_samples$sample[sample_n])
        
        cell.use <- rownames(object@meta.data)[object@meta.data$orig.ident %in% 
                                                       c("Normal",samples)]
        subset.object <- SubsetData(object, cells.use = cell.use)
        subset.object@meta.data$orig.ident %>% unique %>% sort %>% print
        
        SplitSingleFeaturePlot(subset.object, 
                               select.plots = c(1:4,6,5),#c(6:8,1:5)
                               alias = df_markers, 
                               group.by = "ident",split.by = "orig.ident",
                               no.legend = T,label.size=3,do.print =T,
                               markers = markers, threshold = NULL)
}

object@meta.data$orig.ident =gsub("Pt-11-LN-C14","Pt-11-C14",object@meta.data$orig.ident)
