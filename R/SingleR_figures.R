library(SingleR)
library(Seurat)
library(reshape2)
library(pheatmap)
library(kableExtra)
library(dplyr)
library(tidyr)
source("../R/Seurat_functions.R")
source("../R/SingleR_functions.R")
path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)

#====== 3.2 SingleR specifications ==========================================
# Step 1: Spearman coefficient
(load(file="data/Lung_Harmony_3_20190109.Rda"))
(load(file="./output/singler_Lung_3F_20190109.Rda"))

# if singler didn't find all cell labels
if(length(singler$singler[[1]]$SingleR.single$labels) != ncol(object@data)){
        all.cell = object@cell.names;length(all.cell)
        know.cell = rownames(singler$singler[[1]]$SingleR.single$labels);length(know.cell)
        object = SubsetData(object, cells.use = know.cell)
}
object
##############################
# add singleR label to Seurat
###############################
singlerDF = data.frame("singler1sub"=singler$singler[[1]]$SingleR.single$labels,
                       "singler1main"=singler$singler[[1]]$SingleR.single.main$labels,
                       "singler2sub"=singler$singler[[2]]$SingleR.single$labels,
                       "singler2main"=singler$singler[[2]]$SingleR.single.main$labels,
                       "kang" = FineTune(singler$other, main.type = FALSE),
                       row.names = rownames(singler$singler[[1]]$SingleR.single$labels))

table(rownames(singlerDF) %in% object@cell.names)

apply(singlerDF,2,function(x) length(unique(x)))
object <- AddMetaData(object = object,
                   metadata = singlerDF)
object <- SetAllIdent(object = object, id = "singler2sub")

##############################
# check the spearman correlation
###############################
#Or by all cell types (showing the top 50 cell types):
jpeg(paste0(path,"/DrawHeatmap_sub1.jpeg"), units="in", width=10, height=7,
     res=600)
print(SingleR.DrawHeatmap(singler$singler[[1]]$SingleR.single, top.n = 50,normalize = F))
dev.off()
jpeg(paste0(path,"/DrawHeatmap_sub1_N.jpeg"), units="in", width=10, height=7,
     res=600)
print(SingleR.DrawHeatmap(singler$singler[[1]]$SingleR.single,top.n = 50,normalize = T))
dev.off()

#Finally, we can also view the labeling as a table compared to the original identities:

kable(table(singler$singler[[2]]$SingleR.single$labels, singler$meta.data$orig.ident)) %>%
        kable_styling()
singler$meta.data$orig.ident %>% table() %>% kable() %>% kable_styling()
object@meta.data$singler2sub %>% table() %>% kable() %>% kable_styling()

##############################
# process color scheme
##############################
singler_colors <- readxl::read_excel("./doc/singler.colors.xlsx")
singler_colors1 = as.vector(singler_colors$singler.color1[!is.na(singler_colors$singler.color1)])
singler_colors1[duplicated(singler_colors1)]
length(singler_colors1)
apply(object@meta.data[,c("singler2sub","singler2main")],
      2,function(x) length(unique(x)))
object@meta.data[,c("singler2sub")] %>% table() %>% kable() %>% kable_styling()
object <- AddMetaColor(object = object, label= "singler2sub", colors = singler_colors1[1:33])
object <- SetAllIdent(object = object, id = "singler2sub")
TSNEPlot.1(object, colors.use = ExtractMetaColor(object),no.legend = F)
##############################
# draw tsne plot
##############################
object <- SetAllIdent(object = object, id = "singler2sub")
p3 <- TSNEPlot.1(object = object, do.label = T, group.by = "ident",
                 do.return = TRUE, no.legend = T,
                 colors.use = singler.colors,#ExtractMetaColor(object),
                 pt.size = 1,label.size = 3,force = 2)+
  ggtitle("Supervised cell type labeling by Blueprint + Encode")+
  theme(text = element_text(size=10),							
        plot.title = element_text(hjust = 0.5,size = 18, face = "bold")) 

jpeg(paste0(path,"PlotTsne_sub2.jpeg"), units="in", width=10, height=7,
     res=600)
print(p3)
dev.off()

save(object,file="data/Lung_Harmony_3_20190109.Rda")
##############################
# subset Seurat
###############################
table(object@meta.data$orig.ident)
table(object@ident)
object@meta.data$orig.ident = gsub("BH|DJ|MD|NZ","Normal",object@meta.data$orig.ident)

df_samples <- readxl::read_excel("doc/181227_scRNAseq_info.xlsx")
colnames(df_samples) <- tolower(colnames(df_samples))
object <- SetAllIdent(object, id="singler1sub")
tests <- paste0("test",c(3))
for(test in tests){
        sample_n = which(df_samples$tests %in% c("control",test))
        samples <- unique(df_samples$sample[sample_n])

        cell.use <- rownames(object@meta.data)[object@meta.data$orig.ident %in% 
                                                       c("Normal",samples)]
        subset.object <- SubsetData(object, cells.use = cell.use)
        subset.object@meta.data$orig.ident %>% unique %>% sort %>% print
        g <- SplitTSNEPlot(subset.object,group.by = "ident",split.by = "orig.ident",
                           select.plots = c(1,3,2),#c(6:8,1:5)
                           no.legend = T,do.label =F,label.size=3,size=20,
                           return.plots =T, label.repel = T,force=2)
        jpeg(paste0(path,test,"_TSNEPlot.jpeg"), units="in", width=10, height=7,
             res=600)
        print(do.call(plot_grid, c(g, ncol = 2)))
        dev.off()
}

###################################
# creat seurat object
(load(file = 'data/ref_MCL_blue_encode.RData'))
MCL_blue_encode <- CreateSeuratObject(raw.data = ref$data) %>%
  NormalizeData %>% ScaleData

MCL_blue_encode@data = ref$data
ident = ref$main_types
ident[!grepl("B_cells|MCL",ident)] = "others"
names(ident) = MCL_blue_encode@cell.names
MCL_blue_encode@ident = factor(ident, levels = c("others","B_cells","MCL"))
MCL_B_markers <- FindAllMarkers.UMI(MCL_blue_encode, logfc.threshold = 0.1, only.pos = T,
                                   test.use = "MAST")
write.csv(MCL_B_markers, paste0(path,"MCL_B_markers.csv"))
g <- DoHeatmap.1(MCL_blue_encode, MCL_B_markers,Top_n = 50,
                 ident.use = paste("B cells vs. MCL vs. other cell types"),
                 group.label.rot = T,cex.row = 4,remove.key =T)
jpeg(paste0(path,"heatmap_B_MD_others.jpeg"), units="in", width=10, height=7,
     res=600)
print(g)
dev.off()
