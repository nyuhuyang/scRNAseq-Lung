library(SingleR)
library(Seurat)
library(reshape2)
library(pheatmap)
library(kableExtra)
library(dplyr)
library(tidyr)
library(magrittr)
source("../R/Seurat_functions.R")
source("../R/SingleR_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
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
object <- SetAllIdent(object = object, id = "manual")
p3 <- TSNEPlot.1(object = object, do.label = F, group.by = "ident",
                 do.return = TRUE, no.legend = F,
                 colors.use = ExtractMetaColor(object),
                 pt.size = 1,label.size = 3,force = 2)+
  ggtitle("Supervised cell type labeling by Blueprint + Encode")+
  theme(text = element_text(size=10),							
        plot.title = element_text(hjust = 0.5,size = 18, face = "bold")) 

jpeg(paste0(path,"PlotTsne-combine.jpeg"), units="in", width=10, height=7,
     res=600)
print(p3)
dev.off()

save(object,file="data/Lung_Harmony_3_20190109.Rda")
##############################
# subset Seurat
###############################
table(object@meta.data$orig.ident)
table(object@ident)
df_samples <- readxl::read_excel("doc/20190119_scRNAseq_info.xlsx")
colnames(df_samples) <- tolower(colnames(df_samples))
sample_n = which(df_samples$tests %in% paste0("test",0:1))
df <- as.data.frame(df_samples[sample_n,])
samples <- unique(df$sample)

object %<>% SetAllIdent(id = "orig.ident")

lapply(samples,function(sample) {
        SubsetData(object, ident.use = sample) %>%
                SetAllIdent(id = "manual") %>%
                TSNEPlot.1(no.legend = F,do.label =F,label.size=3,size=20,
                           colors.use = ExtractMetaColor(.),
                           return.plots =T, label.repel = T,force=2,
                           do.print = T, do.return = F)
        })

##############################
# subset Seurat
###############################
object %<>% SetAllIdent(id = "res.0.8")
Epi <- SubsetData(object, ident.use = c(2,8,9,10,12))
TSNEPlot(Epi,do.label = T)
remove <- FeaturePlot(Epi,features.plot = "KRT19",do.identify = T)
keep  <- Epi@cell.names[!(Epi@cell.names %in% remove)]
Epi <- SubsetData(object,cells.use = keep)
Epi %<>% SetAllIdent("manual")
g <- TSNEPlot.1(Epi,do.label =F,no.legend = F,label.repel = F,pt.size = 2,
           colors.use = ExtractMetaColor(Epi))+
        #ggtitle("Epithelial cells labeling by Blueprint + Encode")+
        theme(text = element_text(size=10),							
              plot.title = element_text(hjust = 0.5,size = 18, face = "bold")) 
jpeg(paste0(path,"PlotTsne-Epi.jpeg"), units="in", width=10, height=7,res=600)
print(g)
dev.off()
