########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################
# conda activate r4.0
invisible(lapply(c("Seurat","dplyr","cowplot","kableExtra",
                   "magrittr"), function(x) {
                           suppressPackageStartupMessages(library(x,character.only = T))
                   }))
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)

set.seed(101)
#======1.2 load  Seurat =========================
# load files
object = readRDS(file = "data/Lung_30_20200710.rds") 
DefaultAssay(object) = "SCT"
# Extend Data
Idents(object) = "cell_types"

cell_Freq <- table(Idents(object)) %>% as.data.frame
cell_Freq$Percent <- prop.table(cell_Freq$Freq) %>% round(digits = 3) %>% scales::percent()

cell_Freq = cell_Freq[order(cell_Freq$Var1),]
cell_Freq$col = ExtractMetaColor(object)
cell_Freq = cell_Freq[order(cell_Freq$Freq,decreasing = T),]
cell_Freq$Var1 %<>% factor(levels = as.character(cell_Freq$Var1))
colnames(cell_Freq)[1:2] = c("Cell_Type", "Cell_Number")

jpeg(paste0(path,"cell_type_numbers.jpeg"), units="in", width=10, height=6,res=600)
ggbarplot(cell_Freq, "Cell_Type", "Cell_Number",
          fill = "Cell_Type", color = "black",xlab = "",
          palette = cell_Freq$col,x.text.angle = 90,
          ylab = "Cell Number",
          #label = "Percent",
          sort.val = "desc",
          width = 1, size = 0.5,
          lab.pos = c("out", "in"),
          title = "Numbers of cell types in total 30 samples")+NoLegend()+
        theme(plot.title = element_text(hjust = 0.5,size=15),
              axis.text.x = element_text(vjust = 0.5))
dev.off()

#=======
Seurat_list <- SplitObject(object, split.by = "orig.ident")

QC_list <- as.data.frame(df_samples)
QC_list["cell.number"] <- sapply(Seurat_list, function(x) length(colnames(x)))
QC_list["mean.nUMI"] <- sapply(Seurat_list, function(x) mean((x$nCount_SCT)))
QC_list["mean.nGene"] <- sapply(Seurat_list, function(x) mean((x$nFeature_SCT)))
QC_list["mean.percent.mt"] <- sapply(Seurat_list, function(x) mean((x$percent.mt)))

write.csv(QC_list,paste0(path,"QC_list.csv"))
#QC.list %>% kable() %>% kable_styling()

remove(Seurat_list);GC()

metrics_summary <- read.csv("output/20210407/metrics_summary.csv", stringsAsFactors = F)
metrics_summary$Mean.Reads.per.Cell %<>% gsub(",","",.) %>% as.integer()
metrics_summary$submitter = "Rustam"
jpeg(paste0(path,"mean.Reads.per.Cell.jpeg"), units="in", width=5, height=5,res=600)
ggviolin(metrics_summary, y= "Mean.Reads.per.Cell",
         title = "Mean reads per cell in each scRNA-seq",
         xlab = "",ylab = "Mean Reads per Cell",
         add = c("jitter","mean_sd"),
         draw_quantiles = 0.5,
         yscale = "log10"
)+
        #scale_y_continuous(expand = c(0, 0), limits = c(0,350000))+
        theme(plot.title = element_text(hjust = 0.5,size=15),
              axis.title.x=element_blank(),
              axis.text.x=element_blank())
dev.off()

jpeg(paste0(path,"cell.number.jpeg"), units="in", width=5, height=5,res=600)
ggviolin(QC_list, y= "cell.number",
         title = "Cell number in each scRNA-seq",
         xlab = "",ylab = "Cell Number",
         add = c("jitter","mean_sd"),
         #ylim = c(0, max(QC_list$cell.number)+1000),
         draw_quantiles = 0.5,
         yscale = "log10"
)+
        #scale_y_continuous(expand = c(0, 0), limits = c(0,6000))+
        theme(plot.title = element_text(hjust = 0.5,size=15),
              axis.title.x=element_blank(),
              axis.text.x=element_blank())
dev.off()

jpeg(paste0(path,"UMI.per.Cell.jpeg"), units="in", width=5, height=5,res=600)
ggviolin(QC_list, y= "mean.nUMI",
         title = "Mean transcripts per cell in each scRNA-seq",
         xlab = "",ylab = "Mean UMI per Cell",
         add = c("jitter","mean_sd"),
         draw_quantiles = 0.5,
         yscale = "log10"
)+
        #scale_y_continuous(expand = c(0, 0), limits = c(0,12000))+
        theme(plot.title = element_text(hjust = 0.5,size=13),
              axis.title.x=element_blank(),
              axis.text.x=element_blank())
dev.off()

jpeg(paste0(path,"nGene.per.Cell.jpeg"), units="in", width=5, height=5,res=600)
ggviolin(QC_list, y= "mean.nGene",
         title = "Mean genes per cell in each scRNA-seq",
         xlab = "",ylab = "Mean genes per Cell",
         add = c("jitter","mean_sd"),
         draw_quantiles = 0.5,
         yscale = "log10"
)+
        #scale_y_continuous(expand = c(0, 0), limits = c(0,4000))+
        theme(plot.title = element_text(hjust = 0.5,size=13),
              axis.title.x=element_blank(),
              axis.text.x=element_blank())
dev.off()

jpeg(paste0(path,"percent.mt.jpeg"), units="in", width=5, height=5,res=600)
ggviolin(QC_list, y= "mean.percent.mt",
         title = "Mean mitochondrial gene % in each scRNA-seq",
         xlab = "",ylab = "Mean mitochondrial gene percentages",
         add = c("jitter","mean_sd"),
         draw_quantiles = 0.5,
         #yscale = "log10"
)+
        scale_y_continuous(expand = c(0, 0), limits = c(0,15))+
        theme(plot.title = element_text(hjust = 0.5,size=13),
              axis.title.x=element_blank(),
              axis.text.x=element_blank())
dev.off()
