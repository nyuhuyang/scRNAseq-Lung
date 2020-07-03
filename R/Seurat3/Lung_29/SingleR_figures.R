library(SingleR)
library(Seurat)
library(magrittr)
library(pheatmap)
library(kableExtra)
library(dplyr)
library(tidyr)
library(ggpubr)
library(openxlsx)
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat3_functions.R")
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/SingleR_functions.R")

path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
#====== 3.2 SingleR specifications ==========================================
# Step 1: Spearman coefficient
object = readRDS(file = "data/Lung_29_20200617.rds")
singler = readRDS(file="output/singlerT_Lung_20200630.rds")

# if singler didn't find all cell labels
length(singler$singler[[1]]$SingleR.single$labels) == ncol(object)

##############################
# create singleR data frame
###############################

if(length(singler$singler[[1]]$SingleR.single$labels) < ncol(object)){
        all.cell = colnames(object);length(all.cell)
        know.cell = rownames(singler$singler[[1]]$SingleR.single$labels);length(know.cell)
        unknown.cell = all.cell[!(all.cell %in% know.cell)]; length(unknown.cell)
        singlerDF = data.frame("singler1sub" = singler$singler[[1]]$SingleR.single$labels,
                               "sub_types" = singler$singler[[1]]$SingleR.single.main$labels,
                               row.names = rownames(singler$singler[[1]]$SingleR.single$labels),
                               stringsAsFactors = F)
        unknownDf = data.frame("singler1sub" = rep("unknown",length(unknown.cell)),
                               "sub_types" = rep("unknown",length(unknown.cell)),
                               row.names = unknown.cell,
                               stringsAsFactors = F)
        singlerDF = rbind(singlerDF,unknownDf)
}

if(length(singler$singler[[1]]$SingleR.single$labels) > ncol(object)){
}


singlerDF$main_types = gsub("\\ .*","", singlerDF$sub_types)
singlerDF$orig.ident = gsub("_.*","", rownames(singlerDF))
singlerDF = singlerDF[colnames(object),]
table(rownames(singlerDF) == colnames(object))

##############################
# cell type summary
##############################
# - Table: number of cells per cell types (per each sample and total)
opts_pairs <- data.frame(c("singler1sub" , "annotations3"),
                         c("main_types","orig.ident"),
                         c("sub_types" , "orig.ident"),
                         c("singler1sub" , "orig.ident"),
                         c("singler1sub" , "SCT_snn_res.2"),
                         c("singler1sub" , "SCT_snn_res.3.5"),
                         c("singler1sub" , "SCT_snn_res.4.8"),
                         c("singler1sub" , "SCT_snn_res.4.9"),
                       stringsAsFactors = F) %>% t
colnames(opts_pairs) = c("types","clus_sample")
df_list <- list()
for(i in 1:nrow(opts_pairs)){
          type = opts_pairs[,"types"][i]
          clus_sample = opts_pairs[,"clus_sample"][i]
          df <- table(singlerDF[,type], object@meta.data[,clus_sample]) %>% 
            as.data.frame()
          colnames(df) = c(type,clus_sample,"Freq")
          df %<>% spread(clus_sample,"Freq")
          rownames(df) = df[,type]
          df = df[order(df[,type]),]
          df_list[[i]] = df
          names(df_list)[i] = paste0(type, "_vs_",clus_sample)
}

write.xlsx(df_list, file = paste0(path,"cell.type_summary.xlsx"),
           colNames = TRUE, borders = "surrounding",colWidths = c(NA, "auto", "auto"))
##############################
# process color scheme
##############################
#singler_colors <- readxl::read_excel("doc/singler.colors.xlsx")
#singler_colors1 = as.vector(singler_colors$singler.color1[!is.na(singler_colors$singler.color1)])
#singler_colors1[duplicated(singler_colors1)]
#singler_colors2 = as.vector(singler_colors$singler.color2[!is.na(singler_colors$singler.color2)])
object %<>% AddMetaData(metadata = singlerDF)
object <- AddMetaColor(object = object, label= "main_types", colors = Singler.colors)
object <- AddMetaColor(object = object, label= "sub_types", colors = Singler.colors)

lapply(c(FALSE, TRUE), function(label)
    UMAPPlot.1(object = object, label = label, label.repel = T,
               group.by = "main_types",
               no.legend = T,
               pt.size = 0.1,label.size = 4, do.print = T,do.return = F,
           title = "main_types in HIV"))

saveRDS(object,file="data/HIVseurat.rds")

##############################
# draw tsne plot
##############################
object <- subset(object,idents = c("HSC/progenitors","Nonhematopoietic cells"), invert = TRUE)
Idents(object) = "cell.types"
object %<>% sortIdent()
table(Idents(object))


cell_Freq <- table(Idents(object)) %>% as.data.frame
cell_Freq = cell_Freq[order(cell_Freq$Var1),]
cell_Freq$col = ExtractMetaColor(object)
cell_Freq = cell_Freq[order(cell_Freq$Freq,decreasing = T),]
cell_Freq$Var1 %<>% factor(levels = as.character(cell_Freq$Var1))
colnames(cell_Freq)[1:2] = c("Cell_Type", "Cell_Number")

jpeg(paste0(path,"cell_type_numbers.jpeg"), units="in", width=10, height=7,res=600)
ggbarplot(cell_Freq, "Cell_Type", "Cell_Number",
          fill = "Cell_Type", color = "Cell_Type",xlab = "",
          palette = cell_Freq$col,x.text.angle = 90,
          title = "Numbers of major cell types in total 43 samples")+NoLegend()+
  theme(plot.title = element_text(hjust = 0.5,size=15))
dev.off()
