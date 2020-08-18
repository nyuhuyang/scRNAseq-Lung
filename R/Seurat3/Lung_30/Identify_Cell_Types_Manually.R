library(Seurat)
library(dplyr)
library(tidyr)
library(kableExtra)
library(magrittr)
library(readxl)
library(cowplot)
library(stringr)
library(openxlsx)
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)
# load data
object = readRDS(file = "data/Lung_30_20200702.rds")
DefaultAssay(object) = "SCT"
#======== rename ident =================
object %<>% FindClusters(resolution = 4.1)
object %<>% FindClusters(resolution = 4.9)
object %<>% FindClusters(resolution = 5)

df_annotation <- readxl::read_excel("doc/Annotations/20200704_Annotations.xlsx",
                                    sheet = "Sheet1")
resolutions = unique(df_annotation$Resolution) %>% .[complete.cases(.)]
df_res_list <- list()
for(i in seq_along(resolutions)){
        res = resolutions[i]
        temp <- df_annotation[df_annotation$Resolution %in% res,]
        Cluster <- stringr::str_split(string = temp$Cluster, 
                                        pattern = "\\+")
        df_res <- data.frame("Cluster" = as.integer(unlist(Cluster)), stringsAsFactors = F)
        for(m in 1:nrow(df_res)) {
                c = as.character(df_res$Cluster[m])
                c = paste0("\\b",c,"\\b")
                df_res$Celltype[m] = temp[grep(c, Cluster),"Cell type"]
        }
        df_res$Resolution = as.character(res)
        df_res_list[[i]] = df_res
}
df_res = bind_rows(df_res_list);rm(df_res_list);GC()
df_res = df_res[complete.cases(df_res$Cluster),]
#===========================================
meta.data = cbind(object@meta.data,  object@reductions$umap@cell.embeddings)

meta.data$annotations3 = "Un"
for(m in 1:nrow(df_res)){
        resolution = paste0("SCT_snn_res.",df_res$Resolution[m])
        meta.data[meta.data[,resolution] %in% df_res$Cluster[m],
                  "annotations3"] = df_res$Celltype[m]
        Progress(m, nrow(df_res))
}


# label by expression
exp = FetchData(object, vars = c("SFTPC","SCGB1A1","SCGB3A2","SFTPB","MUC5AC",
                                 "KRT5","KRT14","MKI67","UBE2C","SERPINB4","LY6D",
                                 "SERPINB3","KRT15","AGER","7SK.2"))
C11 <- meta.data$SCT_snn_res.2 == 11
meta.data[C11 & exp$SFTPC < 2 & exp$SCGB1A1 >2 & (exp$SCGB3A2 > 1 | exp$SFTPB > 1),
          "annotations3"] = "S-d+"

C23 <- meta.data$SCT_snn_res.2 == 23
meta.data[C23 & (exp$SCGB3A2 > 1 | exp$SFTPB > 1), "annotations3"] = "S-d++"
meta.data[C23 & exp$MUC5AC > 2 , "annotations3"] = "S-Muc+"

C40 <- meta.data$SCT_snn_res.4 == 40
meta.data[C40 & (exp$KRT5 > 0 | exp$KRT14 > 0), "annotations3"] = "MEC"

C60 <- meta.data$SCT_snn_res.2 == 60
BC_p <- C60 & (exp$MKI67 > 0 | exp$UBE2C > 0) & exp$KRT5 > 0
IC2 <- C60 & !BC_p & (exp$SERPINB4 > 2 | exp$LY6D > 2)
IC1 <- C60 & !BC_p & !IC2 & exp$SERPINB3 > 1
BCplus <- C60 & !BC_p & !IC2 & !IC1 & (exp$KRT5 > 1 | exp$KRT15 >1)

meta.data[BC_p, "annotations3"] = "BC-p"
meta.data[IC2, "annotations3"] = "IC2+"
meta.data[IC1, "annotations3"] = "IC1+"
meta.data[BCplus, "annotations3"] = "BC+"

meta.data[meta.data$annotations3 %in% "S" & (exp$SCGB3A2>1 | exp$SFTPB>1),
          "annotations3"] = "S-d+++"

C88 <- meta.data$SCT_snn_res.5 == 88
table(meta.data[C88, "annotations3"])
AT2 <- C88 & meta.data$annotations3 == "AT2"
meta.data[AT2 & exp$SFTPC == 0 & (exp$SCGB3A2>0 | exp$SFTPB>0),
          "annotations3"] = "S-d1+"
meta.data[AT2 & exp$SFTPC < 2 & exp$SCGB1A1>2, "annotations3"] = "S-d2+"
meta.data$annotations3 %<>% gsub("S-d.*","S-d",.)
meta.data$annotations3 %<>% gsub("S-Muc.*","S",.)
meta.data$annotations3 %<>% gsub("BC\\+","BC",.)
meta.data[meta.data$annotations3 == "AT2" & exp$AGER>2, "annotations3"] = "AT2-1"

# label by corrdinates ============
rectangle_name <- function(df, x_left = -Inf, x_right = Inf, y_bottom = -Inf, y_top = Inf,
                           Markers = "", old.name, new.name){
        tmp = df[,"UMAP_1"] > x_left & df[,"UMAP_1"] < x_right & 
                df[,"UMAP_2"] > y_bottom & df[,"UMAP_2"] < y_top & 
                df$annotations3 %in% old.name
        if(Markers != "") tmp = tmp & eval(parse(text = Markers))
        df[rownames(df)[tmp], "annotations3"] = new.name
        return(df)
}
df_annotation1 <- readxl::read_excel("doc/Annotations/20200709_Labeling optimization.xlsx",
                                    sheet = "Sheet1")
df_annotation1$x_left[is.na(df_annotation1$x_left)] <- -Inf
df_annotation1$x_right[is.na(df_annotation1$x_right)] <- Inf
df_annotation1$y_bottom[is.na(df_annotation1$y_bottom)] <- -Inf
df_annotation1$y_top[is.na(df_annotation1$y_top)] <- Inf
df_annotation1$Markers[is.na(df_annotation1$Markers)] = ""

for(m in 1:nrow(df_annotation1)){
        a = df_annotation1[m,]
        meta.data %<>% rectangle_name(x_left = a$x_left,
                                      x_right = a$x_right,
                                      y_bottom = a$y_bottom,
                                      y_top = a$y_top,
                                      Markers = a$Markers,
                                      old.name = a$old.name,
                                      new.name = a$new.name)
}

# label by barcodes 20200713 ============
df_annotation2 <- readxl::read_excel("doc/Annotations/20200713_Cell IDs for re-labeling.xlsx",
                                     sheet = "Sheet1")
(dup_cell <- as.vector(df_annotation2$`Cell ID`[duplicated(df_annotation2$`Cell ID`)]))
#dup_df = df_annotation2[df_annotation2$`Cell ID` %in% dup_cell,] %>% as.data.frame()
#dup_df = dup_df[order(dup_df$`Cell ID`),]
#write.xlsx(dup_df, file = paste0(path,"duplicated_cells.xlsx"),
#           colNames = TRUE, borders = "surrounding",colWidths = c(NA, "auto", "auto"))
df_annotation2 %<>% as.data.frame()
meta.data[df_annotation2$`Cell ID`,"annotations3"] = df_annotation2$`Change to`
meta.data$annotations3 %<>% gsub("M1-M2","M1-2",.)
# label by barcodes 20200719 ============
df_annotation3 <- readxl::read_excel("doc/Annotations/20200719_BC-IC-S-re-labelingRS.xlsx",
                                     sheet = "Sheet1")
(dup_cell <- as.vector(df_annotation3$`Cells`[duplicated(df_annotation3$`Cells`)]))
df_annotation3 %<>% as.data.frame()
meta.data[df_annotation3$`Cells`,"annotations3"] = df_annotation3$`New Label`
# label by barcodes 20200719 ============
df_annotation4 <- readxl::read_excel("doc/Annotations/20200726_Re-labeling after BC-IC-S based on UMAP clustering data.xlsx",
                                     sheet = "Sheet1")
(dup_cell <- as.vector(df_annotation4$`Cells`[duplicated(df_annotation4$`Cells`)]))
df_annotation4 %<>% as.data.frame()
meta.data[df_annotation4$`Cells`,"annotations3"] = df_annotation4$`New Label`
# label by barcodes 20200802 ============
df_annotation5 <- readxl::read_excel("doc/Annotations/20200802_Re-labeling - for RS.xlsx",
                                     sheet = "Sheet1")
df_annotation5 = df_annotation5[!duplicated(df_annotation5$`Cells`),]
(dup_cell <- as.vector(df_annotation5$`Cells`[duplicated(df_annotation5$`Cells`)]))
df_annotation5 %<>% as.data.frame()
meta.data[df_annotation5$`Cells`,"annotations3"] = df_annotation5$`Change to`

# T7 cluster contains some ciliated (C) cells and neutrophils:
meta.data[meta.data$annotations3 == "T7","annotations3"] = "Un"
df_annotation6 <- readxl::read_excel("doc/Annotations/T-int revision.xlsx",
                                     sheet = "Sheet1")
meta.data[df_annotation6$`Cells`,"annotations3"] = df_annotation6$`Change to`

meta.data[exp$`7SK.2` >2 & 
          meta.data$SCT_snn_res.3 == 11 & 
          grepl("^T",meta.data$annotations3),"annotations3"] = "T7"
meta.data = cbind(meta.data,  exp$`7SK.2`)
meta.data[rownames(meta.data) %in% df_annotation6$`Cells` &
                  meta.data[,"exp$`7SK.2`"] >2, "annotations3"] = "T7"
# M0-M-M2 relabeling
df_annotation7 <- readxl::read_excel("doc/Annotations/20200809_M0-M1-M2 revision.xlsx",
                                     sheet = "Sheet1")
(dup_cell <- as.vector(df_annotation7$`Cells`[duplicated(df_annotation7$`Cells`)]))
meta.data[df_annotation7$`Cells`,"annotations3"] = df_annotation7$Final

# immune relabeling
df_annotation8 <- readxl::read_excel("doc/Annotations/20200814_T cell label revision.xlsx",
                                     sheet = "Sheet1")
(dup_cell <- as.vector(df_annotation8$`Cells`[duplicated(df_annotation8$`Cells`)]))
df_annotation8$Final= gsub("T-7","T7", df_annotation8$Final)
meta.data[df_annotation8$`Cells`,"annotations3"] = df_annotation8$Final


#===========================
# prepare UMAP
#===========================

object[["annotations3"]] = meta.data$annotations3

Idents(object)= "annotations3"

UMAPPlot.1(object,group.by = "annotations3", cols = Singler.colors, 
           label = F, label.repel = F, no.legend = T,do.return = F, do.print = T)
UMAPPlot.1(object,group.by = "annotations3", cols = Singler.colors, 
           label = T, label.repel = T, no.legend = T,do.return = F, do.print = T)
samples[!(samples %in% Idents(object))]
Idents(object)[!(Idents(object) %in% samples)]

PrepareShiny(object, samples, Rshiny_path, split.by = "annotations3", 
             reduction = "umap",verbose = T)

saveRDS(object, file = paste0("data/Lung_30_20200702.rds"))



# - Table: number of cells per cell types (per each sample and total)
df <- table(sub_object$annotations3, sub_object$orig.ident) %>% 
        as.data.frame()
colnames(df) = c("cell.types","samples","Freq")
df %<>% tidyr::spread("samples","Freq")
rownames(df) = df$cell.types
cols = c("UNC-48-P","UNC-55-P","UNC-66-P","UNC-69-P","UNC-71-P","UNC-48-D",
         "UNC-55-D","UNC-66-D","UNC-69-D","UNC-71-D","UNC-44-D","UNC-54-D",
         "UNC-57-D","UNC-67-D","UNC-70-D","CU-12-D","CU-12-D-R","VU-27-D",
         "UNC-48-T","UNC-55-T","UNC-66-T","UNC-69-T","UNC-71-T","UNC-44-T",
         "CU-12-T","UNC-51-D","UNC-52-D","UNC-61-D","VU19-D","VU-37-D")
rows = c("AT1","AT2","AT2-1","AT2-p","BC","BC-p","BC-S","IC1","IC2","IC-S","H","p-C",
         "C1","C2","C3","S","S-d","Ion","NEC","SMG-Muc","SMG-Ser","MEC","Cr",
         "F1","F2","F3","F4","Gli","SM1","SM2","SM3","Pr","En-A","En-C","En-C1",
         "En-V","En-p","En-SM","En-L","Nr","Neu","MC","Mon","M0","M1","M2",
         "M1-2","M-p","DC","P-DC","B","PC","T-cn","T-reg","T-rm","T-NK","T7",
         "T-ifn","T-int","T-p","T-un","RBC","Un")

colnames(df)[!(colnames(df) %in% cols)]
rownames(df)[!(rownames(df) %in% rows)]

df = df[rows,cols]
df = df[complete.cases(df),]
#write.csv(df, paste0(path,"Lung_24-",con,"_cell.types_by_samples.csv"))
write.csv(df, paste0(path,"Cell_types_by_samples_doublets.csv"))

df1 = read.csv(paste0(path,"Cell_types_by_samples.csv"),row.names = 1)
df1$cell.types = NULL
df2 = read.csv(paste0("output/20200626/","Cell_types_by_samples.csv"),row.names = 1)
df2$cell.types = NULL
table(rowSums(df1),rowSums(df2))

# recluster BC-IC-S
#sub_object = subset(object, idents = c("BC","BC-p","IC1","IC2","IC-S","S","S-d"))
sub_object = readRDS(file = "data/Lung_30_BC-IC-S_20200714.rds")
DefaultAssay(sub_object) = "SCT"

df_samples <- readxl::read_excel("doc/Annotations/Marker list for BC-IC-S clustering.xlsx",
                                 col_names = FALSE)
marker_genes <- CaseMatch(search = df_samples$...1, match = rownames(sub_object)) %>% as.character()
VariableFeatures(sub_object) = marker_genes
print(length(VariableFeatures(sub_object)))

sub_object %<>% ScaleData()
npcs <- 100
sub_object %<>% RunPCA(verbose = T,npcs = npcs, features = VariableFeatures(sub_object))
sub_object %<>% RunUMAP(reduction = "pca", dims = 1:npcs)
Idents(sub_object)= "annotations3"
UMAPPlot.1(sub_object,group.by = "annotations3", cols = rev(Singler.colors), 
           label = F, label.repel = F, no.legend = T,do.return = F, do.print = T)
UMAPPlot.1(sub_object,group.by = "annotations3", cols = rev(Singler.colors), 
           label = T, label.repel = T, no.legend = T,do.return = F, do.print = T)
saveRDS(sub_object, file = paste0("data/Lung_30_BC-IC-S_20200714.rds"))
sub_object = readRDS(file = paste0("data/Lung_30_BC-IC-S_20200714.rds"))

sub_meta.data = meta.data[colnames(sub_object),]
sub_object@meta.data$annotations3 = sub_meta.data$annotations3

PrepareShiny(sub_object, samples, Rshiny_path, split.by = "annotations3", 
             reduction = "umap",verbose = T)


# serial resolution
DefaultAssay(sub_object) = "SCT"
sub_object %<>% FindNeighbors(reduction = "pca",dims = 1:npcs)
save.path <- paste0(path,"serial_resolutions/")
if(!dir.exists(save.path))dir.create(save.path, recursive = T)
resolutions = c(seq(0.001,0.009, by = 0.001),seq(0.01,0.09, by = 0.01),seq(0.1,5, by = 0.1))
for(i in 1:length(resolutions)){
  sub_object %<>% FindClusters(resolution = resolutions[i])
  Idents(sub_object) = paste0("SCT_snn_res.",resolutions[i])
  UMAPPlot.1(sub_object, group.by=paste0("SCT_snn_res.",resolutions[i]),pt.size = 0.3,label = T,
             label.repel = T,alpha = 0.9,
             do.return = F,
             no.legend = T,label.size = 4, repel = T, 
             title = paste("res =",resolutions[i]),
             do.print = T, save.path = save.path)
  Progress(i,length(resolutions))
}
PrepareShiny(sub_object, samples, Rshiny_path, split.by = "SCT_snn_res.5", 
             reduction = "umap",verbose = T)
saveRDS(object, file = paste0("data/Lung_30_20200710.rds"))


# population expression
sub_object = readRDS(file = "data/Lung_30_BC-IC-S_20200714.rds")
DefaultAssay(sub_object) = "SCT"
exp = FetchData(sub_object, vars = c("annotations3","KRT5","KRT15","DST","LAMB3","S100A2","MIR205HG","SERPINB3",
                                     "GPX2","ALDH3A1","ADH7","UPK1B","FGFBP1","SERPINB4","SFN",
                                     "KRT6A","LY6D","SPRR2A","SPRR1B","SCGB1A1","SCGB3A1",
                                     "TFF3","MUC5B","BPIFB1","TSPAN8","MKI67","UBE2C","KIAA0101",
                                     "CENPW","TK1","TOP2A","SCGB3A2","SFTPB","MUC5AC"))
exp_list <- split(exp, exp$annotations3)
exp_list = lapply(exp_list, function(x) x[,-1])
exp_list[["BC"]] = NULL

write.xlsx(exp_list, file = paste0(path,"Expression matrix BC-p IC1 IC2 S.xlsx"),row.names=TRUE,
           colNames = TRUE, borders = "surrounding",colWidths = c(NA, "auto", "auto"))


# test ==================
cols <- c("SCT_snn_res.2","SCT_snn_res.3","SCT_snn_res.4","SCT_snn_res.4.1","SCT_snn_res.4.9","SCT_snn_res.5")
for(m in seq_along(cols)){
        object@meta.data[,cols[m]] %<>% as.character() %>% as.integer()
}
Idents(object)= "annotations3"
AT2_S <- subset(object, idents = c("AT2","S","S-d+","S-d++","S-d+++"))
df_list <- list()
for(m in seq_along(cols)){
        df_list[[m]] <- table(AT2_S$annotations3, AT2_S@meta.data[,cols[m]]) %>% as.data.frame.matrix
}
names(df_list) = cols
write.xlsx(df_list, file = paste0(path,"AT2_S_cell.types_cluters.xlsx"),
           colNames = TRUE,rowNames = TRUE, borders = "surrounding",colWidths = c(NA, "auto", "auto"))

Idents(object) = "SCT_snn_res.5"
C88 <- subset(object, idents = 88)
# density plot and gene plots
Idents(object) = "SCT_snn_res.2"
C11 <- subset(object, idents = 11)
dir.create(paste0(path,"Cluster 11"), recursive = T)
jpeg(paste0(path,"Cluster 11/DensityPlot_Cluster 11.jpeg"), units="in", width=10, height=10,res=600)
print(RidgePlot(C11, features = c("SFTPC", "SCGB1A1"), ncol = 1,group.by = "annotations3"))
dev.off()
jpeg(paste0(path,"Cluster 11/FeatureScatter_Cluster 11_SFTPC-SCGB1A1.jpeg"), units="in", width=10, height=10,res=600)
FeatureScatter(C11, feature1 = "SFTPC", feature2 = "SCGB1A1",group.by = "annotations3")
dev.off()
jpeg(paste0(path,"Cluster 11/FeatureScatter_Cluster 11_SCGB3A2-SCGB1A1.jpeg"), units="in", width=10, height=10,res=600)
FeatureScatter(C11, feature1 = "SCGB3A2", feature2 = "SCGB1A1",group.by = "annotations3")
dev.off()
jpeg(paste0(path,"Cluster 11/FeatureScatter_Cluster 11_SFTPB-SCGB1A1.jpeg"), units="in", width=10, height=10,res=600)
FeatureScatter(C11, feature1 = "SFTPB", feature2 = "SCGB1A1",group.by = "annotations3")
dev.off()
jpeg(paste0(path,"Cluster 11/FeatureScatter_Cluster 11_SCGB1A1-SCGB3A2.jpeg"), units="in", width=10, height=10,res=600)
FeatureScatter(C11, feature1 = "SCGB1A1", feature2 = "SCGB3A2",group.by = "annotations3")
dev.off()

C23 <- subset(object, idents = 23)
jpeg(paste0(path,"DensityPlot_Cluster 23.jpeg"), units="in", width=10, height=10,res=600)
print(RidgePlot(C23, features = c("SCGB3A2", "SFTPB"),group.by = "annotations3", ncol = 1))
dev.off()
jpeg(paste0(path,"FeatureScatter_Cluster 23_SCGB3A2-MUC5AC.jpeg"), units="in", width=10, height=10,res=600)
FeatureScatter(C23, feature1 = "SCGB3A2", feature2 = "MUC5AC",group.by = "annotations3",)
dev.off()
jpeg(paste0(path,"FeatureScatter_Cluster 23_SFTPB-MUC5AC.jpeg"), units="in", width=10, height=10,res=600)
FeatureScatter(C23, feature1 = "SFTPB", feature2 = "MUC5AC",group.by = "annotations3",)
dev.off()
jpeg(paste0(path,"FeatureScatter_Cluster 23_MUC5AC-MUC5B.jpeg"), units="in", width=10, height=10,res=600)
FeatureScatter(C23, feature1 = "MUC5AC", feature2 = "MUC5B",group.by = "annotations3",)
dev.off()

Idents(object) = "SCT_snn_res.4"
C40 <- subset(object, idents = 40)
jpeg(paste0(path,"DensityPlot_Cluster 40.jpeg"), units="in", width=10, height=10,res=600)
print(RidgePlot(C40, features = c("KRT14", "KRT5"),group.by = "annotations3", ncol = 1))
dev.off()
jpeg(paste0(path,"FeatureScatter_Cluster 40_KRT14-KRT5.jpeg"), units="in", width=10, height=10,res=600)
FeatureScatter(C40, feature1 = "KRT14", feature2 = "KRT5",group.by = "annotations3",)
dev.off()

Idents(object) = "annotations3"
Proliferating <- subset(object, idents = "Proliferating")
# label by corrdinates
g <- UMAPPlot.1(Proliferating, group.by = "annotations3", cols = Singler.colors, 
                label = T, label.repel = T, no.legend = T,do.return = T)+
        rectangle(-4.5, -1, 1, 4.5,colour = "black")+ # BC-p
        rectangle(-13, -6.5, 8, 13.5,colour = "black")+ # AT2-p
        rectangle(3, 12, 2.5, 10.5,colour = "black")+ # T-p
        rectangle(-1.5, 7, -14, -8,colour = "black")+ # En-p
        rectangle(-8.5, -1.75, -12, -2.5,colour = "black")+ # St-p
        rectangle(0, 7, -4.5, 1.75,colour = "black") # M-p

jpeg(paste0(path,"UMAP_all.jpeg"), units="in", width=10, height=7,res=600)
print(g)
dev.off()


sub_object = subset(object, idents = c("BC","IC1","IC2","IC-S","S"))
df <- table(sub_object$annotations3, sub_object$SCT_snn_res.2) %>% 
        as.data.frame()
colnames(df) = c("cell.types","samples","Freq")
df %<>% tidyr::spread("samples","Freq")
rownames(df) = df$cell.types
df=  df[,-1]
df = df[c("BC","IC1", "IC2", "IC-S","S"),]
colnames(df) = paste("cluster",colnames(df))
write.csv(df, paste0(path,"BC_IC_S_res2.csv"))

sub_object = subset(object, idents = c("S","BC","BC-p"))
sub_object$annotations3_res2 = paste(sub_object$annotations3, "cluster", sub_object$SCT_snn_res.2)
PrepareShiny(sub_object, samples, Rshiny_path, split.by = "annotations3_res2", 
             reduction = "umap",verbose = T)

df <- table(object$annotations3, object$SCT_snn_res.2) %>% 
        as.data.frame()
colnames(df) = c("cell.types","samples","Freq")
df %<>% tidyr::spread("samples","Freq")
rownames(df) = df$cell.types
df=  df[,-1]
colnames(df) = paste("cluster",colnames(df))
write.csv(df, paste0(path,"Celltypes_res2.csv"))

Idents(sub_object) = "annotations3_res2"
clusters <- c("S cluster 35", 
              "S cluster 1", 
              "S cluster 2", 
              "S cluster 4", 
              "S cluster 8", 
              "S cluster 24", 
              "S cluster 28", 
              "S cluster 60", 
              "S cluster 61", 
              "BC cluster 9", 
              "BC cluster 25", 
              "BC cluster 27", 
              "BC cluster 35", 
              "BC cluster 1", 
              "BC cluster 4", 
              "BC cluster 31", 
              "BC cluster 44", 
              "BC cluster 59", 
              "BC cluster 60", 
              "BC-p cluster 9", 
              "BC-p cluster 25", 
              "BC-p cluster 27", 
              "BC-p cluster 59", 
              "BC-p cluster 60")
sub_object <- subset(sub_object, idents = clusters)

exp = FetchData(sub_object, vars = c("annotations3_res2","KRT5", "KRT15", "S100A2", 
                                     "SERPINB3", "SERPINB4", "LY6D", "SCGB1A1", "MUC5AC", "MUC5B","MKI67","UBE2C"))
exp_list <- split(exp, exp$annotations3_res2)
exp_list = lapply(exp_list, function(x) x[,-1])

write.xlsx(exp_list, file = paste0(path,"B, BC and BC-p - cluster res 2 expression.xlsx"),row.names=TRUE,
           colNames = TRUE, borders = "surrounding",colWidths = c(NA, "auto", "auto"))

# label by barcodes 20200802 ============
df_annotation7 <- readxl::read_excel("doc/Annotations/20200807_M0-M1-M2 revision.xlsx",
                                     sheet = "Sheet1")
(dup_cell <- as.vector(df_annotation7$`Cells`[duplicated(df_annotation7$`Cells`)]))
dup_df = df_annotation7[df_annotation7$Cells %in% dup_cell,] %>% as.data.frame()
dup_df = dup_df[order(dup_df$`Cells`),]
(N = nrow(dup_df))
table(dup_df[seq(1,N,2),"Cells"] == dup_df[seq(1,N,2)+1,"Cells"]) # all duplicated cells only has at most one copy

ambigous <- which(dup_df[seq(1,N,2),"Change to"] != dup_df[seq(2,N+1,2),"Change to"])
dup_df = dup_df[dup_df$Cells %in% dup_df[ambigous*2,"Cells"],]
rownames(dup_df)= NULL
sub_object = subset(object, cells = unique(dup_df$Cells))
exp = FetchData(sub_object, vars = unique(c("KRT5", "KRT15", "S100A2", "SERPINB3", "SERPINB4",
                                            "LY6D", "SCGB1A1", "MUC5AC", "MUC5B","MKI67","UBE2C",
                                            "KRT5","KRT15","DST","LAMB3","S100A2","MIR205HG","SERPINB3",
                                            "GPX2","ALDH3A1","ADH7","UPK1B","FGFBP1","SERPINB4","SFN",
                                            "KRT6A","LY6D","SPRR2A","SPRR1B","SCGB1A1","SCGB3A1",
                                            "TFF3","MUC5B","BPIFB1","TSPAN8","MKI67","UBE2C","KIAA0101",
                                            "CENPW","TK1","TOP2A","SCGB3A2","SFTPB","MUC5AC")))
exp$Cells = rownames(exp)
dup_df = left_join(dup_df, exp, by = "Cells")

write.xlsx(dup_df, file = paste0(path,"Re-labeling 8-2-20 RS_duplicated.xlsx"),
           colNames = TRUE,rowNames = TRUE, borders = "surrounding",colWidths = c(NA, "auto", "auto"))

# T7 cluster contains some ciliated (C) cells and neutrophils:
T7 <- subset(object, idents= "T7")
exp = FetchData(T7, vars = c("FOXJ1", "CAPS", "TPPP3","SRGN","LAPTM5","DCN","S100A8","S100A9",
                             "CD3D", "CD3E", "IL7R", "CXCR4", "CCL5"))
barcode_list <- list()
barcode_list[[1]] = rownames(exp)[exp$FOXJ1 >0 | exp$CAPS >0 | exp$TPPP3 >0]
barcode_list[[2]] = rownames(exp)[exp$SRGN >0 | exp$LAPTM5 >0 | exp$DCN >0 ]
barcode_list[[3]] = rownames(exp)[exp$S100A8 >0 | exp$S100A9 >0]
barcode_list[[4]] = rownames(exp)[exp$CD3D >0 | exp$CD3E >0 | exp$IL7R >0 | exp$CXCR4 >0 | exp$CCL5 >0]
names(barcode_list) = c("Ciliated","immune and fibroblast","Neutrophils","T7")
euler_df <- eulerr::euler(barcode_list,shape = "circle")
g <- plot(euler_df, quantities = TRUE)
jpeg(paste0(path,"T7~.jpeg"), units="in", width=10, height=10,res=600)
print(g)
dev.off()

write.csv(exp, paste0(path,"T7_barcodes.csv"), row.names = T)

Idents(object) = "SCT_snn_res.3"
C11 <- subset(object, idents = 11)
jpeg(paste0(path,"T7_RidgePlot.jpeg"), units="in", width=10, height=10,res=600)
RidgePlot(T7, features = "7SK.2")
dev.off()

jpeg(paste0(path,"C11_RidgePlot.jpeg"), units="in", width=10, height=10,res=600)
RidgePlot(C11, features = "7SK.2", group.by = "annotations3")+xlim(0,7.5)
dev.off()
exp = FetchData(C11, vars = c("orig.ident","annotations3","7SK.2"))
table(meta.data[df_annotation6$`Cells`,"SCT_snn_res.3"]) %>% as.data.frame()

table(exp$annotations3, exp$`7SK.2` > 0)
table(exp$annotations3, exp$`7SK.2` > 0 & exp$`7SK.2` <= 1 )
table(exp$annotations3, exp$`7SK.2` > 1 & exp$`7SK.2` <= 2 )
table(exp$annotations3, exp$`7SK.2` > 2 & exp$`7SK.2` <= 3)
table(exp$annotations3, exp$`7SK.2` > 3)


table(meta.data[df_annotation6$`Cells`,"annotations3"], 
      ((rownames(meta.data) %in% df_annotation6$Cells) & exp[,"7SK.2"] >2 )) %>% 
        as.data.frame()
