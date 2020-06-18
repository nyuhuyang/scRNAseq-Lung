library(Seurat)
library(dplyr)
library(tidyr)
library(kableExtra)
library(magrittr)
library(readxl)
library(cowplot)
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)
# load data
(load(file="data/Lung_28_Nointeg_20200131.Rda"))
#======== rename ident =================
object %<>% FindClusters(resolution = 3)
object %<>% FindClusters(resolution = 4)
object %<>% FindClusters(resolution = 1)
#resolution = 1
Idents(object) = "SCT_snn_res.1"
object %<>% RenameIdents("0" = "En-V",
                         "3" = "En-C",
                         "4" = "Mon",
                         "6" = "C1",
                         "7" = "Neu",
                         "9" = "B",
                         "11" = "NK",
                         "12" = "SM1",
                         "14" = "BC",
                         "15" = "C3",
                         "19" = "DC",
                         "24" = "MC",
                         "26" = "IC1",
                         "29" = "SM2",
                         "31" = "En-L",
                         "32" = "SM3",
                         "33" = "Pr",
                         "36" = "pre-C",
                         "37" = "Cr",
                         "38" = "PC",
                         "39" = "IC2",
                         "40" = "T-ifn",
                         "42" = "AT1",
                         "43" = "Nr",
                         "44" = "MEC",
                         "46" = "E-SM",
                         "16" = "AT2",
                         "20" = "AT2",
                         "23" = "C2",
                         "28" = "C2",
                         "35" = "C2",
                         "8" = "S",
                         "22" = "S",
                         "25" = "S",
                         "45" = "S",
                         "1" = "T",
                         "2" = "T",
                         "5" = "T",
                         "30" = "T")
#object[["cell.types"]] = as.character(Idents(object)) # Lung_28_harmony_rmD_20200205.Rda
object[["annotations2"]] = as.character(Idents(object))

exp = FetchData(object, vars = c("KRT5","CD163","CDH5","DCN","ACTA2","CD3E","SFTPC"))
c_34 <- object$SCT_snn_res.1 %in% 34
object@meta.data[c_34 & exp$KRT5 >0,"annotations2"] = "BC-p"
object@meta.data[c_34 & exp$CD163 >0,"annotations2"] = "M-p"
object@meta.data[c_34 & exp$CDH5 >0,"annotations2"] = "En-p"
object@meta.data[c_34 & exp$DCN >0,"annotations2"] = "F-p"
object@meta.data[c_34 & exp$ACTA2 >0,"annotations2"] = "SM-p"
object@meta.data[c_34 & exp$CD3E >0,"annotations2"] = "T-p"
object@meta.data[c_34 & exp$SFTPC >0,"annotations2"] = "AT-p"

#resolution = 3
object@meta.data[object$SCT_snn_res.3 %in% 31,"annotations2"] = "En-C1"
object@meta.data[object$SCT_snn_res.3 %in% c(10,16),"annotations2"] = "F1"
object@meta.data[object$SCT_snn_res.3 %in% 36,"annotations2"] = "F3"
object@meta.data[object$SCT_snn_res.3 %in% 46,"annotations2"] = "F2"
object@meta.data[object$SCT_snn_res.3 %in% 54,"annotations2"] = "En-A"
#resolution = 4
object@meta.data[object$SCT_snn_res.4 %in% c(9,45),"annotations2"] = "C1"
object@meta.data[object$SCT_snn_res.4 %in% 34,"annotations2"] = "SMG-Ser"
object@meta.data[object$SCT_snn_res.4 %in% 87,"annotations2"] = "SMG-Muc"
object@meta.data[object$SCT_snn_res.4 %in% 26,"annotations2"] = "M1"
object@meta.data[object$SCT_snn_res.4 %in% 38,"annotations2"] = "M1/2"
object@meta.data[object$SCT_snn_res.4 %in% c(23,74),"annotations2"] = "M2"
object@meta.data[object$SCT_snn_res.4 %in% 71,"annotations2"] = "S-d"
object@meta.data[object$SCT_snn_res.4 %in% 56,"annotations2"] = "T-reg"
object@meta.data[object$SCT_snn_res.4 %in% 7,"annotations2"] = "DC"
object@meta.data[object$SCT_snn_res.4 %in% 91,"annotations2"] = "P-DC2"
object@meta.data[object$SCT_snn_res.4 %in% 94,"annotations2"] = "Gli"
object@meta.data[object$SCT_snn_res.4 %in% c(5,39),"annotations2"] = "NK"
object@meta.data[object$SCT_snn_res.4 %in% c(15,24,25,30,31),"annotations2"] = "T-rm"
object@meta.data[object$SCT_snn_res.4 %in% c(12,20,77,18,56),"annotations2"] = "T-cn"
object@meta.data[object$SCT_snn_res.4 %in% 44,"annotations2"] = "T-0"
object@meta.data[object$SCT_snn_res.4 %in% 19,"annotations2"] = "T-7"

Idents(object)= "annotations2"

# select by corrdinates
jpeg(paste0(path,"UMAP.jpeg"), units="in", width=10, height=10,res=600)
UMAPPlot.1(object,group.by = "annotations2", cols = Singler.colors, 
           label = T, label.repel = T, no.legend = T,do.return = T)+
        rectangle(1, 3, -12.5, -11,colour = "black")+
        rectangle(1.5, 4, -15, -13,colour = "black")+
        rectangle(-4.5, -3, 0, 1,colour = "black")
dev.off()

meta.data = cbind(object@meta.data, object@reductions$umap@cell.embeddings)

Ion = meta.data[,"UMAP_1"] > 1 & meta.data[,"UMAP_1"] < 3 & 
        meta.data[,"UMAP_2"] > -12.5 & meta.data[,"UMAP_2"] < -11 & 
        meta.data[,"annotations2"] == 41
NEC = meta.data[,"UMAP_1"] > 1.5 & meta.data[,"UMAP_1"] < 4 & 
        meta.data[,"UMAP_2"] > -15 & meta.data[,"UMAP_2"] < -13 & 
        meta.data[,"annotations2"] == 41
PDC1 = meta.data[,"UMAP_1"] > -4.5 & meta.data[,"UMAP_1"] < -3 & 
        meta.data[,"UMAP_2"] > 0 & meta.data[,"UMAP_2"] < 1 & 
        meta.data[,"SCT_snn_res.4"] == 91
object@meta.data[rownames(meta.data)[Ion], "annotations2"] = "Ion"
object@meta.data[rownames(meta.data)[NEC], "annotations2"] = "NEC"
object@meta.data[rownames(meta.data)[PDC1], "annotations2"] = "P-DC1"


PrepareShiny(sub_object, samples, Rshiny_path, split.by = "annotations2", 
             reduction = "umap",verbose = T)

Idents(object) = "annotations2"
object %<>% sortIdent()
object <- AddMetaColor(object = object, label= "annotations2", colors = Singler.colors)
lapply(c(TRUE, FALSE), function(lab)
        UMAPPlot.1(sub_object, group.by="annotations2",pt.size = 0.5,label = lab,
                   cols = ExtractMetaColor(object),
                   label.repel = T, alpha= 0.9,
                   no.legend = T,label.size = 4, repel = T, title = "Annotations",
                   do.print = T, do.return = F))
sub_object <- subset(object, idents = c(10,13,17,18,21,27,34,41))

saveRDS(object, file = paste0("data/Lung_Nointeg_20200601.rds"))

# annotate unknown ==========
object <- readRDS("data/Lung_Nointeg_20200601.rds")
table(object$annotations2)
unknown <- subset(object, idents = c(10,13,17,18,21,27,34,41))
UMAPPlot.1(unknown, no.legend = T)
# select by corrdinates
jpeg(paste0(path,"unknown_UMAP.jpeg"), units="in", width=10, height=7,res=600)
UMAPPlot.1(unknown,group.by = "annotations2", cols = Singler.colors, 
           label = T, label.repel = T, no.legend = T,do.return = T)+
        rectangle(-9, -6.5, 8, 11,colour = "red")+ #Neu
        rectangle(-11, -9.5, 2.7, 4,colour = "purple")+ #T-p
        rectangle(-14, -11.5, -1, 3,colour = "dodgerblue")+ #NK
        rectangle(-11.5, -9.5, -1.5, 0.2,colour = "green")+#Trm
        rectangle(-11.5, -8, 0.2, 2.5,colour = "green")+ #Trm
        rectangle(-8, -6, 0.2, 2.5,colour = "red")+ #T7
        rectangle(-9.5, -6, -3,0.2,colour = "blue")+ #Tcn
        rectangle(-6, -3, 0, 4,colour = "lightskyblue")+ #B
        rectangle(-4.2, -2, -3, 0,colour = "red")+ #MC
        rectangle(-6, -2.5, -7.2, -3,colour = "darkgreen")+ #Mon
        rectangle(-13, -10.8, -9, -6,colour = "red")+ #PC
        rectangle(-7, -3.8, -9, -7.2,colour = "red")+ #M1/2
        rectangle(-7, -5.7, -9, -8,colour = "purple")+ #M-p
        rectangle(-3.8, -1, -12, -7.2,colour = "khaki")+ #M1
        rectangle(-7, -3.8, -12, -9,colour = "firebrick")+ #M2
        rectangle(-1, 1.8, 9.4, 14,colour = "red")+ #En-V
        rectangle(0.8, 1.8, 11.5, 13,colour = "purple")+ #En-p
        rectangle(-2, 1.8, 7.3, 9.4,colour = "red")+ #En-C
        rectangle(-2, 1, 6, 7.3,colour = "red")+ #En-C1
        rectangle(1.8, 3, 7, 14.5,colour = "darkgreen")+ #En-SM
        rectangle(5, 7, 16, 18,colour = "red")+ #Nr
        rectangle(3,7, 7.5, 13,colour = "red")+ #SM1
        rectangle(7,9 ,8.3, 11,colour = "red")+ #SM2
        rectangle(7,10 ,6.5, 8.3,colour = "red")+ #SM3
        rectangle(3,5.5, 4.5, 7.5,colour = "red")+ #F2
        rectangle(5.5,7, 4.5, 7.5,colour = "darkblue")+ #F1
        rectangle(3,7, 2, 4.5,colour = "darkblue")+ #F1
        rectangle(7,9, 2, 6.5,colour = "darkgreen")+ #F3
        rectangle(11,13, 4.5, 7.5,colour = "red")+ #En-L
        rectangle(14, 16, 4, 6,colour = "red")+ #Cr
        rectangle(0, 3, -2, 1.5, colour = "khaki")+ #BC
        rectangle(1, 2.8, -4.5, -1.8, colour = "red")+ #IC1
        rectangle(2.8, 4.2, -4, -1, colour = "green")+ #IC2
        rectangle(4.2, 8, -4, -0.5, colour = "red")+ #H
        rectangle(4, 6, -6, -4,colour = "red")+ #S
        rectangle(6, 8, -4.5, -3.5,colour = "blue")+ #S-d
        rectangle(7.5, 8.5, -6, -4.5,colour = "red")+ #H
        rectangle(9, 16, -5, 2,colour = "red")+ #AT2
        rectangle(1, 3.5, -12, -4.5,colour = "red")+ #SMG-Ser
        rectangle(5, 7.5, -12, -5.8,colour = "red")+ #C1
        rectangle(7.5, 9, -10, -6,colour = "darkblue")+ #C2
        rectangle(9, 12, -9, -5,colour = "darkblue")+ #C2
        rectangle(9,10.5, -11, -9,colour = "lightskyblue")+ #Pre-C
        rectangle(1, 3.5, -15, -13,colour = "red") #NEC

dev.off()


opts = data.frame(c(-9, -6.5, 8, 11, "Neu"),
                  c(-11, -9.5, 2.7, 4, "T-p"),
                  c(-14, -11.5, -1, 3, "NK"),
                  c(-11.5, -9.5, -1.5, 0.2,"Trm"),
                  c(-11.5, -8, 0.2, 2.5, "Trm"),
                  c(-8, -6, 0.2, 2.5, "T7"),
                  c(-9.5, -6, -3,0.2, "Tcn"),
                  c(-6, -3, 0, 4, "B"),
                  c(-4.2, -2, -3, 0, "MC"),
                  c(-6, -2.5, -7.2, -3, "Mon"),
                  c(-13, -10.8, -9, -6, "PC"),
                  c(-7, -3.8, -9, -7.2, "M1/2"),
                  c(-7, -5.7, -9, -8, "M-p"),
                  c(-3.8, -1, -12, -7.2, "M1"),
                  c(-7, -3.8, -12, -9, "M2"),
                  c(-1, 1.8, 9.4, 14, "En-V"),
                  c(0.8, 1.8, 11.5, 13, "En-p"),
                  c(-2, 1.8, 7.3, 9.4, "En-C"),
                  c(-2, 1, 6, 7.3, "En-C1"),
                  c(1.8, 3, 7, 14.5, "En-SM"),
                  c(5, 7, 16, 18, "Nr"),
                  c(3,7, 7.5, 13, "SM1"),
                  c(7,9 ,8.3, 11, "SM2"),
                  c(7,10 ,6.5, 8.3, "SM3"),
                  c(3,5.5, 4.5, 7.5, "F2"),
                  c(5.5,7, 4.5, 7.5, "F1"),
                  c(3,7, 2, 4.5, "F1"),
                  c(7,9, 2, 6.5, "F3"),
                  c(11,13, 4.5, 7.5, "En-L"),
                  c(14, 16, 4, 6, "Cr"),
                  c(0, 3, -2, 1.5, "BC"),
                  c(1, 2.8, -4.5, -1.8,  "IC1"),
                  c(2.8, 4.2, -4, -1,  "IC2"),
                  c(4.2, 8, -4, -0.5,  "H"),
                  c(4, 6, -6, -4, "S"),
                  c(6, 8, -4.5, -3.5, "S-d"),
                  c(7.5, 8.5, -6, -4.5, "H"),
                  c(9, 16, -5, 2, "AT2"),
                  c(1, 3.5, -14, -4.5, "SMG-Ser"),
                  c(5, 7.5, -12, -5.8, "C1"),
                  c(7.5, 9, -10, -6, "C2"),
                  c(9, 12, -9, -5, "C2"),
                  c(9,10.5, -11, -9, "Pre-C"),
                  c(1, 3.5, -15, -13, "NEC"),
                  stringsAsFactors = F) %>% t %>% as.data.frame()

colnames(opts) = c("x_left", "x_right", "y_bottom","y_top", "cell.types")
opts[,1:4] %<>% apply(2, as.numeric)

# rename unknown
meta.data = unknown@meta.data
for(i in 1:nrow(opts)){
        opt = opts[i,]
        cell.embeddings = unknown@reductions$umap@cell.embeddings
        keep_cells = cell.embeddings[,"UMAP_1"] > opt$x_left & cell.embeddings[,"UMAP_1"] < opt$x_right & 
                cell.embeddings[,"UMAP_2"] > opt$y_bottom & cell.embeddings[,"UMAP_2"] < opt$y_top
        
        keep_cells = names(keep_cells)[keep_cells]
        meta.data[keep_cells,"annotations2"] = as.character(opt$cell.types)
        Progress(i, nrow(opts))
}
unknown@meta.data = meta.data
Idents(unknown) = "annotations2"
unknown %<>% AddMetaColor(label= "annotations2", colors = Singler.colors)

object@meta.data[rownames(meta.data),"annotations2"] = meta.data$annotations2
Idents(object) = "annotations2"
object %<>% AddMetaColor(label= "annotations2", colors = Singler.colors)
UMAPPlot.1(object, label = T, label.repel = T, do.print = T, no.legend = T, do.return = F)
saveRDS(object, file = paste0("data/Lung_Nointeg_20200601.rds"))

# check unseen cell at different resolution,
opts <- list(list("SCT_snn_res.1", 0 , "En-V"),
          list("SCT_snn_res.1", 3 , "En-C"),
          list("SCT_snn_res.1", 4 , "Mon"),
          list("SCT_snn_res.1", 6 , "C1"),
          list("SCT_snn_res.1", 7 , "Neu"),
          list("SCT_snn_res.1", 9 , "B"),
          list("SCT_snn_res.1", 11 , "NK"),
          list("SCT_snn_res.1", 12 , "SM1"),
          list("SCT_snn_res.1", 14 , "BC"),
          list("SCT_snn_res.1", 15 , "C3"),
          list("SCT_snn_res.1", 19 , "DC"),
          list("SCT_snn_res.1", 24 , "MC"),
          list("SCT_snn_res.1", 26 , "IC1"),
          list("SCT_snn_res.1", 29 , "SM2"),
          list("SCT_snn_res.1", 31 , "En-L"),
          list("SCT_snn_res.1", 32 , "SM3"),
          list("SCT_snn_res.1", 33 , "Pr"),
          list("SCT_snn_res.1", 36 , "pre-C"),
          list("SCT_snn_res.1", 37 , "Cr"),
          list("SCT_snn_res.1", 38 , "PC"),
          list("SCT_snn_res.1", 39 , "IC2"),
          list("SCT_snn_res.1", 40 , "T-ifn"),
          list("SCT_snn_res.1", 42 , "AT1"),
          list("SCT_snn_res.1", 43 , "Nr"),
          list("SCT_snn_res.1", 44 , "MEC"),
          list("SCT_snn_res.1", 46 , "E-SM"),
          list("SCT_snn_res.1", 16 , "AT2"),
          list("SCT_snn_res.1", 20 , "AT2"),
          list("SCT_snn_res.1", 23 , "C2"),
          list("SCT_snn_res.1", 28 , "C2"),
          list("SCT_snn_res.1", 35 , "C2"),
          list("SCT_snn_res.1", 8 , "S"),
          list("SCT_snn_res.1", 22 , "S"),
          list("SCT_snn_res.1", 25 , "S"),
          list("SCT_snn_res.1", 45 , "S"),
          list("SCT_snn_res.1", 1 , "T"),
          list("SCT_snn_res.1", 2 , "T"),
          list("SCT_snn_res.1", 5 , "T"),
          list("SCT_snn_res.1", 30 , "T"),
          #resolution = 3
          list("SCT_snn_res.3", 31, "En-C1"),
          list("SCT_snn_res.3", c(10,16), "F1"),
          list("SCT_snn_res.3", 36, "F3"),
          list("SCT_snn_res.3", 46, "F2"),
          list("SCT_snn_res.3", 54, "En-A"),
          #resolution = 4
          list("SCT_snn_res.4", c(9,45), "C1"),
          list("SCT_snn_res.4", 34, "SMG-Ser"),
          list("SCT_snn_res.4", 87, "SMG-Muc"),
          list("SCT_snn_res.4", 26, "M1"),
          list("SCT_snn_res.4", 38, "M1-M2"),
          list("SCT_snn_res.4", c(23,74), "M2"),
          list("SCT_snn_res.4", 71, "S-d"),
          list("SCT_snn_res.4", 56, "T-reg"),
          list("SCT_snn_res.4", 7, "DC"),
          list("SCT_snn_res.4", 91, "P-DC2"),
          list("SCT_snn_res.4", 94, "Gli"),
          list("SCT_snn_res.4", c(5,39), "NK"),
          list("SCT_snn_res.4", c(15,24,25,30,31), "T-rm"),
          list("SCT_snn_res.4", c(12,20,77,18,56), "T-cn"),
          list("SCT_snn_res.4", 44, "T-0"),
             list("SCT_snn_res.4", 19, "T-7"))
for(m in 49:length(opts)){
        opt <- opts[[m]]
        Idents(object) = opt[[1]]
        sub_object <- subset(object, idents = opt[[2]])
        Idents(sub_object) = "annotations2"
        g <- UMAPPlot.1(sub_object, cols = ExtractMetaColor(sub_object),
                        no.legend = T, label = T, label.repel = T,
                        title = paste(opt[[3]], "in",opt[[1]], "cluster =",paste(opt[[2]], collapse = ",")),
                        do.print = F, do.return = T)
        g = g+ xlim(-15,17)+ylim(-17,20)
        
        jpeg(paste0(path,opt[[3]],"_",opt[[1]],"_",paste(opt[[2]], collapse = ","),".jpeg"), units="in", width=10, height=10,res=600)
        print(g)
        dev.off()
        Progress(m, length(opts))
}
