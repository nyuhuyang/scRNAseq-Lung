library(Seurat)
library(dplyr)
library(tidyr)
library(kableExtra)
library(magrittr)
library(readxl)
library(cowplot)
library(stringr)
library(openxlsx)
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)
# load data
object = readRDS(file = "data/Lung_SCT_time6_20210908.rds")
DefaultAssay(object) = "SCT"
meta.data = object@meta.data
#======== rename ident =================
df_annotation <- readxl::read_excel("doc/Annotations/Exp 12 time-course d0-d28 - annotations.xlsx",
                                    sheet = "Sheet1")
colnames(df_annotation) %<>% gsub("res ","SCT_snn_res.",.)
(resolutions <- grep("SCT_snn_res.",colnames(df_annotation),value = T))
res = resolutions[1]
keep = !is.na(pull(df_annotation[,res]))
meta.data[,"Cell_subtype"] = plyr::mapvalues(meta.data[,res],
                                             from = pull(df_annotation[keep,res]),
                                             to = pull(df_annotation[keep,"Subtype"])
                                             )
meta.data[,"Cell_subtype"] %<>% as.character()
for(i in 2:length(resolutions)){
    keep = which(!is.na(pull(df_annotation[,resolutions[i]])))
        for(m in keep){
            cl = pull(df_annotation[m,resolutions[i]])
            change_to = pull(df_annotation[m,"Subtype"])
            meta.data[meta.data[,resolutions[i]] %in% cl,"Cell_subtype"] = change_to
            print(paste (resolutions[i],"at",cl,"------->",change_to))
            }
}

# add color
df_color = t(data.frame(
    c("AT1","#C946D4"),
    c("AT2","#A794D7"),
    c("B","#B4C7E7"),
    c("BC","#70AD47"),
    c("C-s","#FFE699"),
    c("C1","#FFC000"),
    c("CD8-T1","#5B9BD5"),
    c("cDC","#0070C0"),
    c("Cr","#B8C6DA"),
    c("En-a","#00B0F0"),
    c("En-c1","#A7E5BC"),
    c("En-ca","#BDD7EE"),
    c("En-l","#7CAFDD"),
    c("En-SM","#FFF2CC"),
    c("En-v","#C7BCE7"),
    c("Fb1","#FFA919"),
    c("Fb2","#C00000"),
    c("Fb3","#FA7FA9"),
    c("Fb4","#77D900"),
    c("G-Muc","#00B050"),
    c("G-Ser","#ADCDEA"),
    c("Gli","#8FAADC"),
    c("H","#EBE621"),
    c("IC","#FBE5D6"),
    c("Ion","#0F23FF"),
    c("M1","#FF6600"),
    c("M1-2","#D29F26"),
    c("M2","#BFBFBF"),
    c("MC","#2FE6F9"),
    c("ME","#E169CD"),
    c("Mon","#FBC8EF"),
    c("NE","#FF0000"),
    c("Neu","#FF7C88"),
    c("NK","#FF8B6A"),
    c("p-C","#EB6C11"),
    c("PC","#ED7D31"),
    c("pDC","#E31A1C"),
    c("Pr","#747474"),
    c("S-Muc","#B38600"),
    c("S1","#FFBAD1"),
    c("SM1","#C8D8E0"),
    c("SM2","#0087B6"),
    c("SM3","#FFFF00"),
    c("T-ifn","#FD5A00"),
    c("TASC","#FF3990"),
    c("Tcn","#9EF971"),
    c("T-NK","#FF42A4"),
    c("Trm","#f1c232"),#fc8e66
    c("T-un","#b1bcc5"),
    c("Un","#DEEBF7"),
    c("BC-m","#A6CEE3"),
    c("BC-p","#ce7e00"),
    c("BC-s","#7570B3"),
    c("IC-sq","#72eb3b"),
    c("IC1","#FBE5D5")
))#
df_color %<>% as.data.frame
rownames(df_color) = NULL
colnames(df_color) = c("Cell_subtype","Cell_subtype.colors")
table(duplicated(df_color$Cell_subtype.colors))
df_color$Cell_subtype.colors[duplicated(df_color$Cell_subtype.colors)]
meta.data$Cell_subtype.colors =  plyr::mapvalues(meta.data$Cell_subtype,
                                            from = df_color$Cell_subtype,
                                            to = df_color$Cell_subtype.colors)


Cell_types <- c("Subtype","Cell_type","Major_Cell_type")

df_annotation = df_annotation[order(df_annotation$Subtype),Cell_types]
df_annotation = df_annotation[!duplicated(df_annotation$Subtype),]
for(Cell_type in Cell_types[2:3]){
    meta.data[,Cell_type] = plyr::mapvalues(meta.data$Cell_subtype,
                                            from = pull(df_annotation[,"Subtype"]),
                                            to = pull(df_annotation[,Cell_type]))
}

#=========
table(rownames(object@meta.data) == rownames(meta.data))
#colnames(object@meta.data) %<>% gsub("cell_types","old_cell_types",.)
#object@meta.data %<>% cbind(meta.data)
object$Cell_subtype = meta.data$Cell_subtype
object$Cell_subtype.colors = meta.data$Cell_subtype.colors
object$Cell_type = meta.data$Cell_type
object$Major_Cell_type = meta.data$Major_Cell_type

meta.data = object@meta.data
saveRDS(meta.data, "output/20211020/meta.data_Cell_subtype_time6.rds")
object@meta.data = meta.data
saveRDS(object, file = "data/Lung_SCT_time6_20210908.rds")
