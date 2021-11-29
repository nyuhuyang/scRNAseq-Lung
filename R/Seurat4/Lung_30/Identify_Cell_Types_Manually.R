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
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_differential_expression.R")
source("shinyApp/Lung_30_hg38/util_palette.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)
# load data
object = readRDS(file = "data/Lung_SCT_30_20210831.rds")
DefaultAssay(object) = "SCT"
#meta.data = object@meta.data
#umap = object[["umap"]]@cell.embeddings
meta.data_exp = FetchData(object, vars = c("SFTPC","SCGB1A1","SCGB3A2","SFTPB","CCL19","SFRP2","GJA5","DKK2"))
#saveRDS(umap, "output/20210901/umap_dist.0.3_spread.1.rds")
#saveRDS(meta.data_exp, "output/20210901/meta.data_exp.rds")


min.dist = 0.5
spread = 1.4
file.name = paste0("dist.",min.dist,"_spread.",spread)
meta.data = readRDS(paste0("output/20210901/meta.data_",file.name,".rds"))
colnames(meta.data) %<>% gsub(file.name,"SCT_snn_res",.)
barcode = rownames(meta.data)
meta.data %<>% sapply(function(x) as.character(x)) %>% as.data.frame
rownames(meta.data) = barcode

umap = readRDS("output/20210901/umap_dist.0.3_spread.1_orig.rds")
table(rownames(meta.data) == rownames(umap$umap@cell.embeddings))
meta.data %<>% cbind(umap$umap@cell.embeddings)
meta.data %<>% cbind(meta.data_exp)

#======== rename ident =================
df_annotation <- readxl::read_excel("doc/Annotations/20210917_20UMAP res0.8 annotations.xlsx",
                                    sheet = "Sheet1")
resolutions = paste0("SCT_snn_res.",c(0.8,2,4,5))
res = resolutions[1]
keep = !is.na(pull(df_annotation[,res]))
meta.data[,"Cell_subtype"] = plyr::mapvalues(meta.data[,res],
                                             from = pull(df_annotation[keep,res]),
                                             to = pull(df_annotation[keep,"Cell_subtype"])
                                             )
for(i in 1:length(resolutions)){
    keep = which(!is.na(pull(df_annotation[,resolutions[i]])))
        for(m in keep){
            cl = pull(df_annotation[m,resolutions[i]])
            change_to = pull(df_annotation[m,"Cell_subtype"])
            meta.data[meta.data[,resolutions[i]] %in% cl,"Cell_subtype"] = change_to
            print(paste (resolutions[i],"at",cl,"------->",change_to))
            }
}

keep = !is.na(df_annotation$modify_condition)
df_annotation_UMAP = df_annotation[keep,c("SCT_snn_res.0.8","modify_SCT_snn_res.2","modify_condition","modify_label")]
colnames(df_annotation_UMAP) %<>% gsub("modify_","",.)

for(i in 1:length(resolutions[1:2])){
    keep = which(!is.na(pull(df_annotation_UMAP[,resolutions[i]])))
    for(m in keep){
        cl = pull(df_annotation_UMAP[m,resolutions[i]])
        change_to = pull(df_annotation_UMAP[m,"label"])
        
        select_id = meta.data %>% dplyr::filter(!!as.name(resolutions[i]) %in% cl) %>% 
                              dplyr::filter(eval(parse(text = df_annotation_UMAP$condition[m])))
        meta.data[rownames(select_id),"Cell_subtype"] = change_to
        print(paste (resolutions[i],"at",cl,df_annotation_UMAP$condition[m],"------->",change_to))
    }
}

df_annotation = df_annotation %>% dplyr::filter(!is.na(Cell_subtype)) %>% dplyr::filter(!duplicated(Cell_subtype))

Cell_types <- c("Cell_type","UMAP_land","Family","Superfamily")
for(Cell_type in Cell_types){
    meta.data[,Cell_type] = plyr::mapvalues(meta.data$Cell_subtype,
                                            from = pull(df_annotation[,"Cell_subtype"]),
                                            to = pull(df_annotation[,Cell_type]))
}

#======== rename ident =================
df_annotation <- readxl::read_excel("doc/Annotations/Annotation adjustments 10-4-21 RS.xlsx",
                                    sheet = "Sheet1")
colnames(df_annotation) %<>% gsub("modify_","",.)
colnames(meta.data) %<>% gsub("SCT_snn_res.","X20UMAP_res_",.)
for(i in 1:nrow(df_annotation)){
    change_from = pull(df_annotation[i,"Cell_subtype"])
    change_to = pull(df_annotation[i,"label"])
    print(paste ("If",df_annotation$Cell_subtype[i],df_annotation$condition[i],"------->",change_to))
    
    select_id = meta.data %>% dplyr::filter(Cell_subtype == change_from) %>% 
        dplyr::filter(eval(parse(text = df_annotation$condition[i])))
    meta.data[rownames(select_id),"Cell_subtype"] = change_to
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
    c("Un","#DEEBF7")))#
df_color %<>% as.data.frame
rownames(df_color) = NULL
colnames(df_color) = c("Cell_subtype","Cell_subtype.colors")
table(duplicated(df_color$Cell_subtype.colors))

meta.data$Cell_subtype.colors =  plyr::mapvalues(meta.data$Cell_subtype,
                                                 from = df_color$Cell_subtype,
                                                 to = df_color$Cell_subtype.colors)

df_annotation <- readxl::read_excel("doc/Annotations/20210917_20UMAP res0.8 annotations.xlsx",
                                    sheet = "Sheet1")
Cell_types <- c("Cell_subtype","Cell_type","UMAP_land","Family","Superfamily")

df_annotation = df_annotation[order(df_annotation$Cell_subtype),Cell_types]
df_annotation = df_annotation[!duplicated(df_annotation$Cell_subtype),]
for(Cell_type in Cell_types[2:5]){
    meta.data[,Cell_type] = plyr::mapvalues(meta.data$Cell_subtype,
                                            from = pull(df_annotation[,"Cell_subtype"]),
                                            to = pull(df_annotation[,Cell_type]))
}

#=========
table(rownames(object@meta.data) == rownames(meta.data))
#colnames(object@meta.data) %<>% gsub("cell_types","old_cell_types",.)
#object@meta.data %<>% cbind(meta.data)
object$Cell_subtype = meta.data$Cell_subtype
object$Cell_subtype.colors = meta.data$Cell_subtype.colors
object$Cell_type = meta.data$Cell_type
object$UMAP_land = meta.data$UMAP_land
object$Family = meta.data$Family
object$Superfamily = meta.data$Superfamily
meta.data = object@meta.data
saveRDS(meta.data, "output/20211004/meta.data_Cell_subtype.rds")
object@meta.data = meta.data
saveRDS(object, file = "data/Lung_SCT_30_20210831.rds")


# run Mann-Whitney test to test Differences in TASC% between D and COPD
df_TASC <- readxl::read_excel("output/20210917/Region specific 30 normal lung TASC percentage.xlsx")
D <- c("CU_12_D", "CU_12_D_R", "UNC_44_D", "UNC_48_D", "UNC_54_D", "UNC_55_D", "UNC_57_D", "UNC_66_D","UNC_67_D", "UNC_69_D", "UNC_70_D", "UNC_71_D", "VU_27_D")
COPD <- c("UNC_51_D", "UNC_52_D", "UNC_61_D", "VU_19_D", "VU_37_D")
table(c(D,COPD) %in% colnames(df_TASC))

y = as.vector(t(df_TASC[3,D]))
x = as.vector(t(df_TASC[3,COPD]))
(ASE <- wilcox.test(y,x,correct=FALSE))

y = as.vector(t(df_TASC[8,D]))
x = as.vector(t(df_TASC[8,COPD]))
(ASE <- wilcox.test(y,x,correct=FALSE))

y = as.vector(t(df_TASC[13,D]))
x = as.vector(t(df_TASC[13,COPD]))
(ASE <- wilcox.test(y,x,correct=FALSE))


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
object = readRDS(file = "data/Lung_SCT_30_20210831.rds")
DefaultAssay(object) = "SCT"
object %<>% subset(subset = Doublets %in% "Singlet")



SCG_exp = FetchData(object, vars = c("SCGB1A1","SCGB3A2"))
SCG_exp$SCGB1A1 = plyr::mapvalues(as.character(SCG_exp$SCGB1A1 >= 2),
                                  from = c("TRUE","FALSE"),
                                  to = c("SCGB1A1-hi","SCGB1A1-lo"))
SCG_exp$SCGB3A2 = plyr::mapvalues(as.character(SCG_exp$SCGB3A2 > 0),
                                  from = c("TRUE","FALSE"),
                                  to = c("SCGB3A2+","SCGB3A2-"))
colnames(SCG_exp) %<>% paste0("_lvl")
#object@meta.data %<>% cbind(SCG_exp)
SCG_exp$SCGB1A1_lvl %<>% paste0(" ")
SCG_exp$SCGB3A2_lvl %<>% paste0(" ")
object$Cell_subtype.SCGB1A1 = object$Cell_subtype
object$Cell_subtype.SCGB3A2 = object$Cell_subtype


SCG_exp[object$Cell_subtype != "TASC","SCGB1A1_lvl"] = ""
SCG_exp[object$Cell_subtype != "TASC","SCGB3A2_lvl"] = ""
object$Cell_subtype.SCGB1A1 %<>% paste0(SCG_exp$SCGB1A1_lvl,.)
object$Cell_subtype.SCGB3A2 %<>% paste0(SCG_exp$SCGB3A2_lvl,.)

TASC <- subset(object, subset = Cell_subtype %in% "TASC" & 
                   Doublets %in% "Singlet" &
                   Regions %in% c("distal","terminal"))
TASC %<>% FindNeighbors(reduction = "umap",dims = 1:2)
TASC %<>% FindClusters(resolution = 0.06)
UMAPPlot.1(TASC,do.print = T)
TASC_exp = FetchData(TASC, vars = c("SCGB1A1","SCGB3A2"))
TASC_exp$SCGB1A1 = plyr::mapvalues(as.character(TASC_exp$SCGB1A1 >= 2),
                                  from = c("TRUE","FALSE"),
                                  to = c("SCGB1A1-hi","SCGB1A1-lo"))
TASC_exp$SCGB3A2 = plyr::mapvalues(as.character(TASC_exp$SCGB3A2 > 0),
                                  from = c("TRUE","FALSE"),
                                  to = c("SCGB3A2+","SCGB3A2-"))
table(TASC_exp$SCGB1A1 > 5, TASC_exp$SCGB3A2 > 2)

FeaturePlot.1(TASC,features = "SCGB3A2",do.print = T)
TASC_marker <- FindAllMarkers_UMI(TASC, assay = "SCT",
                                  group.by = "SCT_snn_res.0.06",
                                  logfc.threshold = 0.5,
                                  only.pos = T,
                                  test.use = "MAST",
                                  latent.vars = "nFeature_SCT")
openxlsx::write.xlsx(TASC_marker, file =  paste0(path,"20211124_TASC_4clusters.xlsx"),
                     colNames = TRUE,row.names = T,borders = "surrounding",colWidths = c(NA, "auto", "auto"))
Top_n = 25
TASC_marker %<>% filter(!grepl("^RPL.*|^RPS.*|^MT-.*", gene))
top = TASC_marker %>% group_by(cluster) %>% top_n(Top_n, avg_log2FC)
table(top$cluster)
TASC %<>% ScaleData(features=top$gene)
featuresNum <- make.unique(top$gene, sep = ".")
TASC %<>% MakeUniqueGenes(features = top$gene)

TASC$SCT_snn_res.0.06 %<>% factor(levels = 0:3)
Idents(TASC) = "SCT_snn_res.0.06"
table(Idents(TASC))
TASC$Regions %<>% droplevels()
DoHeatmap.2(TASC, features = featuresNum, Top_n = Top_n,
            do.print=T, angle = 0,
            group.by = c("SCT_snn_res.0.06","Regions"),group.bar = T,
            group1.colors = hue_pal()(4),
            group2.colors = color_generator("NEJM",4)[4:3],
            title.size = 16, no.legend = F,size=10,hjust = 0.5,
            group.bar.height = 0.02, label=F, cex.row= 7,
            width=10, height=13,
            pal_gsea = FALSE,
            file.name = "Heatmap_top25_X4clusters~.jpeg",
            title = "Top 25 DE genes in 4 TASC clusters",
            save.path = path)
