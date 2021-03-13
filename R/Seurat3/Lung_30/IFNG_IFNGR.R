invisible(lapply(c("Seurat","dplyr","magrittr","tidyr","openxlsx",#"Seurat","MAST","future",
                   "gplots"), function(x) {
                           suppressPackageStartupMessages(library(x,character.only = T))
                   }))
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)

set.seed(101)

"""
Violin plots:

IFNG expression in immune cell types (order as in our dot-plot you generated for single-cell markers) 
IFNGR1 expression in surface airway epithelial cell types (order as listed in my email below for GSEA)
IFNGR2 expression in surface airway epithelial cell types (order as listed in my email below for GSEA)

4 variants (v1-v4) for each plot:

v1 - for all samples not divided into groups other than cell types (highlight cell types with significant enrichment; adj p<0.05 vs other cell types) 
v2 - grouped per P and D side-by-side for each cell type listed above (highlight cell types with significant difference; adj p<0.05 P vs D) 
v3 - grouped per D and COPD for each cell type listed above (highlight cell types with significant difference; adj p<0.05 D vs COPD) 
v4 - grouped per D-young and D-old for each cell type listed above (highlight cell types with significant difference; adj p<0.05 D-young vs D-old) 

Excel data: statistics for each of the 4 variants above: 

average UMI, Pct1 % positive cells, total number of cells of a given cell type 
log2FC and adj p value for comparisons:

v1 - compared to other cell types in the family shown in the plot 
v2 - between P and D
v3 - between COPD and D
v4 - between D-old and D-young
"""

object = readRDS(file = "data/Lung_30_20200710.rds")
SAE_cells = c("BC1","BC2","BC-p","IC1","IC2","IC3","S","d-S",
           "H","p-C","C1","C2","C3","Ion","NE")
Immune_cells = c("MC","Neu","Mon","M0","M1","M1-2","M2","M-p","DC","p-DC",
           "B","PC","T-cn","T-reg","T-int","T-rm","T-NK","T-ifn","T-p")
age_df = data.frame("orig.ident" = c("UNC-54-D", "UNC-57-D", "UNC-66-D", "UNC-70-D",
                                     "UNC-44-D", "UNC-48-D", "UNC-55-D", "UNC-67-D", "UNC-69-D", "UNC-71-D", "VU-27-D"),
                    "age" = c(rep("older",4),rep("younger",7)))

Idents(object) = "cell_types"
Immune <- subset(object, idents = Immune_cells)
Immune$cell_types %<>% factor(levels = Immune_cells)
Immune %<>% ScaleData("IFNG")
SAE <- subset(object, idents = SAE_cells)
SAE$cell_types %<>% factor(levels = SAE_cells)
SAE %<>% ScaleData(features = c("IFNGR1","IFNGR2"))
# ============= Immune ==================
jpeg(paste0(path,"IFNG_Immune_v1.jpeg"), units="in", width=10, height=7,res=600)
VlnPlot(Immune, features = "IFNG",group.by = "cell_types",
        assay = "RNA",pt.size = 0.01)  + 
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

##################  subset ################## 
Idents(Immune) = "conditions"
P_D <- subset(Immune,idents = c("proximal","distal"))
P_D$conditions %<>% factor(levels = c("proximal","distal"))
D_COPD <- subset(Immune,idents = c("distal","COPD"))
D_COPD$conditions %<>% factor(levels = c("distal","COPD"))
D <- subset(Immune,idents = "distal")
D$age <- plyr::mapvalues(D$orig.ident, from = age_df$orig.ident, to =age_df$age)
Idents(D) = "age"
D %<>% subset(idents = c("younger","older"))
D$age %<>% factor(levels = c("younger","older"))

jpeg(paste0(path,"IFNG_Immune_v2.jpeg"), units="in", width=10, height=7,res=600)
VlnPlot(P_D, group.by = "cell_types",features = "IFNG",
        split.by = "conditions",assay = "RNA",pt.size = 0.01)  + 
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

jpeg(paste0(path,"IFNG_Immune_v3.jpeg"), units="in", width=10, height=7,res=600)
VlnPlot(D_COPD, group.by = "cell_types",features = "IFNG",
        split.by = "conditions",assay = "RNA",pt.size = 0.01)  + 
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

jpeg(paste0(path,"IFNG_Immune_v4.jpeg"), units="in", width=10, height=7,res=600)
VlnPlot(D, group.by = "cell_types",features = "IFNG",split.by = "age",
        assay = "RNA",pt.size = 0.01)  + 
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
AverageExpression_pct <- function(object, group.by = "cell_types", 
                                  features = "IFNG", assays = "SCT"){
        UMI_df <- FetchData(object = object, vars = c(features,group.by),)
        UMI_df$IFNG %<>% expm1
        UMI <- group_by(UMI_df,cell_types) %>% 
                summarise_at(vars(features),list(average_UMI = mean,
                                                 p_positive_cells = function(x) sum(x >0)/length(x),
                                                 n_positive_cells = function(x) sum(x >0)))
        UMI$average_UMI %<>% log1p()
        rownames(UMI) = UMI$cell_types
        return(UMI)
}

Idents(Immune) = "conditions"
P <- subset(Immune,idents = c("proximal"))
D <- subset(Immune,idents = c("distal"))
COPD <- subset(Immune,idents = c("COPD"))
D$age <- plyr::mapvalues(D$orig.ident, from = age_df$orig.ident, to =age_df$age)
Idents(D) = "age"
younger<- subset(D,idents = "younger")
older<- subset(D,idents = "older")
res_list <- list()
UMI_list1 <- pbapply::pblapply(list(Immune,P,D,younger,older,COPD), 
                              FUN = AverageExpression_pct,features = "IFNG")
res_list[["IFNG immune"]] = bind_cols(UMI_list1)
################## surface airway epithelial ################## 
Idents(SAE) = "cell_types"
jpeg(paste0(path,"IFNGR12_SAE_v1.jpeg"), units="in", width=10, height=7,res=600)
VlnPlot(SAE, features = c("IFNGR1","IFNGR2"),group.by = "cell_types", 
        assay = "RNA",pt.size = 0)
dev.off()

# subset
Idents(SAE) = "conditions"
P_D <- subset(SAE,idents = c("proximal","distal"))
P_D$conditions %<>% factor(levels = c("proximal","distal"))
D_COPD <- subset(SAE,idents = c("distal","COPD"))
D_COPD$conditions %<>% factor(levels = c("distal","COPD"))
D <- subset(SAE,idents = "distal")
D$age <- plyr::mapvalues(D$orig.ident, from = age_df$orig.ident, to =age_df$age)
Idents(D) = "age"
D %<>% subset(idents = c("younger","older"))
D$age %<>% factor(levels = c("younger","older"))


jpeg(paste0(path,"IFNGR12_SAE_v2.jpeg"), units="in", width=10, height=7,res=600)
VlnPlot(P_D, group.by = "cell_types",features = c("IFNGR1","IFNGR2"),
        assay = "RNA",split.by = "conditions",pt.size = 0)
dev.off()

jpeg(paste0(path,"IFNGR12_SAE_v3.jpeg"), units="in", width=10, height=7,res=600)
VlnPlot(D_COPD, group.by = "cell_types",features = c("IFNGR1","IFNGR2"),
        assay = "RNA",split.by = "conditions",pt.size = 0)
dev.off()

jpeg(paste0(path,"IFNGR12_SAE_v4.jpeg"), units="in", width=10, height=7,res=600)
VlnPlot(D, group.by = "cell_types",features = c("IFNGR1","IFNGR2"),
        assay = "RNA",split.by = "age",pt.size = 0)
dev.off()

Idents(SAE) = "conditions"
P <- subset(SAE,idents = c("proximal"))
D <- subset(SAE,idents = c("distal"))
COPD <- subset(SAE,idents = c("COPD"))
D$age <- plyr::mapvalues(D$orig.ident, from = age_df$orig.ident, to =age_df$age)
Idents(D) = "age"
younger<- subset(D,idents = "younger")
older<- subset(D,idents = "older")

UMI_list2 <- pbapply::pblapply(list(SAE,P,D,younger,older,COPD), 
                              FUN = AverageExpression_pct,features = "IFNGR1")
UMI_list3 <- pbapply::pblapply(list(SAE,P,D,younger,older,COPD), 
                              FUN = AverageExpression_pct,features = "IFNGR2")
res_list[["IFNGR1 surface airway epi"]] = list2df(UMI_list2)
res_list[["IFNGR2 surface airway epi"]] = list2df(UMI_list3)

openxlsx::write.xlsx(res_list, file =  "Yang/Lung_30/Cell_Phone_DB/IFNG_IFNGR1/IFNG-IFNGR data raw.xlsx",
                     colNames = TRUE,row.names = F,borders = "surrounding",colWidths = c(NA, "auto", "auto"))


######################### p-value ###########################
deg_list <- pbapply::pblapply(Immune_cells, function(x){
        deg <- readxl::read_excel("Yang/Lung_30/DE_analysis/groups/DE_results_Immune.xlsx",
                                  sheet = x)
        deg[deg$gene %in% "IFNG",]
})
IFNG <- bind_rows(deg_list)

table(Immune$cell_types)%>% prop.table *100 %>% round(digits = 3)

deg_list <- pbapply::pblapply(SAE_cells, function(x){
        deg <- readxl::read_excel("Yang/Lung_30/DE_analysis/groups/DE_results_Surface Airway Epithelial.xlsx",
                                  sheet = x)
        deg[deg$gene %in% c("IFNGR1","IFNGR2"),]
})
IFNGR <- bind_rows(deg_list)
IFNGR$cell_type = gsub(" .*","",IFNGR$cluster)
table(SAE$cell_types) %>% prop.table *100 %>% round(digits = 3)

