########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################
invisible(lapply(c("Seurat","dplyr","cowplot",
                   "magrittr","data.table","future","ggplot2","tidyr"), function(x) {
                           suppressPackageStartupMessages(library(x,character.only = T))
                   }))

# ######################################################################
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)

#======1.2 load  Seurat =========================
# load files
object = readRDS(file = "data/Lung_SCT_30_20200710.rds")
DefaultAssay(object) = "SCT"
Idents(object) = "Doublets"
object <- subset(object, idents = "Singlet")
Idents(object) = "cell_types"
TASC <- subset(object, idents = "TASC")
Idents(TASC) = "conditions"
TASC <- subset(TASC, idents =c("distal","terminal","COPD"))
#=====================
gene = "SCGB3A2"#"SCGB1A1"
gene_lvl = 0#2.5
df_gene <- FetchData(TASC, vars = test_gene)
df_gene$barcode = rownames(df_gene)
df_gene %<>% mutate(gene_lvl = ifelse(df_gene[,gene] > gene_lvl, "SCGB3A2+","SCGB3A2-"))#"high", "low"))
rownames(df_gene) = df_gene$barcode
head(df_gene)
TASC@meta.data %<>% cbind(df_gene)
Idents(TASC) = "gene_lvl"
system.time(gene_markers <- FindAllMarkers.UMI(TASC, 
                                                  logfc.threshold = 0.1, only.pos = T,
                                                  test.use = "MAST",
                                                  latent.vars = "nFeature_RNA",
                                                  p.adjust.methods = "fdr"))
colnames(gene_markers)[8] = paste(gene,"level")
write.csv(gene_markers,paste0(path,"TASC_",gene,"_degs.csv"))

(lvls = unique(TASC$gene_lvl))
markers_list <- list()
for(i in seq_along(lvls)){
        lvl = lvls[i]
        sub_TASC <- subset(TASC,idents = lvl)
        Idents(sub_TASC) = "conditions"
        tmp <- FindPairMarkers(object = sub_TASC,
                                           ident.1 = c("COPD","terminal"),
                                           ident.2 = c("distal","distal"),
                                           logfc.threshold = 0.1, only.pos = F,
                                           test.use = "MAST",
                                           latent.vars = "nFeature_RNA",
                                           p.adjust.methods = "fdr")
        tmp$cluster1.vs.cluster2 %<>% paste(lvl,.)
        tmp %<>% split(f = tmp$cluster1.vs.cluster2)
        markers_list[[i]] = tmp
}
markers = c(markers_list[[1]],markers_list[[2]])
names(markers) %<>% gsub("/","vs",.)
openxlsx::write.xlsx(markers, file =  paste0(path,"TASC_",gene,"_degs.xlsx"),
                     colNames = TRUE,row.names = T,borders = "surrounding",colWidths = c(NA, "auto", "auto"))

sub_TASC <- subset(TASC,idents = "low")
Idents(sub_TASC) = "conditions"
sub_TASC <- subset(sub_TASC,idents = c("distal","terminal"))
sub_TASC@meta.data %<>% select(!c("SCGB1A1","barcode","SCGB1A1_lvl"))
sub_TASC$cell_types = "SCGB1A1-low-TASC"

AT2 <- subset(object, idents = "AT2")
Idents(AT2) = "conditions"
AT2 <- subset(AT2, idents =c("distal","terminal"))

new_object <- merge(AT2, sub_TASC)

Idents(new_object) = "cell_types"
markers <- FindAllMarkers.UMI(new_object,
                              logfc.threshold = 0.1, only.pos = T,
                              test.use = "MAST",
                              latent.vars = "nFeature_RNA",
                              p.adjust.methods = "fdr")

colnames(markers)[8] = "cell.types"
write.csv(markers,paste0(path,"AT2_vs_TASC_low_SCGB1A1_degs.csv"))
