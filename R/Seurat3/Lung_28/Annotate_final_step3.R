########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################
invisible(lapply(c("Seurat","dplyr","kableExtra","cowplot",
                   "magrittr","data.table","future","ggplot2","tidyr"), function(x) {
                           suppressPackageStartupMessages(library(x,character.only = T))
                   }))
source("../R/Seurat3_functions.R")
path <- paste0("Yang/proximal_distal_terminal_COPD/Harmony/Expression/")
if(!dir.exists(path))dir.create(path, recursive = T)

#======== read csv =================
(csv_files <-list.files(path, pattern = "20200219.csv"))
csv_files = csv_files[-grep("ambigous",csv_files)]
anno_files <- lapply(paste0(path, csv_files), read.csv, stringsAsFactors = F)

anno <- rbindlist(anno_files) %>% as.data.frame()

df_ambigous = read.csv(paste0(path,"Lung_28_21-ambigous-pca_20200219.csv"),
                                 stringsAsFactors = F)
df_clean <- readxl::read_excel("doc/Ambiguous-cleaning-2-19-20.xlsx")
df_clean %<>% .[,c("barcodes", "Annotation")] %>% as.data.frame()
colnames(df_clean)[2] = "cell.labels"
#======== adjust cell labels =================
df_unambigous <- anno[!anno$barcodes %in% df_ambigous$barcodes,];dim(df_unambigous)
df_all = rbind(df_unambigous, df_ambigous);dim(df_all)
df_all <- df_all[!df_all$barcodes %in% df_clean$barcodes,];dim(df_all)

df_clean %<>% .[.[,"cell.labels"] != "impurity",];dim(df_clean)
df_all %<>% rbind(df_clean)
table(tuple::triplicated(df_all$barcodes))


df1 = df_all[!duplicated(df_all$barcodes),];dim(df1)
df2 = df_all[duplicated(df_all$barcodes),];dim(df2)
df_duplicated  = inner_join(df1, df2, by ="barcodes");dim(df_duplicated)
df_unique = df_all[!df_all$barcodes %in% df_duplicated$barcodes,];dim(df_unique)

# ambigous have two labels
ambigous <- df_duplicated[df_duplicated[,"cell.labels.x"] != df_duplicated[,"cell.labels.y"],];dim(ambigous)
# unambigous have one label
unambigous <- df_duplicated[df_duplicated[,"cell.labels.x"] == df_duplicated[,"cell.labels.y"],];dim(unambigous)
unambigous = unambigous[,-ncol(unambigous)]
colnames(unambigous)[2] = "cell.labels"
df_unique %<>% rbind(unambigous);dim(df_unique)
table(duplicated(df_unique$barcodes))
rownames(df_unique) = df_unique$barcodes
df_unique %<>% .[.["cell.labels"] != "Un-Prox-6",]
#======== load Seurat =================
object = readRDS(file = "data/Lung_28_Global_20200206.rds") 
table(colnames(object) %in% df_unique$barcodes)
meta.data = object@meta.data
meta.data$barcodes = rownames(meta.data)
meta.data = left_join(meta.data, df_unique, by = "barcodes")
#meta.data[is.na(meta.data$cell.labels),"cell.labels"] = ""
meta.data %<>% cbind(object@reductions$umap@cell.embeddings)
meta.data[(meta.data$UMAP_1 > 1.5 & meta.data$cell.labels %in% "En-SM"),
          "cell.labels"] = "En-V"

rownames(meta.data) = meta.data$barcodes
meta.data %<>% .[colnames(object),]

object[["cell.labels"]] = meta.data$cell.labels
object$cell.labels[is.na(object$cell.labels)] = "zImpurity"
Idents(object) = "cell.labels"
object %<>% sortIdent()
object %<>% AddMetaColor(label= "cell.labels",colors = Singler.colors)
object$cell.labels.colors[is.na(meta.data$cell.labels)] = "#000000"

UMAPPlot.1(object, group.by = "cell.labels",
           label = T, cols = ExtractMetaColor(object),
           label.repel = T, no.legend = T,
           do.print = T, do.return = F)
object %<>% subset(idents = "zImpurity", invert =T)
UMAPPlot.1(object, group.by = "cell.labels",
           label = T, cols = ExtractMetaColor(object),
           label.repel = T, no.legend = T,
           do.print = T, do.return = F)
saveRDS(object, file = "data/Lung_28_Global_20200219.rds")

##################
path <- paste0("Yang/proximal_distal_terminal_COPD/Harmony/table/")
if(!dir.exists(path))dir.create(path, recursive = T)
object = readRDS(file = "data/Lung_28_Global_20200219.rds") 
# - Table: number of cells per cell types (per each sample and total)
df <- table(object$cell.labels, object$orig.ident) %>% 
        as.data.frame()
colnames(df) = c("cell.labels","samples","Freq")

df1 <- spread(df,"samples","Freq")
rownames(df1) = df1$cell.labels
df1 = df1[order(df1$cell.labels),]

write.csv(df1, paste0(path,"Cell_types_by_samples.csv"))

# - Table: number of expressed genes per cluster (per each sample and total) ==============
for(i in seq_len(nrow(df))) {
        cells <- object$cell.labels %in% df[i,"cell.labels"] & object$orig.ident %in% df[i,"samples"]
        df[i,"nGene"] = as.integer(mean(object@meta.data[cells,"nFeature_SCT"]))
        svMisc::progress(i/nrow(df)*100)
}
df2 <- df[,-grep("Freq",colnames(df))]
df2 <- spread(df2, key = "samples", value = "nGene")
df2[is.na(df2)] = 0
rownames(df2) = df2$cell.labels

write.csv(df2,paste0(path,"nGene_by_samples.csv"))

# 4. Expression data for all cells in the S-d, S, SMG-Muc, SMG-Ser, IC clusters per samples (per region P, D, T, COPD; with sample ID) for the following genes: SCGB3A2, SFTPB, SCGB1A1, SFTPC, AGER, PRB4, AZGP1, SERPINB3, SERPINB4, MUC5B, MUC5AC, FOXJ1, KRT5. 
object = readRDS(file = "data/Lung_28_Global_20200219.rds") 
Idents(object) = "cell.labels"
groups <- list(c("S-d", "S", "SMG-Muc", "SMG-Ser", "IC"),
               c("AT1","AT2","AT-p"),
               "MEC",
               c("BC", "BC-p", "Sq"),
               c("Ion", "NEC"),
               "H",
               c("C1", "C2", "C3", "C4"))
unlist(groups) %in% Idents(object) %>% table
for(i in seq_along(groups)){
        print(groups[i])
        sub_object <- subset(object, idents = groups[[i]])
        genes = FilterGenes(sub_object,marker.genes = c("SCGB3A2", "SFTPB", "SCGB1A1", "SFTPC", "AGER", "PRB4", "AZGP1", 
                                                        "SERPINB3", "SERPINB4", "MUC5B", "MUC5AC", "FOXJ1", "KRT5",
                                                        "KRT14", "ACTA2",
                                                        "KRT15","MIR205HG","MKI67","KRT19","BPIFB1","LCN2","AGR2",
                                                        "WFDC2","TMEM45A","NAPSA ","SFTPD","CLDN18","SCEL","AQP4",
                                                        "EMP2","PIP","LTF","SAA1","SAA2","PRB3","TFF3","BPIFB2",
                                                        "LY6D","SPRR2D","SPRR3","SPRR2A","SPRR1A","S100A9","S100A14",
                                                        "KRT6A","SFN","UPK1B","GPX2","ADH7","FOXI1","HEPACAM2","ATP6V0B",
                                                        "TMEM61","CFTR","GRP","SEC11C","CPE","PCSK1N","CHGB","CHGA",
                                                        "CDH19","L1CAM","DCN","SRGN","LAPTM5","CD3E","CDH5","CAPS","TPPP3"))
        exp = sub_object[["SCT"]]@data[genes,]
        meta.data = sub_object@meta.data[,c("orig.ident","conditions","cell.labels")]
        table(rownames(meta.data) == colnames(exp))
        meta.data %<>% cbind(Matrix::t(exp))
        meta.data = meta.data[order(meta.data$cell.labels),]
        write.csv(meta.data,paste0(path,"Expression ",paste(groups[[i]],collapse = ","),".csv"))
        }





Idents(object) = "cell.labels"
sub_object <- subset(object, idents = c("MEC"))
genes = FilterGenes(sub_object,marker.genes = c("KRT14","ACTA2","SCGB3A2", "SFTPB", "SCGB1A1", "SFTPC", "AGER", "PRB4", "AZGP1", 
                                                "SERPINB3", "SERPINB4", "MUC5B", "MUC5AC", "FOXJ1", "KRT5"))
exp = sub_object[["SCT"]]@data[genes,]
meta.data = sub_object@meta.data[,c("orig.ident","conditions","cell.labels")]
table(rownames(meta.data) == colnames(exp))
meta.data %<>% cbind(Matrix::t(exp))
meta.data = meta.data[order(meta.data$cell.labels),]
write.csv(meta.data,paste0(path,"Expression MEC.csv"))
