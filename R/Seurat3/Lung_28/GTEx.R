library(CePa)
library(biomaRt)
library(dplyr)
library(magrittr)
library(Seurat)
library(kableExtra)
source("../R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)
set.seed(101)
# on linux
GTE_meta.data = read.csv("data/RNA-seq/GTEx Portal - lung sample info.csv",
                         stringsAsFactors = F)
rownames(GTE_meta.data) = GTE_meta.data$Tissue.Sample.ID

counts = read.gct("data/RNA-seq/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct")
format(object.size(counts),unit = "GB")

rownames(GTE_meta.data) = gsub("-","\\.", rownames(GTE_meta.data))
select = paste(rownames(GTE_meta.data),collapse = "|")
select = grep(select,colnames(counts), value = T)

counts = counts[,select]
counts = counts[rowSums(counts) > 0,]
write.csv(counts, file="data/RNA-seq/GTEx-Lung")

# on Mac
GTE_meta.data = read.csv("data/RNA-seq/GTEx Portal - lung sample info.csv",
                         stringsAsFactors = F)
rownames(GTE_meta.data) = GTE_meta.data$Tissue.Sample.ID
counts = read.csv("data/RNA-seq/GTEx-Lung-counts.csv",row.names = 1)
df_counts = data.frame(ensembl_gene_id = gsub("\\..*","",rownames(counts)),
                       row.names = rownames(counts))

ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
attributes <- listAttributes(ensembl)
head(attributes,200)
symbols <- getBM(attributes=c("ensembl_gene_id",'hgnc_symbol'), 
                filters = 'ensembl_gene_id', 
                values = df_counts$ensembl_gene_id, 
                mart = ensembl)
df_counts$gene <- symbols[match(df_counts$ensembl_gene_id, symbols$ensembl_gene_id),
                      "hgnc_symbol"]
df_counts %<>% cbind(counts)
df_counts = subset(df_counts, gene != "")
df_counts = df_counts[!duplicated(df_counts$gene),]
rownames(df_counts) = df_counts$gene
df_counts = df_counts[, -grep("ensembl_gene_id|gene",colnames(df_counts))]
df_counts %<>% as.matrix()

colnames(df_counts) %<>% gsub("\\.SM.*","",.)
rownames(GTE_meta.data) %<>% gsub("-",".",.)
GTE_meta.data = GTE_meta.data[colnames(df_counts),]
table(rownames(GTE_meta.data) == colnames(df_counts) )
object <- CreateSeuratObject(counts = df_counts, meta.data = GTE_meta.data,
                             min.cells = 3, min.features = 200)
object %<>% SCTransform
object %<>% FindVariableFeatures(selection.method = "vst",
                               num.bin = 20,nfeatures = 2000,
                               mean.cutoff = c(0.1, 8), 
                               dispersion.cutoff = c(1, Inf))
object %<>%  RunPCA
object %<>%  FindNeighbors(dims = 1:20)
object %<>%  FindClusters
object %<>%  RunUMAP(dims = 1:20)
object@meta.data$Age = gsub("-.*","",object$Age.Bracket) %>% as.integer()
object@meta.data$SCT_snn_res.0.8 %<>% as.numeric()
object@meta.data %>%
        group_by(SCT_snn_res.0.8, Age, Sex) %>%
        summarize(n())
object@meta.data %>%
        group_by(SCT_snn_res.0.8) %>%
        summarize(mean = mean(Age))
UMAPPlot.1(object, group.by = "Age.Bracket", do.print = T, cols = Singler.colors)#Sex Age.Bracket
FeaturePlot.1(object,features = "nCount_RNA")
FeaturePlot.1(object,features = "Age")
save(object, file = "data/Lung_GTEx_20200307.Rda")

df_counts  = as.matrix(object[["RNA"]]@counts)
df_counts = cbind(object[["Age.Bracket"]],t(df_counts))
df_counts = df_counts[order(df_counts[,"Age.Bracket"]),]
write.csv(t(df_counts), file = "data/RNA-seq/GTEx-Lung-counts.csv")
df_counts = read.csv("data/RNA-seq/GTEx-Lung-counts.csv",row.names = 1)

#  ======== DE ============
(load(file = "data/Lung_GTEx_20200307.Rda"))
Idents(object) = "Age.Bracket"
object %<>% sortIdent()
table(Idents(object))
Young_ages <- list(c("20-29","30-39"),
                   c("20-29","30-39"),
                   c("20-29","30-39"),
                   c("20-29"),
                   c("20-29"),
                   c("20-29"),
                   c("20-29"))
old_ages <- list(c("60-69","70-79"),
                 c("50-59","60-69","70-79"),
                 c("40-49","50-59","60-69","70-79"),
                 c("60-69","70-79"),
                 c("50-59","60-69","70-79"),
                 c("40-49","50-59","60-69","70-79"),
                 c("30-39","40-49","50-59","60-69","70-79"))
age_markers <- FindPairMarkers(object, 
                               ident.1 = Young_ages,
                               ident.2 = old_ages,
                               p.adjust.methods = "BH",
                               logfc.threshold = 1,only.pos = F, 
                               min.pct = 0.1)

age_markers$FC = 2^(age_markers$avg_logFC)
write.csv(age_markers,paste0(path,"age_markers_FC0.csv"))
age_markers = age_markers[age_markers$p_val_adj<0.05,]
write.csv(age_markers,paste0(path,"age_markers_FC0-p_val_adj0.05.csv"))

age_markers <- read.csv(paste0("Yang/GTEx/","age_markers_FC0.csv"),row.names = 1)

age_markers %>% group_by()
(top <-  age_markers %>% 
        group_by(cluster1.vs.cluster2) %>% 
        top_n(2, avg_logFC))
grep("SCGB3A2", age_markers$gene)


GCT <- c("BC-BC-p-Sq-IC-output","S-d-SMG-S-morpheus-output","SMG-MEC-morpheus-output")
for(i in seq_along(GCT)){
        counts = read.gct(paste0("data/RNA-seq/",GCT[i],".gct"))
        colnames(counts) = gsub("\\.","-",colnames(counts))
        write.csv(counts,paste0(path,GCT[i],".csv"))
}

gene = "SCGB3A2"
object[[gene]] = as.vector(object@assays$SCT[gene,])
df = cbind(object[["Age.Bracket"]], object[[gene]]) %>% as.data.frame()
df %>%
        group_by(Age.Bracket) %>%
        summarize(SCGB3A2 = mean(SCGB3A2)) %>% kable %>% kable_styling()

SCGB3A2_marker <- list()
for(i in seq_along(Young_ages)){
        SCGB3A2_marker[[i]]<- FindMarkers.UMI(object, ident.1 = Young_ages[[i]],
                                              ident.2 = old_ages[[i]],
                                              features = "SCGB3A2",
                                              p.adjust.methods = "fdr",
                                              return.thresh = 1,
                                              logfc.threshold = 0,only.pos = F, 
                                              min.pct = 0.1)
        ident.1vs2 <- paste(ident.1[i],"/", ident.2[i])
        SCGB3A2_marker[[i]]$cluster1.vs.cluster2 <- ident.1vs2
        
}
SCGB3A2_marker = bind_rows(SCGB3A2_marker)
SCGB3A2_marker$gene = "SCGB3A2"
SCGB3A2_marker$FC = 2^(SCGB3A2_marker$avg_logFC)
write.csv(SCGB3A2_marker,paste0(path,"SCGB3A2_marker.csv"))
