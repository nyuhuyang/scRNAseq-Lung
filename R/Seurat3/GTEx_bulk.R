library(data.table)
library(dplyr)
library(magrittr)
library(xbioc)
library(cowplot)
library(Seurat)
library(Biobase)
library(org.Hs.eg.db)
source("../R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)

TPM_path = "~/Downloads/GTEx_Analysis_2016-01-15_v7_RSEMv1.2.22_transcript_tpm.txt"
TPM = fread(TPM_path)
dim(TPM)
head(TPM[,1:5],)

TPM[, gene_id:= gsub("\\..*","",TPM$gene_id)]
TPM[,transcript_id := NULL]
#remove duplicate rownames with lower rowsumns
#' @param DT input as data.table with gene name
#' @export DT matrix with gene as rownames, no duplicated genes
RemoveDup <- function(DT){
        gene_id <- DT[[1]]
        DT <- DT[,-1]
        DT[is.na(DT)] = 0
        DT[, rowSums := rowSums(DT)]
        DT[, rowNames := 1:nrow(DT)]
        setorder(DT, -rowSums)
        gene_id <- gene_id[as.numeric(DT$rowNames)]
        remove_index <- duplicated(gene_id)
        DT <- DT[!remove_index,]
        DT[,gene_id:= gene_id[!remove_index]]
        DT[,c("rowSums","rowNames") := NULL]
        
        return(DT)
}

# convert to gene symbol
TPM <- RemoveDup(TPM) # takes long time
TPM[,gene_name:= mapIds(org.Hs.eg.db, keys=TPM$gene_id, column="SYMBOL", 
                        keytype="ENSEMBL", multiVals="first")]
TPM <- na.omit(TPM, cols="gene_name")

# annotation ==========================
annotation_path <- "~/Downloads/GTEx_v7_Annotations_SampleAttributesDS.txt"
annotation = fread(annotation_path)
head(annotation[,1:4])
annotation = annotation[annotation$SMTSD %in% "Lung"]
(keep.column = complete.cases(t(annotation)))
annotation = annotation[,keep.column,with=FALSE]
head(annotation)
(samples <- annotation$SAMPID) %>% head

keep.samples <- samples %in% colnames(TPM)
table(keep.samples)
samples <- samples[keep.samples]
TPM <- TPM[, c(samples,"gene_name"), with=FALSE]

TPM <- as.data.frame(TPM)
TPM <- TPM[!duplicated(TPM$gene_name),]
rownames(TPM) = TPM$gene_name
head(TPM[,(ncol(TPM)-1):ncol(TPM)])
TPM = TPM[,-ncol(TPM)]

head(TPM[,1:5],2) %>% t
annotation[keep.samples,1:2] %>% head
annotation = annotation[keep.samples,]
annotation = as.data.frame(annotation)
rownames(annotation) = annotation$SAMPID

# create Seurat object
object <- CreateSeuratObject(counts = TPM, project = "Lung_bulk",
                                meta.data = annotation)
save(object, file = paste0("output/Lung_bulk_",length(samples),"_",
                              gsub("-","",Sys.Date()),".Rda"))

# Visualize QC metrics as a violin plot
jpeg(paste0(path,"S1_nGene_nUMI.jpeg"), units="in", width=10, height=7,res=600)
VlnPlot(object, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
dev.off()
######################################

# run sctransform
object <- SCTransform(object, vars.to.regress = "nCount_RNA", verbose = FALSE)

# Identify the 20 most highly variable genes
top20 <- head(VariableFeatures(object), 20)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(object)
plot2 <- LabelPoints(plot = plot1, points = top20, repel = TRUE)
jpeg(paste0(path,"VariableFeaturePlot_LungBulk.jpeg"), units="in", width=10, height=7,res=600)
print(plot2)
dev.off()
hvf.info <- HVFInfo(object = object)
hvf.info = hvf.info[VariableFeatures(object),]
write.csv(hvf.info, file = paste0(path,"LungBulk_high_variable_genes.csv"))

object %<>% RunPCA(features = VariableFeatures(object),verbose =F,npcs = 100)
object %<>% JackStraw(num.replicate = 20,dims = 100)
object %<>% ScoreJackStraw(dims = 1:85)

jpeg(paste0(path,"JackStrawPlot_LungBulk.jpeg"), units="in", width=10, height=7,res=600)
JackStrawPlot(object, dims = 30:40)+
        ggtitle("JackStrawPlot")+
        theme(text = element_text(size=15),							
              plot.title = element_text(hjust = 0.5,size = 18)) 
dev.off()
# These are now standard steps in the Seurat workflow for visualization and clustering
npcs =38
object %<>% FindNeighbors(reduction = "pca",dims = 1:npcs)
object %<>% FindClusters(reduction = "pca",resolution = 1.2,
                       dims.use = 1:npcs,print.output = FALSE)
object %<>% RunTSNE(reduction = "pca", dims = 1:npcs)
object %<>% RunUMAP(reduction = "pca", dims = 1:npcs)

UMAPPlot.1(object = object, label = F, group.by = "SMNABTCHT", 
           do.return = TRUE, no.legend = T, title = "UMAP plot for all clusters",
           pt.size = 1,label.size = 6, do.print = T)

save(object, file = paste0("output/Lung_bulk_",length(samples),"_", gsub("-","",Sys.Date()),".Rda"))
