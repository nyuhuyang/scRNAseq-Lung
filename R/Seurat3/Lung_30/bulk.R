invisible(lapply(c("Seurat","dplyr","magrittr","tibble","edgeR","Biobase"), function(x) {
                           suppressPackageStartupMessages(library(x,character.only = T))
                   }))
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat3_functions.R")
source("R/Seurat3/differential_expression.R")

path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)

FPKM <- readxl::read_excel("data/RNA-seq/P-D paired samples FPKM revised.xlsx", sheet = "FPKM")
meta.data <- readxl::read_excel("data/RNA-seq/P-D paired samples FPKM revised.xlsx", sheet = "meta.data")
meta.data %<>% column_to_rownames(var = "...1") %>% t %>% as.data.frame()
meta.data$patient %<>% as.factor()
meta.data$region %<>% as.factor()

# remove dup gene by mean expression
FPKM = FPKM[order(rowMeans(FPKM[,-1]),decreasing = TRUE),]
FPKM = FPKM[!duplicated(FPKM$gene),]
FPKM %<>% column_to_rownames(var = "gene") 
TPM = pbapply::pbapply(FPKM,2, function(x) x/sum(x)* 10^6) 
head(TPM)
table(rownames(meta.data) == colnames(TPM) )
write.csv(TPM,file = "data/RNA-seq/P-D paired samples TPM revised.csv")

##========= remove Pseudogene genes and others=============
rm_genes <- grep("^AC[0-9+]|^AL[0-9+]|^LINC[0-9+]",rownames(TPM),value = T)
rm_genes %<>% c(grep("^CTD-|^CTB-|^CTC-",rownames(TPM),value = T))
rm_genes %<>% c(grep("^RP[0-9+]|^RPL|^RPS",rownames(TPM),value = T))
rm_genes %<>% c(grep("^MT-|^MIR",rownames(TPM),value = T))
rm_genes %<>% sort
length(rm_genes)
write.csv(rm_genes,file = "Yang/RNA-seq/rm_genes.csv",quote = F,
          row.names = F)
rm_genes1 = read.csv(file = "data/RNA-seq/gene_remove_Renat.csv")
rm_genes = unique(c(rm_genes, rm_genes1$x))
rm = which(rownames(TPM) %in% rm_genes)
length(rm)
##========= remove low expression genes =============
drop <- which(apply(TPM, 1, function(x) sum(as.numeric(x > 1))) <= 3)
drop = unique(c(drop, rm))
length(drop)

TPM <- TPM[-drop,]
FPKM <- FPKM[-drop,] 
dim(TPM);dim(FPKM)

log2FC_table <- function(data){
        df <- cbind(data, log2(data+1)) %>% as.data.frame()
        n = ncol(data)/2
        for(i in 1:nrow(data)) {
                df[i,"ttest"] = t.test(x = log2(data[i,1:n]+1), y = log2(data[i,(n+1):(2*n)]+1),paired = T)$p.value
                df[i,"wilcox"] = wilcox.test(x = log2(data[i,1:n]+1), y = log2(data[i,(n+1):(2*n)]+1),paired = T)$p.value
                svMisc::progress(i/nrow(data)*100)
        }
        df$ttest.adj = p.adjust(df$ttest, method = "BH")
        df$wilcox.adj = p.adjust(df$wilcox, method = "BH")
        
        FC_list <- list()
        for(i in 1:n) {
                FC_list[[i]] =  log2(TPM[,i+6]+1) - log2(TPM[,i]+1)
        }
        names(FC_list) = paste0("log2FC_D_vs_P_", unique(meta.data$patient))
        FC <- bind_cols(FC_list)
        FC$Mean_log2FC_D_vs_P  = rowMeans(FC)
        df %<>% cbind(FC)
        df = df[order(df$ttest,decreasing = F),]
        return(df)
}

df1 = log2FC_table(TPM)
write.csv(df1,file = "Yang/RNA-seq/TPM_log2FC.csv",row.names = T)
FPKM %<>% as.matrix()

df2 = log2FC_table(FPKM)
write.csv(df2,file = "Yang/RNA-seq/FPKM_log2FC.csv",row.names = T)


FPKM_FC <- readxl::read_excel("~/Downloads/D-P log2FC paired t test RS 3-9-21.xlsx")
renat_deg <- FPKM_FC$`remove non-coding`
which(renat_deg %in% c("SFTPB","DGCR11","PSKH1","COX7A1"))

library(eulerr)
green_tpm = rownames(df1)[df1$ttest <0.05 & df1$Mean_log2FC_D_vs_P >0]
blue_tpm = rownames(df1)[df1$ttest <0.05 & df1$Mean_log2FC_D_vs_P <0]
geneNames <- list("green" = renat_deg[5:881],
                  "blue" = renat_deg[882:1962],
                  "green_tpm" = green_tpm,
                  "blue_tpm" = blue_tpm
                  )
plot(euler(geneNames, shape = "circle"), quantities = TRUE)

# 2. Determine intersection of genes up-regulated in distal (green highlight in column BC) and proximal (blue highlight in column BC) with cell type DEG marker genes in airway surface epithelium in the single-cell data (you did recently after p-C revision). For overlapped genes, please keep the columns for log2 fold-change and paired t test p value (for bulk RNA-seq data) and log2 fold-change and adj p value for single-cell marker genes.

sheets = c("BC1","BC2","BC-p","IC1","IC2","IC3","S","d-S",
  "H","p-C","C1","C2","C3","Ion","NE")
file_path = "Yang/Lung_30/DE_analysis/groups/DE_results_Surface Airway Epithelial.xlsx"
DEG_list <- lapply(sheets, function(x) {
        tmp = readxl::read_excel(file_path,sheet = x)
        tmp = tmp[tmp$p_val_adj <0.05,]
        return(tmp) 
        }
        )
names(DEG_list) = sheets
DEG = bind_rows(DEG_list)

geneNames <- list("bulk_RNA" = rownames(TPM),
                  "Surface_Airway_Epi" = unique(DEG$gene)
)

plot(euler(geneNames[c("bulk_RNA","Surface_Airway_Epi")], shape = "circle"), quantities = TRUE)

#=== limma =====
# https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/thursday/DE.html
exp0 <- ExpressionSet(assayData=TPM,
                     phenoData = AnnotatedDataFrame(data=meta.data))
cutoff <- 1
drop <- which(apply(cpm(exp0), 1, max) < cutoff)
exp <- exp0[-drop,] 
dim(exp) # number of genes left

plotMDS(exp, col = as.numeric(exp$patient))

patient = exp$patient
region = exp$region
design <- model.matrix(~0+patient+region)
head(design)

y <- voom(exp, design, plot = T)
fit <- lmFit(y, design)
head(coef(fit))

contr <- makeContrasts(~ region)
contr

tmp <- contrasts.fit(fit, contrasts=contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)
# ======== wilcox =========
object <- CreateSeuratObject(counts = TPM, meta.data = as.data.frame(meta.data))
min(object[["RNA"]]@counts)
object[["RNA"]]@data = as(log(as.matrix(object[["RNA"]]@counts+1)), "sparseMatrix")
range(Matrix::rowSums(object[["RNA"]]@counts))
VlnPlot(object = object, features = c("nFeature_RNA", "nCount_RNA"))

object %<>% FindVariableFeatures(selection.method = "vst",
                                 num.bin = 20,nfeatures = 2000,
                                 mean.cutoff = c(0.1, 8), 
                                 dispersion.cutoff = c(1, Inf))
object %<>%  ScaleData
dims = ncol(FPKM)-1
object %<>%  RunPCA(npcs=dims)
#object %<>%  FindNeighbors(dims = 1:dims)
#object %<>%  FindClusters
#object %<>%  RunUMAP(dims = 1:dims)
PCAPlot.1(object, group.by = "region", do.print = T)
PCAPlot.1(object, group.by = "patient", do.print = T, cols = Singler.colors)

FeaturePlot.1(object,features = "nCount_RNA")
save(object, file = "data/Lung_bulk_20210226.Rda")
load("data/Lung_bulk_20210226.Rda")
# differential anlaysis
Idents(object) = "region"
sub_object <- subset(object, idents = c("Prox","Distal"))
marker <- FindAllMarkers.UMI(sub_object,
                            test.use = "wilcox",
                            paired = TRUE,
                            p.adjust.methods = "BH",
                            logfc.threshold = 0.1,
                            return.thresh = 0.05,
                            only.pos = T, verbose = T,
                            min.pct = 0)
write.csv(marker,paste0(path,"P_vs_D_wilcox.csv"))

object$region1 = gsub("Distal|Terminal","Distal,Terminal",object$region)
Idents(object) = "region1"
marker <- FindAllMarkers.UMI(object,
                             test.use = "wilcox",
                             paired = F,
                             p.adjust.methods = "BH",
                             logfc.threshold = 0.1,
                             return.thresh = 0.01,
                             only.pos = T, 
                             min.pct = 0)
write.csv(marker,paste0(path,"P_vs_D+T_paired.csv"))
