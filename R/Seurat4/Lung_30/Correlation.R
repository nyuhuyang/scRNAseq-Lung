# ######################################################################
invisible(lapply(c("Seurat","dplyr","ggpubr","openxlsx","Hmisc","gplots","tibble",
                   "magrittr"), function(x) {
                           suppressPackageStartupMessages(library(x,character.only = T))
                   }))
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)

set.seed(101)
#======1.2 load  Seurat =========================
# load data
object = readRDS(file = "data/Lung_SCT_30_20210831.rds")
DefaultAssay(object) = "SCT"

meta.data = readRDS("output/20211004/meta.data_Cell_subtype.rds")
table(rownames(object@meta.data) == rownames(meta.data))
table(object$barcode ==meta.data$barcode)
object@meta.data = meta.data

object %<>% subset(subset = Doublets %in% "Singlet")
object %<>% subset(subset = Cell_subtype != "Un")

Idents(object)= "Cell_subtype"
object %<>% sortIdent()
table(Idents(object))

Cell_subtype.colors =ExtractMetaColor(object)
# dendrogram based on top 2000 variable genes ==================
object_exp <- AverageExpression(object,assays = "SCT", features = VariableFeatures(object))
exp = object_exp$SCT %>% t %>% scale %>% t

opts = data.frame(hclust_methods = rep(c("complete", "average","ward.D2"),2),
                  cor_methods = rep(c("spearman","pearson"), each = 3),
                  stringsAsFactors = F)
exp_list <- list()
Top_n = 2000
for(i in 1:nrow(opts)){
        hclust_method = opts[i,"hclust_methods"]
        cor_method = opts[i,"cor_methods"]
        
        print(paste(cor_method, "correlation with",hclust_method,"linkage"))
        hc <- hclust(as.dist(1-cor(exp, method=cor_method)), method=hclust_method)
        cc = Cell_subtype.colors[hc$labels]
        
        hm <- heatmap.2(as.matrix(exp),
                  Colv = as.dendrogram(hc),
                  Rowv= TRUE,
                  dendrogram = "column")
        exp_list[[i]] = exp[rev(hm$rowInd), hm$colInd] # revese rowInD to match heatmap.2 output
        
        jpeg(paste0(path,"Heatmap2_",cor_method,"_",hclust_method,".jpeg"), units="in", width=10, height=6.5,res=600)
        heatmap.2(as.matrix(exp),
                  Colv = as.dendrogram(hc),
                  Rowv= TRUE,
                  ColSideColors = cc, 
                  trace ="none",
                  dendrogram = "column",
                  key.xlab = "scale log nUMI",
                  cexRow = 0.5,
                  margins = c(12,5),
                  breaks = seq(-3,3,length.out = 101),
                  col = bluered,
                  main = paste("Clustering dendrogram for all cell types\n based on top",Top_n, 
                               "genes\n by",cor_method, "correlation with",hclust_method,"linkage"))
        dev.off()
}
names(exp_list) = paste(opts$cor_methods, opts$hclust_methods)
write.xlsx(exp_list, file = paste0(path,"Lung_30_correlation_clustering.xlsx"),
           colNames = TRUE, rowNames = TRUE,borders = "surrounding",colWidths = c(NA, "auto", "auto"))

# - Expression data for all 2000 variable genes for each of the cell type 
# (average UMI per cell type; using the same order of cell types as in Spearman complete)
template <- readxl::read_excel("doc/75 genes for dendrogram dot plot.xlsx")
cell.types <- colnames(template) %>% .[-1]
x = as.matrix(exp)
Rowv <- rowMeans(x, na.rm = TRUE)
distr <- dist(x)
hcr <- hclust(distr)
ddr <- as.dendrogram(hcr)
reorderfun = function(d, w) reorder(d, w)
ddr <- reorderfun(ddr, Rowv)
rowInd <- order.dendrogram(ddr)

write.xlsx(object_exp$SCT[rowInd,cell.types], file = paste0(path,"Lung_30_UMI_2000_variable_genes.xlsx"),
           colNames = TRUE, rowNames = TRUE,borders = "surrounding",colWidths = c(NA, "auto", "auto"))

# dendrogram based on top 75 variable genes ==================
#================== DE on group - C - Cell types ================
read.path = "Yang/Lung_30/DE_analysis/"
list_files_C <- list.files(path = paste0(read.path,"C_Cell_types"),
                           pattern ="Lung_30-",full.names = T)
deg_list <- list()
for(i in seq_along(list_files_C)) {
        df <- read.csv(list_files_C[i],header = TRUE, row.names = 1)
        if(nrow(df) > 75) {
                df1 <- head(df, 75)
                df2 <- df[df$avg_logFC >= 2,]
                df <- rbind(df1,df2)
                df <- df[!duplicated(df$gene),]
        }
        deg_list[[i]] = df
        Progress(i, length(list_files_C))
}
deg <- bind_rows(deg_list)
deg_genes <- unique(deg$gene);length(deg_genes)
table(Idents(object))
object_exp <- AverageExpression(object,assays = "SCT", features = deg_genes)
exp = object_exp$SCT %>% t %>% scale %>% t

opts = data.frame(hclust_methods = rep(c("complete", "average","ward.D2"),2),
                  cor_methods = rep(c("spearman","pearson"), each = 3),
                  stringsAsFactors = F)
exp_list <- list()
Top_n = 75
for(i in 1:nrow(opts)){
        hclust_method = opts[i,"hclust_methods"]
        cor_method = opts[i,"cor_methods"]
        
        print(paste(cor_method, "correlation with",hclust_method,"linkage"))
        hc <- hclust(as.dist(1-cor(exp, method=cor_method)), method=hclust_method)
        cc = as.character(as.numeric(as.factor(hc$labels)))
        
        hm <- heatmap.2(as.matrix(exp),
                        Colv = as.dendrogram(hc),
                        Rowv= TRUE,
                        dendrogram = "column")
        exp_list[[i]] = exp[rev(hm$rowInd), hm$colInd] # revese rowInD to match heatmap.2 output
        
        jpeg(paste0(path,"Heatmap2_",cor_method,"_",hclust_method,".jpeg"), units="in", width=10, height=6.5,res=600)
        heatmap.2(as.matrix(exp),
                  Colv = as.dendrogram(hc),
                  Rowv= TRUE,
                  ColSideColors = cc, 
                  trace ="none",
                  dendrogram = "column",
                  key.xlab = "scale log nUMI",
                  cexRow = 0.5,
                  margins = c(12,5),
                  breaks = seq(-3,3,length.out = 101),
                  col = bluered,
                  main = paste("Clustering dendrogram\n based on top",Top_n, 
                               "cell types markers\n by",cor_method, "correlation with",hclust_method,"linkage"))
        dev.off()
}
names(exp_list) = paste(opts$cor_methods, opts$hclust_methods)
write.xlsx(exp_list, file = paste0(path,"Lung_30_correlation_clustering.xlsx"),
           colNames = TRUE, rowNames = TRUE,borders = "surrounding",colWidths = c(NA, "auto", "auto"))
# ======= specify the order =====================

# ==============================
## Column clustering (adjust here distance/linkage methods to what you need!)
y = object[["SCT"]]@data[VariableFeatures(object),]
system.time(cor_res <- Hmisc::rcorr(t(as.matrix(y)), type="spearman"))
cor_res$r[is.na(cor_res$r)] = 0
jpeg(paste0(path,"heatmap-cor-spearman.jpeg"), units="in", width=10, height=7,res=600)
heatmap.2(cor_res$r,trace="none")
dev.off()
hm <- heatmap.2(cor_res$r,trace="none")
write.csv(cor_res$r[rev(hm$rowInd), hm$colInd], file = paste0(path,"Lung_30_correlation_spearman.csv"))
saveRDS(cor_res, file = "data/Lung_30_cor_spearman.rds")


y = object[["SCT"]]@data[VariableFeatures(object),]
system.time(cor_res <- Hmisc::rcorr(t(as.matrix(y)), type="pearson"))
cor_res$r[is.na(cor_res$r)] = 0
jpeg(paste0(path,"heatmap-cor-pearson.jpeg"), units="in", width=10, height=7,res=600)
heatmap.2(cor_res$r,trace="none")
dev.off()
hm <- heatmap.2(cor_res$r,trace="none")
write.csv(cor_res$r[rev(hm$rowInd), hm$colInd], file = paste0(path,"Lung_30_correlation_pearson.csv"))

saveRDS(cor_res, file = "data/Lung_30_cor_pearson.rds")

Methods <- c("complete", "average","ward.D2")
for(method in Methods){
        print(method)
        object <- BuildClusterTree(object, method = method)
        
        jpeg(paste0(path,"PlotClusterTree_Lung30_",method,".jpeg"), units="in", width=14, height=8.5,res=600)
        PlotClusterTree(object)
        title(paste("PlotClusterTree using",method,"algorithm"))
        dev.off()
}
