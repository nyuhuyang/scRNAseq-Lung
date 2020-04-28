# ######################################################################
invisible(lapply(c("Seurat","dplyr","ggpubr","openxlsx","Hmisc",
                   "magrittr"), function(x) {
                           suppressPackageStartupMessages(library(x,character.only = T))
                   }))
source("../R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)

# SLURM_ARRAY_TASK_ID
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
if (length(slurm_arrayid)!=1)  stop("Exact one argument must be supplied!")
# coerce the value to an integer
args <- as.numeric(slurm_arrayid)
print(paste0("slurm_arrayid=",args))


#======1.2 load  Seurat =========================
object = readRDS(file = "data/Lung_28_Global_20200219.rds") 
table(object$cell.labels)
(cell.labels <- sort(unique(object$cell.labels)))
(label <- cell.labels[args])
sub_object <- subset(object, idents = label)
rm(object);GC()
covid19_genes <- FilterGenes(sub_object,c("ACE2","ANPEP","BSG","CTSL","TMPRSS2"))
## Column clustering (adjust here distance/linkage methods to what you need!)
sub_object <- FindVariableFeatures(object = sub_object, selection.method = "vst",
                               num.bin = 20,nfeatures = 3000,
                               mean.cutoff = c(0.1, 8), 
                               dispersion.cutoff = c(1, Inf))
test_genes <- unique(c(VariableFeatures(sub_object),covid19_genes))
y = sub_object[["SCT"]]@data[test_genes,]
system.time(cor_res <- Hmisc::rcorr(t(as.matrix(y)), type="spearman"))
cor_res$r[is.na(cor_res$r)] = 0
jpeg(paste0(path,"heatmap-cor-",label,".jpeg"), units="in", width=10, height=7,res=600)
heatmap.2(cor_res$r,trace="none")
dev.off()
#write.csv(exp,paste0(path,"exp1.csv"))
df_list <- list()
for(i in seq_along(covid19_genes)){
        gene <- covid19_genes[i]
        df <- data.frame("correlation" = cor_res$r[gene,],
                         "p.value" = cor_res$P[gene,])
        df = df[!is.na(df$correlation),]
        df = df[(rownames(df) != gene),]
        df[,"log10.p.value"] = -log10(df$`p.value`)
        df = df[order(df$correlation),]
        df[,"genes"] = rownames(df)
        n = 10
        low_cor <- head(df$correlation,n)[n]
        high_cor <- tail(df$correlation,n)[1]
        jpeg(paste0(path,"cor-pvalue-",gene,"-",label,".jpeg"), units="in", width=10, height=7,res=600)
        g <- ggline(df, x = "correlation", y = "log10.p.value",
                    numeric.x.axis = TRUE,
                    ylab = "-log10(p-value)",
                    xlab = "Spearman Correlation",
                    label = "genes",             # Add point labels
                    label.select = list(criteria = paste("`x` <=",low_cor,"| `x` >=",high_cor)),           # show only labels for the top 2 points
                    repel = TRUE, 
                    title = paste("Top 3000 variable correlated with",gene,"in",label))+
                theme_minimal() + TitleCenter()
        print(g)
        dev.off()
        colnames(df)[3] = " -log10(p-value)"
        df_list[[i]] = df
        Progress(i,length(covid19_genes))
}
names(df_list) = covid19_genes
write.xlsx(df_list, file = paste0(path,"cor-pvalue-genes-",label,".xlsx"),
           colNames = TRUE, borders = "surrounding",colWidths = c(NA, "auto", "auto"))
