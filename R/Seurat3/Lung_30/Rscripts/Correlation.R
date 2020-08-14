# ######################################################################
invisible(lapply(c("Seurat","dplyr","ggpubr","openxlsx","Hmisc",
                   "magrittr"), function(x) {
                       suppressPackageStartupMessages(library(x,character.only = T))
                   }))
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)

# SLURM_ARRAY_TASK_ID
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
if (length(slurm_arrayid)!=1)  stop("Exact one argument must be supplied!")
# coerce the value to an integer
args <- as.integer(as.character(slurm_arrayid))
print(paste0("slurm_arrayid=",args))

#======1.2 load  Seurat =========================
object = readRDS(file = "data/Lung_30_20200702.rds") 
object %<>% sortIdent()
table(Idents(object))

types = c("spearman", "pearson")
type = types[args]

## Column clustering (adjust here distance/linkage methods to what you need!)
y = object[["SCT"]]@data[VariableFeatures(object),]
system.time(cor_res <- Hmisc::rcorr(t(as.matrix(y)), type=type))
cor_res$r[is.na(cor_res$r)] = 0
jpeg(paste0(path,"heatmap-cor-",type,".jpeg"), units="in", width=10, height=7,res=600)
heatmap.2(cor_res$r,trace="none")
dev.off()
hm <- heatmap.2(cor_res$r,trace="none")
write.csv(cor_res$r[rev(hm$rowInd), hm$colInd], 
          file = paste0(path,"Lung_30_correlation_",type,".csv"))
saveRDS(cor_res, file = paste0("data/Lung_30_cor_",type,".rds"))