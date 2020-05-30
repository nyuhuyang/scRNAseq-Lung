####################################
invisible(lapply(c("Seurat","dplyr","magrittr","tidyr","SCINA"), function(x) {
        suppressPackageStartupMessages(library(x,character.only = T))
}))
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat3_functions.R")
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/SingleR_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)

# SLURM_ARRAY_TASK_ID
set.seed=101
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
if (length(slurm_arrayid)!=1)  stop("Exact one argument must be supplied!")
# coerce the value to an integer
args <- as.numeric(slurm_arrayid)
print(paste0("slurm_arrayid=",args))

# load Seurat object
(load(file = "data/Lung_28_20200103.Rda"))
# load DEGs and prepare signatures
DEGs <- read.csv("Yang/proximal_distal_terminal/Non-Integration/DEGs/cell_types/Lung_24-FC0_cell_types.csv",
                 stringsAsFactors = F)
DEGs$X = make.unique(DEGs$gene)
DEGs = DEGs[DEGs$cluster != "Unknown",]
top <-  DEGs %>% 
        group_by(cluster) %>% 
        top_n(20, avg_logFC) %>%
        .[,c("X","cluster","gene")]
df_signatures <- top %>% spread(cluster, gene)
df_signatures =  df_signatures[,-grep("X",colnames(df_signatures))]
signatures <- df2list(df_signatures)
sapply(signatures,length)

print(system.time(results <- BigSCINA(exp = object@assays$SCT@data, signatures,
                                N = args, max_iter = 100, convergence_n = 10, 
                                convergence_rate = 0.999, sensitivity_cutoff = 0.9, rm_overlap=FALSE, allow_unknown=FALSE, 
                                log_file=paste0(path,'SCINA.log'))))
table(results$cell_labels)
save(results, file = paste0("output/SCINA_",gsub("-","",Sys.Date()),".Rda"))

# Prepare UMAP
object@meta.data$SCINA = results$cell_labels
df_cell_types <- readxl::read_excel("doc/Cell type abbreviation.xlsx")
object$SCINA %<>% plyr::mapvalues(
        from = df_cell_types$`Cell types`,
        to = df_cell_types$Abbreviation)
Idents(object) = "Doublets"
object %<>% subset(idents = "Singlet")
Idents(object) = "SCINA"
object %<>% sortIdent()
object %<>% AddMetaColor(label= "SCINA", colors = c(Singler.colors,Singler.colors))
UMAPPlot.1(object, group.by = "SCINA",cols = ExtractMetaColor(object),label = T,
           label.repel = T, pt.size = 0.5,label.size = 3, repel = T,no.legend = T,
           unique.name = "conditions",
           do.print = T,do.return = F,
           title = "Cell types in all 28 samples")
TSNEPlot.1(object, group.by = "SCINA",cols = ExtractMetaColor(object),label = T,
           label.repel = T, pt.size = 0.5,label.size = 3, repel = T,no.legend = T,
           unique.name = "conditions",
           do.print = T,do.return = F,title = "Cell types in all 28 samples")

Idents(object) = "conditions"
conditions <- c("proximal", "distal", "terminal", "COPD")
Idents(object) = "conditions"
for(i in seq_along(conditions)){
        sub_object <- subset(object, idents = conditions[i])
        Idents(sub_object) = "SCINA"
        
        UMAPPlot.1(sub_object, group.by = "SCINA",cols = ExtractMetaColor(sub_object),label = T,
                   label.repel = T, pt.size = 0.5,label.size = 3, repel = T,no.legend = T,
                   unique.name = "conditions",
                   do.print = T,do.return = F,
                   title = paste("Cell types in",conditions[i]))
        TSNEPlot.1(sub_object, group.by = "SCINA",cols = ExtractMetaColor(sub_object),label = T,
                   label.repel = T, pt.size = 0.5,label.size = 3, repel = T,no.legend = T,
                   unique.name = "conditions",
                   do.print = T,do.return = F,
                   title = paste("Cell types in",conditions[i]))
        Progress(i,length(conditions))
}

