####################################
invisible(lapply(c("Seurat","dplyr","magrittr","tidyr","SCINA"), function(x) {
        suppressPackageStartupMessages(library(x,character.only = T))
}))
source("../R/Seurat3_functions.R")
source("../R/SingleR_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)

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
                                N = 10000, max_iter = 100, convergence_n = 10, 
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

lapply(c(T,F), function(x) {
        UMAPPlot.1(object, group.by = "SCINA",cols = ExtractMetaColor(object),label = x,
                   label.repel = T, pt.size = 0.5,label.size = 3, repel = T,no.legend = x,
                   unique.name = "conditions",
                   do.print = T,do.return = F,title = "Cell types in all 28 samples")}
)
