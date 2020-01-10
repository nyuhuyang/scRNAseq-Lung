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

# Predict cell types with SCINA
# COPD only
Idents(object) = "conditions"
COPD <- subset(object, idents = "COPD")
system.time(results_COPD <- SCINA(exp = as.matrix(COPD@assays$SCT@data), signatures, max_iter = 100, convergence_n = 10, 
                convergence_rate = 0.999, sensitivity_cutoff = 0.9, rm_overlap=FALSE, allow_unknown=FALSE,
                log_file='SCINA.log'))
results_COPD$cell_labels %<>% gsub("unknown","Unknown",.)
table(results_COPD$cell_labels)
table(COPD@meta.data$cell.types,results_COPD$cell_labels )
COPD@meta.data$SCINA_allow_unknown = results_COPD$cell_labels
df_cell_types <- readxl::read_excel("doc/Cell type abbreviation.xlsx")
COPD$SCINA_allow_unknown %<>% plyr::mapvalues(
        from = df_cell_types$`Cell types`,
        to = df_cell_types$Abbreviation)
Idents(COPD) = "Doublets"
COPD %<>% subset(idents = "Singlet")
Idents(COPD) = "SCINA"
COPD %<>% sortIdent()
COPD %<>% AddMetaColor(label= "SCINA", colors = c(Singler.colors,Singler.colors))


lapply(c(T,F), function(x) {
        UMAPPlot.1(COPD, group.by = "SCINA_allow_unknown",cols = ExtractMetaColor(COPD),label = x,
                   label.repel = T, pt.size = 0.5,label.size = 3, repel = T,no.legend = x,
                   unique.name = "conditions",
                   do.print = T,do.return = F,title = "Cell types in COPD samples")}
)
save(results_COPD, file = paste0("output/SCINA_COPD_",gsub("-","",Sys.Date()),".Rda"))



system.time(results <- BigSCINA(exp = object@assays$SCT@data, signatures,
                                N = 5000, max_iter = 100, convergence_n = 10, 
                convergence_rate = 0.999, sensitivity_cutoff = 0.9, rm_overlap=FALSE, allow_unknown=FALSE, 
                log_file=paste0(path,'SCINA.log')))
table(results$cell_labels)

write.csv(object@reductions$pca@feature.loadings, paste0(path,"PCA_matrix.csv"))
mt <- object@reductions$pca@feature.loadings
feature.loadings = data.frame(gene = rownames(mt),
                              PC = sapply(1:nrow(mt),function(x) which.max(mt[x,])))
feature.loadings = feature.loadings[order(feature.loadings$PC),]
write.csv(feature.loadings, paste0(path,"feature.loadings.csv"))
