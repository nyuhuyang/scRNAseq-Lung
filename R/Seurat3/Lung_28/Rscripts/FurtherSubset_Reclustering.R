########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################
invisible(lapply(c("Seurat","dplyr","kableExtra","cowplot",
                   "magrittr","MAST","future"), function(x) {
                       suppressPackageStartupMessages(library(x,character.only = T))
                   }))
source("../R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)

# change the current plan to access parallelization
plan("multiprocess", workers = 4)
plan()

# SLURM_ARRAY_TASK_ID
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
if (length(slurm_arrayid)!=1)  stop("Exact one argument must be supplied!")
# coerce the value to an integer
args <- as.numeric(slurm_arrayid)
print(paste0("slurm_arrayid=",args))

groups <- c("AT","B","D","En","F","Mon","SAE","SMG","SMP","T")
(g <- groups[args])
save.path <- paste0(path,g,"/")
if(!dir.exists(save.path))dir.create(save.path, recursive = T)

# ========== cell types and markers ==============
# En (CDH5) 
# SAE:C (CAPS, FOXJ1) 
# F (DCN, LUM)
# SAE (KRT5, SFTPC,SCGB1A1)
# SMG (KRT5, SCGB1A1)
# SAE:Basal cells (KRT5) 
# AT (SFTPC,SCGB1A1+)
# B, T, Mon (SRGN, LAPTM5)
# SAE:Sec (SCGB1A1-very high)
step = 1
if(step == 1){
        markers <- list(c("En","CDH5"),
                c("En","SRGN"),
                c("SAE:C","CAPS"),
                c("SAE:C","FOXJ1"),
                c("F","DCN"),
                c("F","LUM"),
                c("SAE","KRT5"),
                c("SAE","SFTPC"),
                c("SAE","SCGB1A1"),
                c("SMG","KRT5"),
                c("SMG","SCGB1A1"),
                #c("SMP","TAGLN"),
                c("SAE:BC","KRT5"),
                c("AT","SFTPC"),
                c("AT","SCGB1A1"),
                c("B","SRGN"),
                c("B","LAPTM5"),
                c("B","CD19"),
                c("T","SRGN"),
                c("T","LAPTM5"),
                c("T","CD3G"),
                c("Mon","SRGN"),
                c("Mon","LAPTM5"),
                c("SAE:Sec","SCGB1A1"))
        df_markers <- unlist(markers) %>% 
                matrix(nrow=length(markers), byrow=T) %>% 
                as.data.frame()
        colnames(df_markers) = c("cell_types", "markers")
        df_markers$cell_types %<>% gsub(":.*","",.)
        # load data 
        load(file = paste0("data/Lung_28_",g,"_20200122.Rda"))
        DefaultAssay(object) = "SCT"
        object@meta.data %<>% cbind(object@reductions$umap@cell.embeddings)
        # test features
        g_marker_genes <- df_markers[df_markers$cell_types %in% g,"markers"]
        all_marker_genes <- FilterGenes(object,unique(df_markers$markers))
        rm_genes = all_marker_genes[!all_marker_genes %in% g_marker_genes]
        FeaturePlot.1(object, features = all_marker_genes,
                      border = T, do.print = T, do.return = F,
                      unique.name = "groups",  save.path = save.path)
        object %<>% AddModuleScore(features = list(rm_genes),
                                   name = paste0(g,"_impurity"))
        fig <- FeaturePlot.1(object, features = paste0(g,"_impurity1"), 
                             title = paste(g,"impurity"), do.return = T)
        if(g == "AT") {
                fig = fig +
                geom_vline(xintercept = -6)+
                geom_vline(xintercept = 5)+
                geom_hline(yintercept = -7)+
                geom_hline(yintercept = 7)
                object %<>% subset(UMAP_1 > -6 & UMAP_1 < 5 &
                                   UMAP_2 > -7 & UMAP_2 < 7)
        }
        if(g == "B") {
                fig = fig +
                        geom_vline(xintercept = 1.5)+
                        geom_vline(xintercept = 6)+
                        geom_hline(yintercept = 4.5)
                object %<>% subset(UMAP_1 > 1.5 & UMAP_2 > 4.5, invert = T)
                object %<>% subset(UMAP_2 < 6)
        }
        if(g == "En") {
                fig = fig +
                        geom_vline(xintercept = -2)+
                        geom_vline(xintercept = 2.7)+
                        geom_vline(xintercept = 5)+
                        geom_hline(yintercept = 2.6)+
                        geom_hline(yintercept = 4.7)
                object %<>% subset(UMAP_1 > -2 & UMAP_2 > 4.7, invert = T)
                object %<>% subset(UMAP_1 > 2.7 & UMAP_2 > 2.6, invert = T)
                object %<>% subset(UMAP_1 < 5)
        }
        if(g == "F") {
                fig = fig +
                        geom_vline(xintercept = 5)+
                        geom_hline(yintercept = 6.5)+
                        geom_hline(yintercept = 1)
                object %<>% subset(UMAP_1 > 5 & UMAP_2 < 1, invert = T)
                object %<>% subset(UMAP_1 < 6.5)
        }
        if(g == "F") {
                fig = fig +
                        geom_vline(xintercept = 5)+
                        geom_hline(yintercept = 6.5)+
                        geom_hline(yintercept = 1)
                object %<>% subset(UMAP_1 > 5 & UMAP_2 < 1, invert = T)
                object %<>% subset(UMAP_1 < 6.5)
        }
        if(g == "Mon") {
                fig = fig +
                        geom_vline(xintercept = -9)+
                        geom_vline(xintercept = -2)+
                        geom_vline(xintercept = 4)+
                        geom_hline(yintercept = 6)+
                        geom_hline(yintercept = 1)
                object %<>% subset(UMAP_1 > -2 & UMAP_2 > 6, invert = T)
                object %<>% subset(UMAP_1 > 4 & UMAP_2 > 1, invert = T)
                object %<>% subset(UMAP_1 > -9)
        }
        if(g == "SAE") {
                fig = fig +
                        geom_vline(xintercept = -5)+
                        geom_vline(xintercept = 0)+
                        geom_hline(yintercept = -1)+
                        geom_hline(yintercept = 0)+
                        geom_hline(yintercept = 2.2)
                object %<>% subset(UMAP_1 < -5 & UMAP_2 < -1, invert = T)
                object %<>% subset(UMAP_1 > 0 & UMAP_2 > 0 & UMAP_2 < 2.2,
                                   invert = T)
        }
        if(g == "SMG") {
                fig = fig +
                        geom_vline(xintercept = -0.8)+
                        geom_vline(xintercept = 0.5)+
                        geom_hline(yintercept = 2.2)
                object %<>% subset(UMAP_1 > -0.8 & UMAP_1 < 0.5 &UMAP_2 > 2.2, invert = T)
        }
        if(g == "SMP") {
                fig = fig +
                        geom_vline(xintercept = -1.5)+
                        geom_vline(xintercept = 0)+
                        geom_hline(yintercept = -2.5)+
                        geom_hline(yintercept = -6)+
                        geom_hline(yintercept = -7)
                object %<>% subset(UMAP_1 > -1.5 & UMAP_2 > -6 & UMAP_2 < -2.5, invert = T)
                object %<>% subset(UMAP_1 > 0 & UMAP_2 > -7 & UMAP_2 < -2.5, invert = T)
        }
        if(g == "T") {
                fig = fig +
                        geom_vline(xintercept = -7)+
                        geom_vline(xintercept = -5)+
                        geom_vline(xintercept = -4)+
                        geom_hline(yintercept = 3.2)+
                        geom_hline(yintercept = 1.5)+
                        geom_hline(yintercept = 1)+
                        geom_hline(yintercept = -1)
                object %<>% subset(UMAP_1 < -5 & UMAP_2 > 3.2, invert = T)
                object %<>% subset(UMAP_1 < -7 & UMAP_2 > 1.5, invert = T)
                object %<>% subset(UMAP_1 > -7 & UMAP_1 < -5 & UMAP_2 < -1, invert = T)
                object %<>% subset(UMAP_1 > -5 & UMAP_1 < -4 & UMAP_2 < 1, invert = T)
        }
        jpeg(paste0(save.path,"UMAPPlot_",g,".jpeg"), 
             units="in", width=10, height=7,res=600)
        print(fig)
        dev.off()
        FeaturePlot.1(object, features = paste0(g,"_impurity1"),
                      title = paste(g,"impurity"), do.return = F,
                      do.print = T, unique.name = "groups",  
                      save.path = save.path)
        save(object, file = paste0("data/Lung_28_",g,"_20200126.Rda"))
}