invisible(lapply(c("Seurat","monocle","dplyr","scales",
                   "magrittr"), function(x) {
                           suppressPackageStartupMessages(library(x,character.only = T))
                   }))
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat3_functions.R")
path <- "Yang/Lung_30/Monocle2/"
if(!dir.exists(path)) dir.create(path, recursive = T)
#SBATCH --mem=64G
# SLURM_ARRAY_TASK_ID
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
if (length(slurm_arrayid)!=1)  stop("Exact one argument must be supplied!")
# coerce the value to an integer
i <- as.numeric(slurm_arrayid)
print(paste0("slurm_arrayid=",i))
run_differentialGeneTest = TRUE        
# load data
step = 1
if(step == 1){
        object = readRDS(file = "data/Lung_30_20200710.rds")
        DefaultAssay(object) = "SCT"
        Idents(object) = "Doublets"
        object <- subset(object, idents = "Singlet")
        Idents(object) = "annotations3"
        object <- subset(object, idents = c("BC","BC-S","BC-p","IC1","IC2",
                                            "IC-S","S","S-d","C1","C2","C3",
                                            "NEC","Ion","p-C","H","AT1","AT2",
                                            "AT2-1","AT2-p"))
        object %<>% sortIdent()
        object %<>% AddMetaColor(label= "annotations3", colors = Singler.colors)
        
        Idents(object) = "conditions"
        sample_pairs = list("distal",
                            "terminal",
                            "proximal",
                            "COPD",
                            c("distal","terminal","proximal","COPD"))
        print(sample <- sample_pairs[[i]])
        object %<>% subset(idents = sample)
        GC()
        save.path = paste0(path, paste(sample, collapse = "_"),"/")
        if(!dir.exists(save.path)) dir.create(save.path, recursive = T)
        
        Idents(object) = "annotations3"
        annotations3_color = ExtractMetaColor(object)
        # Store Data in a CellDataSet Object
        pd <- new("AnnotatedDataFrame", data = object@meta.data)
        fd <- data.frame(gene_short_name = rownames(object),
                         row.names = rownames(object))
        fd <- new("AnnotatedDataFrame", data = fd)
        expr_matrix <- as.matrix(object@assays$SCT@counts)
        cds <- newCellDataSet(expr_matrix, phenoData = pd,featureData = fd,
                              expressionFamily = negbinomial())
        #Estimate size factors and dispersions
        cds <- estimateSizeFactors(cds)
        cds <- estimateDispersions(cds)
        table(pData(cds)$cell.types)
        #Filtering low-quality cells
        cds <- detectGenes(cds, min_expr = 0.1)
        print(head(fData(cds)))
        expressed_genes <- row.names(subset(fData(cds),
                                            num_cells_expressed >= 10))
        length(expressed_genes)
        #######################
        #Trajectory step 1: choose genes that define a cell's progress
        if(run_differentialGeneTest){
                clustering_DEG_genes <- differentialGeneTest(cds[expressed_genes,],
                                                             fullModelFormulaStr = "~annotations3",
                                                             cores = detectCores()/2)
                saveRDS(clustering_DEG_genes, paste0(save.path,"monocle2_",paste(sample, collapse = "-"),"_DE.rds"))
                
        } else clustering_DEG_genes = readRDS(paste0(save.path,"monocle2_",paste(sample, collapse = "-"),"_DE.rds"))
        
        print("clustering_DEG_genes")
        cds_ordering_genes <-row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][
                1:min(nrow(clustering_DEG_genes),1000)]
        cds %<>% setOrderingFilter(ordering_genes = cds_ordering_genes)
        cds %<>% reduceDimension(method = 'DDRTree')
        cds %<>% orderCells()
        
        saveRDS(cds, paste0(save.path,"monocle2_",paste(sample, collapse = "-"),"_cds.rds"))
        #cds = readRDS(paste0(save.path,"monocle2_",paste(sample, collapse = "-"),"_cds.rds"))
        #Trajectory step 2: generate Trajectory plot for all samples
        group_by <- c("orig.ident","conditions","group", "Pseudotime","annotations3")
        
        for(k in seq_along(group_by)){
                jpeg(paste0(save.path,"trajectory_",paste(sample, collapse = "-"),"_",group_by[k],".jpeg"),
                     units="in", width=7, height=7,res=600)
                g <- plot_cell_trajectory(cds, color_by = group_by[k],cell_size = 3,
                                          show_branch_points = FALSE)
                if(group_by[k] == "orig.ident") g = g + 
                        scale_color_manual(values=hue_pal()(length(unique(cds$orig.ident))))
                if(group_by[k] == "conditions") g = g + 
                        scale_color_manual(values=hue_pal()(length(unique(cds$conditions))))
                if(group_by[k] == "group") g = g + 
                        scale_color_manual(values=hue_pal()(length(unique(cds$group))))
                if(group_by[k] == "cell.types") g = g + 
                        scale_color_manual(values=sort(unique(cds$annotations3.colors),decreasing = T))
                
                print(g)
                dev.off()
                Progress(k, length(group_by))
        }
}

#Trajectory step 2: generate Trajectory plot by each sample
if(step == 2){
        save.path.sub = paste0(save.path, "subset/")
        if(!dir.exists(save.path.sub)) dir.create(save.path.sub, recursive = T)
        
        cds = readRDS(paste0(save.path,"monocle2_",paste(sample, collapse = "-"),"_cds.rds"))
        samples <- unique(cds$orig.ident)
        group_by <- c("group", "Pseudotime","annotations3")
        
        for(s in seq_along(samples)){
                valid_cells <- row.names(subset(pData(cds), orig.ident == samples[s]))
                subset_cds <- cds[,valid_cells]
                subset_cds@reducedDimS <- subset_cds@reducedDimS[,valid_cells]
                for(k in seq_along(group_by)){
                        jpeg(paste0(save.path.sub,"trajectory_",samples[s],"_",group_by[k],".jpeg"),
                             units="in", width=7, height=7,res=600)
                        g <- plot_cell_trajectory(subset_cds, cell_size = 3,
                                                  color_by = group_by[k],show_branch_points = FALSE)
                        if(group_by[k] == "group") g = g + 
                                scale_color_manual(values=hue_pal()(length(unique(cds$orig.ident))))
                        if(group_by[k] == "annotations3") g = g + 
                                scale_color_manual(values=sort(unique(cds$annotations3.colors),decreasing = T))
                        g = g + ggtitle(paste(group_by[k], "in", samples[s])) + TitleCenter()
                        print(g)
                        dev.off()
                        Progress((s-1)*length(samples)+k, length(samples)*length(group_by))
                }
        }
        
}