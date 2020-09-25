invisible(lapply(c("Seurat","monocle","dplyr","scales",
                   "magrittr"), function(x) {
                           suppressPackageStartupMessages(library(x,character.only = T))
                   }))
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
#SBATCH --mem=128G
# SLURM_ARRAY_TASK_ID
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
if (length(slurm_arrayid)!=1)  stop("Exact one argument must be supplied!")
# coerce the value to an integer
i <- as.numeric(slurm_arrayid)
print(paste0("slurm_arrayid=",i))
Get_DE_genes <- c("run_DE","read_DE","use_variable_genes")[2]

# load data
step = 1
if(step == 1){
        object = readRDS(file = "data/Lung_30_20200710.rds")
        DefaultAssay(object) = "SCT"
        Idents(object) = "Doublets"
        object <- subset(object, idents = "Singlet")
        object %<>% AddMetaColor(label= "annotations3", colors = Singler.colors)
        
        Idents(object) = "conditions"
        sample_pairs = list("distal",#64GB
                            "terminal",#32GB
                            "proximal",#32GB
                            "COPD",#32GB
                            c("distal","terminal","proximal","COPD"))#256GB
        print(sample <- sample_pairs[[i]])
        object %<>% subset(idents = sample)
        
        Idents(object) = "annotations3"
        wanted_cells = c("BC","BC-S","BC-p","IC1","IC2",
                       "IC-S","S","S-d","C1","C2","C3",
                       "NEC","Ion","p-C","H","AT1","AT2",
                       "AT2-1","AT2-p")
        avaible_cells <- as.data.frame(table(Idents(object)))
        avaible_cells = as.character(avaible_cells[avaible_cells$Freq > 10,"Var1"])
        keep_cells = wanted_cells[wanted_cells %in% avaible_cells]
        print(paste("remove cell types:", 
                    paste(wanted_cells[!(wanted_cells %in% avaible_cells)],
                          collapse = ", ")))
        object <- subset(object, idents = keep_cells)
        GC()
        object %<>% sortIdent()
        object %<>% AddMetaColor(label= "annotations3", colors = Singler.colors)
        

        GC()
        save.path = paste0(path,Get_DE_genes,'/', paste(sample, collapse = "_"),"/")
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
        table(pData(cds)$annotations3)

        #Filtering low-quality cells
        cds <- detectGenes(cds, min_expr = 0.1)
        print(head(fData(cds)))
        expressed_genes <- row.names(subset(fData(cds),
                                            num_cells_expressed >= 10))
        length(expressed_genes)
        
        #######################
        #Trajectory step 1: choose genes that define a cell's progress
        DE_genes <- switch(Get_DE_genes,
                           "run_DE" = {clustering_DEG_genes <- differentialGeneTest(cds[expressed_genes,],
                                                                                    fullModelFormulaStr = "~annotations3",
                                                                                    cores = detectCores()/2)
                           saveRDS(clustering_DEG_genes, paste0(save.path,"monocle2_",paste(sample, collapse = "_"),"_DE.rds"))
                           return(row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][
                                   1:min(nrow(clustering_DEG_genes),1000)])
                           },
                           "read_DE" = {clustering_DEG_genes = readRDS(paste0(save.path,"monocle2_",paste(sample, collapse = "_"),"_DE.rds"))
                           return(row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][
                                   1:min(nrow(clustering_DEG_genes),1000)])
                           },
                           "use_variable_genes" = VariableFeatures(object)[1:1000])

        cds %<>% setOrderingFilter(ordering_genes = DE_genes)
        #Do dimensionality reduction
        cds %<>% reduceDimension(norm_method = 'vstExprs', 
                                 reduction_method='DDRTree', 
                                 verbose = F, max_components = 7)
        cds %<>% orderCells()
        
        saveRDS(cds, paste0(save.path,"monocle2_",paste(sample, collapse = "-"),"_cds.rds"))
        #===============
        #Trajectory step 2: generate Trajectory plot for all samples
        #cds = readRDS(paste0(save.path,"monocle2_",paste(sample, collapse = "-"),"_cds.rds"))
        cds$annotations4 = gsub("-.*$","",(cds$annotations3)) %>% gsub("[0-9]$","",.)
        group_by <- c("orig.ident","conditions","group", "Pseudotime","annotations3","annotations4")
        
        for(k in seq_along(group_by)){
                jpeg(paste0(save.path,"trajectory_",paste(sample, collapse = "-"),"_",group_by[k],".jpeg"),
                     units="in", width=7, height=7,res=600)
                g <- plot_cell_trajectory(cds, color_by = group_by[k],cell_size = 1,
                                          show_branch_points = FALSE)
                if(group_by[k] == "orig.ident") g = g + 
                        scale_color_manual(values=hue_pal()(length(unique(cds[[group_by[k]]]))))
                if(group_by[k] == "conditions") g = g + 
                        scale_color_manual(values=hue_pal()(length(unique(cds[[group_by[k]]]))))
                if(group_by[k] == "group") g = g + 
                        scale_color_manual(values=hue_pal()(length(unique(cds[[group_by[k]]]))))
                if(group_by[k] == "annotations3") g = g + 
                        scale_color_manual(values=sort(unique(cds$annotations3.colors),decreasing = T))
                if(group_by[k] == "annotations4") g = g + 
                        scale_color_manual(values=Singler.colors[1:length(unique(cds[[group_by[k]]]))])
                
                print(g)
                dev.off()
        }
        #custom_colour_map = Singler.colors[1:length(unique(pData(cds)$annotations3))]
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