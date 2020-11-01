# conda activate r4.0
invisible(lapply(c("Seurat","monocle","dplyr","scales","plyr","stringr",
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

#Trajectory step 1: generate Trajectory plot within regions
step = 1
if(step == 1){
        sample_string = c("distal",#128GB
                           "terminal",#32GB
                           "proximal",#32GB
                           "COPD",#32GB
                           "distal,terminal,proximal,COPD")#256GB
        opts = data.frame(methods = rep(c("UseVariableGenes","ReadDE"),each = 10),
                          samples = rep(sample_string, 4),
                          root = rep(c(NA,"BC",NA,"BC"), each = 5),
                          stringsAsFactors = F)
        print(args <- opts[i,])
        
        Get_DE_genes <- args$methods
        sample <- str_split(args$samples, pattern = ",")[[1]]
        root = args$root
        
        object = readRDS(file = "data/Lung_30_20200710.rds")
        DefaultAssay(object) = "SCT"
        Idents(object) = "Doublets"
        object <- subset(object, idents = "Singlet")
        object$annotations4 <- mapvalues(object$annotations3,
                                         from = c("BC",
                                                  "BC-S",
                                                  "IC1",
                                                  "IC2",
                                                  "BC-p",
                                                  "H",
                                                  "p-C",
                                                  "C1",
                                                  "C2",
                                                  "C3",
                                                  "AT1",
                                                  "AT2",
                                                  "AT2-1",
                                                  "AT2-p"),
                                           to = c("BC",
                                                  "BC",
                                                  "IC",
                                                  "IC",
                                                  "IC",
                                                  "H-p-C",
                                                  "H-p-C",
                                                  "C",
                                                  "C",
                                                  "C",
                                                  "AT",
                                                  "AT",
                                                  "AT",
                                                  "AT"),
                                           warn_missing = TRUE)
        object %<>% AddMetaColor(label= "annotations3", colors = Singler.colors)
        object %<>% AddMetaColor(label= "annotations4", colors = Singler.colors)
        
        Idents(object) = "conditions"
        object %<>% subset(idents = sample)
        
        # remove small group
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
        object = subset(object, idents = keep_cells)
        object = FindVariableFeatures(object, selection.method = "vst", nfeatures = 1000)
        GC()

        save.path = paste0(path,i,"-",Get_DE_genes,'-', paste(sample, collapse = "."),"-root=",root,"/")
        if(!dir.exists(save.path)) dir.create(save.path, recursive = T)
        
        # Store Data in a CellDataSet Object
        pd <- new("AnnotatedDataFrame", data = object@meta.data)
        fd <- data.frame(gene_short_name = rownames(object),
                         row.names = rownames(object))
        fd <- new("AnnotatedDataFrame", data = fd)
        cds <- newCellDataSet(object@assays$SCT@counts,
                              phenoData = pd,
                              featureData = fd,
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
        print(length(expressed_genes))
        
        #######################
        #Trajectory step 1: choose genes that define a cell's progress
        DE_genes <- switch(Get_DE_genes,
                           "run_DE" = {clustering_DEG_genes <- differentialGeneTest(cds[expressed_genes,],
                                                                                    fullModelFormulaStr = "~annotations3",
                                                                                    cores = detectCores()/2)
                           saveRDS(clustering_DEG_genes, paste0(save.path,basename(save.path),"_cds.rds"))
                           row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][
                                   1:min(nrow(clustering_DEG_genes),1000)]
                           },
                           "ReadDE" = {
                                   clustering_DEG_genes = read.csv("Yang/Lung_30/DE_analysis/C_Cell_types/top1000_epi_genes.csv")
                           as.vector(clustering_DEG_genes[,1])
                           },
                           "UseVariableGenes" = VariableFeatures(object)[1:1000]
                           )

        cds %<>% setOrderingFilter(ordering_genes = DE_genes)
        #Do dimensionality reduction
        cds %<>% reduceDimension(norm_method = 'vstExprs', 
                                 reduction_method='DDRTree', 
                                 verbose = F, max_components = 7)
        cds %<>% orderCells()
        
        if(identical(root,"BC")) {
                #Find the correct root state the corresponds to the 'BC' 
                tab1 <- table(pData(cds)$State, pData(cds)$annotations3)
                id = which(colnames(tab1) == 'BC')
                root_name = names(which.max(tab1[,id]))
                
                #Run a second time to get the correct root state that overlaps with Stem cells
                cds <- orderCells(cds, root_state=root_name)
        } else root_name = NULL
        
        saveRDS(cds, paste0(save.path,basename(save.path),"_cds.rds"))
        #===============
        #Trajectory step 2: generate Trajectory plot for all samples
        #cds = readRDS(paste0(save.path,basename(save.path),"_cds.rds"))
        group_by <- c("orig.ident","group", "Pseudotime","annotations3","annotations4")
        
        for(k in seq_along(group_by)){
                g <- plot_cell_trajectory(cds, color_by = group_by[k],cell_size = 1,
                                          show_branch_points = FALSE)
                
                jpeg(paste0(save.path,"trajectory_",paste(sample, collapse = "-"),"_",group_by[k],".jpeg"),
                     units="in", width=7, height=7,res=600)
                if(group_by[k] == "orig.ident") g = g + 
                        scale_color_manual(values=hue_pal()(length(unique(cds[[group_by[k]]]))))
                if(group_by[k] == "group") g = g + 
                        scale_color_manual(values=hue_pal()(length(unique(cds[[group_by[k]]]))))
                if(group_by[k] == "annotations3") g = g + 
                        scale_color_manual(values=sort(unique(cds$annotations3.colors),decreasing = T))
                if(group_by[k] == "annotations4") g = g + 
                        scale_color_manual(values=sort(unique(cds$annotations4.colors),decreasing = T))
                
                print(g)
                dev.off()
        }
        #Trajectory step 3: generate complex_cell_trajectory
        #Get a nice colour map

        for(k in seq_along(group_by)){
                g <-  plot_complex_cell_trajectory(cds, color_by = group_by[k], show_branch_points = T, 
                                                   cell_size = 2, cell_link_size = 1, root_states = c(root_name)) +
                        scale_size(range = c(0.2, 0.2)) +
                        theme(legend.position="right", legend.title=element_blank(), legend.text=element_text(size=rel(1.5))) +
                        guides(colour = guide_legend(override.aes = list(size=6)))
                
                jpeg(paste0(save.path,"complex_cell_trajectory_",paste(sample, collapse = "-"),"_",group_by[k],".jpeg"),
                     units="in", width=7, height=7,res=600)

                if(group_by[k] == "orig.ident") g = g + 
                        scale_color_manual(values=hue_pal()(length(unique(cds[[group_by[k]]]))))
                if(group_by[k] == "group") g = g + 
                        scale_color_manual(values=hue_pal()(length(unique(cds[[group_by[k]]]))))
                if(group_by[k] == "annotations3") g = g + 
                        scale_color_manual(values=sort(unique(cds$annotations3.colors),decreasing = T))
                if(group_by[k] == "annotations4") g = g + 
                        scale_color_manual(values=sort(unique(cds$annotations4.colors),decreasing = T))
                print(g)
                dev.off()
        }
        lib_info_with_pseudo <- pData(cds)
        data_df <- t(monocle::reducedDimS(cds)) %>% as.data.frame() %>% 
                select_(dim_1 = 1, dim_2 = 2) %>% 
                tibble::rownames_to_column("barcodes") %>% 
                left_join(lib_info_with_pseudo %>% tibble::rownames_to_column("barcodes"), 
                          by = "barcodes")
        write.csv(data_df, file = paste0(save.path,basename(save.path),"_coordinates.csv"))
        
}

#Trajectory step 2: generate Trajectory plot across regions
if(step == 2){
        object = readRDS(file = "data/Lung_30_20200710.rds")
        DefaultAssay(object) = "SCT"
        Idents(object) = "Doublets"
        object <- subset(object, idents = "Singlet")
        object$annotations4 <- mapvalues(object$annotations3,
                                         from = c("BC",
                                                  "BC-S",
                                                  "IC1",
                                                  "IC2",
                                                  "BC-p",
                                                  "H",
                                                  "p-C",
                                                  "C1",
                                                  "C2",
                                                  "C3",
                                                  "AT1",
                                                  "AT2",
                                                  "AT2-1",
                                                  "AT2-p"),
                                         to = c("BC",
                                                "BC",
                                                "IC",
                                                "IC",
                                                "IC",
                                                "H-p-C",
                                                "H-p-C",
                                                "C",
                                                "C",
                                                "C",
                                                "AT",
                                                "AT",
                                                "AT",
                                                "AT"),
                                         warn_missing = TRUE)
        object %<>% AddMetaColor(label= "annotations3", colors = Singler.colors)
        object %<>% AddMetaColor(label= "annotations4", colors = Singler.colors)
        Idents(object) = "annotations3"
        wanted_cells = c("BC","BC-S","BC-p","IC1","IC2",
                         "IC-S","S","S-d","C1","C2","C3",
                         "NEC","Ion","p-C","H","AT1","AT2",
                         "AT2-1","AT2-p")
        object %<>% subset(idents = wanted_cells)
        
        Idents(object) = "conditions"
        sample_list = list(c("distal","terminal","proximal"),#128GB
                            c("distal","terminal","proximal","COPD"))#128GB
        print(sample <- sample_list[[i]])
        object %<>% subset(idents = sample)

        save.path = paste0(path,'Acorss_regions/', paste(sample, collapse = "_"),"/")
        if(!dir.exists(save.path)) dir.create(save.path, recursive = T)

        object$`annotations4_conditions` = paste0(object$annotations4," (", object$conditions,")")
        object %<>% AddMetaColor(label= "annotations4_conditions", colors = Singler.colors)
        
        # Store Data in a CellDataSet Object
        pd <- new("AnnotatedDataFrame", data = object@meta.data)
        fd <- data.frame(gene_short_name = rownames(object),
                         row.names = rownames(object))
        fd <- new("AnnotatedDataFrame", data = fd)
        cds <- newCellDataSet(object@assays$SCT@counts,
                              phenoData = pd,
                              featureData = fd,
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
        clustering_DEG_genes = read.csv("Yang/Lung_30/DE_analysis/A_Sample_types/top1000_epi_genes.csv")
        DE_genes <- as.vector(clustering_DEG_genes[,1])

        cds %<>% setOrderingFilter(ordering_genes = DE_genes)
        #Do dimensionality reduction
        cds %<>% reduceDimension(norm_method = 'vstExprs', 
                                 reduction_method='DDRTree', 
                                 verbose = F, max_components = 7)
        cds %<>% orderCells()
        #Find the correct root state the corresponds to the 'BC' 
        tab1 <- table(pData(cds)$State, pData(cds)$annotations4)
        id = which(colnames(tab1) == 'BC')
        root_name = names(which.max(tab1[,id]))
        
        #Run a second time to get the correct root state that overlaps with Stem cells
        cds <- orderCells(cds, root_state=root_name)
        
        saveRDS(cds, paste0(save.path,"monocle2_",paste(sample, collapse = "-"),"_cds.rds"))
        #===============
        #Trajectory step 2: generate Trajectory plot for all samples
        #cds = readRDS(paste0(save.path,"monocle2_",paste(sample, collapse = "-"),"_cds.rds"))
        group_by <- c("conditions", "Pseudotime","annotations3","annotations4",
                      "annotations4_conditions")
        
        for(k in seq_along(group_by)){
                g <- plot_cell_trajectory(cds, color_by = group_by[k],cell_size = 1,
                                          show_branch_points = FALSE)
                
                jpeg(paste0(save.path,"trajectory_",paste(sample, collapse = "-"),"_",group_by[k],".jpeg"),
                     units="in", width=7, height=7,res=600)
                if(group_by[k] == "orig.ident") g = g + 
                        scale_color_manual(values=hue_pal()(length(unique(cds[[group_by[k]]]))))
                if(group_by[k] == "group") g = g + 
                        scale_color_manual(values=hue_pal()(length(unique(cds[[group_by[k]]]))))
                if(group_by[k] == "annotations3") g = g + 
                        scale_color_manual(values=sort(unique(cds$annotations3.colors),decreasing = T))
                if(group_by[k] == "annotations4") g = g + 
                        scale_color_manual(values=sort(unique(cds$annotations4.colors),decreasing = T))
                if(group_by[k] == "annotations4_conditions") g = g + 
                        scale_color_manual(values=sort(unique(cds$annotations4_conditions.colors),decreasing = T))
                
                print(g)
                dev.off()
        }
        #Trajectory step 3: generate complex_cell_trajectory
        #Get a nice colour map
        
        for(k in seq_along(group_by)){
                g <-  plot_complex_cell_trajectory(cds, color_by = group_by[k], show_branch_points = T, 
                                                   cell_size = 2, cell_link_size = 1, root_states = c(root_name)) +
                        scale_size(range = c(0.2, 0.2)) +
                        theme(legend.position="right", legend.title=element_blank(), legend.text=element_text(size=rel(1.5))) +
                        guides(colour = guide_legend(override.aes = list(size=6)))
                
                jpeg(paste0(save.path,"complex_cell_trajectory_",paste(sample, collapse = "-"),"_",group_by[k],".jpeg"),
                     units="in", width=7, height=7,res=600)
                
                if(group_by[k] == "orig.ident") g = g + 
                        scale_color_manual(values=hue_pal()(length(unique(cds[[group_by[k]]]))))
                if(group_by[k] == "group") g = g + 
                        scale_color_manual(values=hue_pal()(length(unique(cds[[group_by[k]]]))))
                if(group_by[k] == "annotations3") g = g + 
                        scale_color_manual(values=sort(unique(cds$annotations3.colors),decreasing = T))
                if(group_by[k] == "annotations4") g = g + 
                        scale_color_manual(values=sort(unique(cds$annotations4.colors),decreasing = T))
                if(group_by[k] == "annotations4_conditions") g = g + 
                        scale_color_manual(values=sort(unique(cds$annotations4_conditions.colors),decreasing = T))
                
                print(g)
                dev.off()
        }
}