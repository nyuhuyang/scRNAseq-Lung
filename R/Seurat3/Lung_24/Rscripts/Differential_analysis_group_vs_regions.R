########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################
invisible(lapply(c("Seurat","dplyr","tidyr","magrittr","gplots","MAST"), function(x) {
        suppressPackageStartupMessages(library(x,character.only = T))
}))
source("../R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/Pairwise_SCT/")
if(!dir.exists(path))dir.create(path, recursive = T)

# SLURM_ARRAY_TASK_ID
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
if (length(slurm_arrayid)!=1)  stop("Exact one argument must be supplied!")
# coerce the value to an integer
args <- as.numeric(slurm_arrayid)
print(paste0("slurm_arrayid=",args))

Pairs = c("P.vs.D","P.vs.T","P.vs.(D+T)","D.vs.T","D.vs.(P+T)","T.vs.(P+D)")
labels = c("all.cells","cell.types","group1","group2","group3","group4","group5")
labels_Pairs=paste0(rep(labels, each=6),"_", Pairs)
(p <- labels_Pairs[args])

# load data
(load(file = "data/Lung_24_20191206.Rda"))
DefaultAssay(object) = "SCT"
Idents(object) = "Doublets"
object %<>% subset(idents = "Singlet")
# group cell types
object@meta.data$group1 = gsub("cells:Artery","cells-Artery",object@meta.data$cell.types)
object@meta.data$group1 = gsub("cells:Capillary","cells-Capillary",object@meta.data$group1)
object@meta.data$group1 = gsub("Endothelial cells:.*","Endothelial cells-Others",object@meta.data$group1)

object@meta.data$group1 = gsub(":.*","",object@meta.data$group1)

object$group1[grep("Secretory cells|Mucus-producing cells",object$group1)] = "Surface secretory cells"
object$group1[grep("NK/T cytotoxic cells",object$group1)] = "T cells"

object$group2 = object$group1
object$group2[grep("Basal.*|Ciliated.*|Secretory.*|Intermediate.*|Squamous|Pre-ciliated cells|Hybrid|Neuroendocrine.*|Ionocytes",
                          object$group2)] = "Surface Airway epithelial cells"
object$group2[grep("Smooth.*|Pericytes",object$group2)] = "Smooth muscle cells & Pericytes" 
object$group2[grep("Endothelial cells-.*",object$group2)] = "Endothelial cells" 
object$group2[grep("Macrophages|Dendritic cells|T cells|B cells|Monocytes|Neutrophils|Mast cells|Plasma cells",
                   object$group2)] = "Immune cells"

object$group3 = object$group2
object$group3[grep("Surface Airway epithelial cells|Myoepithelial cells|Submucosal gland",
                   object$group3)] = "Airway epithelial cells"
object$group3[grep("Fibroblasts|Smooth muscle cells & Pericytes|Stromal cells",
                   object$group3)] = "Stromal cells"

object$group4 = object$group3
object$group4[grep("Surface secretory cells|Submucosal gland",object$group4)] = "All secretory cells" 
object$group4[grep("Stromal cells|Chondrocytes",object$group4)] = "Stromal cells & Chondrocytes" 
object$group4[grep("Airway epithelial cells|Alveolar type.*",object$group4)] = "All epithelial cells" 

object$group5 = object$group4
object$group5[grep("Stromal cells & Chondrocytes|Endothelial cells|Immune cells",object$group5)] = "All non-epithelial cells" 

# subset
label = sub("\\_.*","",p)
pair = sub(".*\\_","",p)

Idents(object) = "conditions"
if(pair == "P.vs.D"){
        sub_object <- subset(object, idents = c("proximal","distal"))
}
if(pair == "P.vs.T"){
        sub_object <- subset(object, idents = c("proximal","terminal"))
}
if(pair == "P.vs.(D+T)"){
        sub_object <- object
        sub_object@meta.data$conditions %<>% gsub("distal|terminal","distal+terminal",.)
}
if(pair == "D.vs.T"){
        sub_object <- subset(object, idents = c("distal","terminal"))
}
if(pair == "D.vs.(P+T)"){
        sub_object <- object
        sub_object@meta.data$conditions %<>% gsub("proximal|terminal","proximal+terminal",.)
}
if(pair == "T.vs.(P+D)"){
        sub_object <- object
        sub_object@meta.data$conditions %<>% gsub("proximal|distal","proximal+distal",.)
}

if(label == "group1"){
        Idents(sub_object) = "group1"
        sub_object <- subset(sub_object, idents = c("Basal cells",
                                                    "Ciliated cells",
                                                    "Dendritic cells",
                                                    "Fibroblasts",
                                                    "Smooth muscle cells",
                                                    "Endothelial cells-Artery",
                                                    "Endothelial cells-Capillary",
                                                    "Macrophages",
                                                    "Surface secretory cells"))
}
if(label == "group2"){
        Idents(sub_object) = "group2"
        sub_object <- subset(sub_object, idents = c("Surface Airway epithelial cells",
                                                    "Smooth muscle cells & Pericytes",
                                                    "Endothelial cells",
                                                    "Immune cells"))
}
if(label == "group3"){
        Idents(sub_object) = "group3"
        sub_object <- subset(sub_object, idents = c("Airway epithelial cells",
                                                    "Stromal cells"))
}
if(label == "group4"){
        Idents(sub_object) = "group4"
        sub_object <- subset(sub_object, idents = c("All epithelial cells",
                                                    "All secretory cells",
                                                    "Stromal cells & Chondrocytes"))
}
if(label == "group5"){
        Idents(sub_object) = "group5"
        sub_object <- subset(sub_object, idents = c("All non-epithelial cells"))
}

if(label == "all.cells"){
        Idents(sub_object) = "conditions"
        Lung_markers <- FindAllMarkers.UMI(sub_object,
                                           logfc.threshold = 0.5, only.pos = F,
                                           test.use = "MAST")
}

if(label != "all.cells"){
        cell.type = unique(sub_object@meta.data[,label]) %>% sort
        sub_object@meta.data[,paste0(label,"_conditions")] = paste0(sub_object@meta.data[,label],"_",
                                                    as.character(sub_object@meta.data$conditions))
        Idents(sub_object) = paste0(label,"_conditions")
        
        # remove cluster with less than 3 cells======
        label_conditions <- as.factor(sub_object@meta.data[,paste0(label,"_conditions")])
        levels = paste0(cell.type,"_",rep(unique(sub_object$conditions),
                                          each =length(cell.type)))
        label_conditions %<>% factor(levels = levels)
        table_subset <- table(label_conditions) %>% as.data.frame
        (remove <- table_subset[table_subset$Freq < 3,"label_conditions"] %>% 
                        gsub("\\_.*","",.) %>% unique() %>% as.character())
        (cell.type <- cell.type[!(cell.type %in% remove)])
        Idents(sub_object) = label
        sub_object %<>% subset(idents= remove, invert = T)
        
        # Differential analysis
        DefaultAssay(object) = "SCT"
        Idents(sub_object) = paste0(label,"_conditions")
        sub_object %<>% sortIdent
        
        ident1 =switch(sub(".*vs\\.","",pair),
                       "P" = "_proximal",
                       "D" = "_distal",
                       "T" = "_terminal",
                       "(D+T)" = "_distal+terminal",
                       "(P+T)" = "_proximal+terminal",
                       "(P+D)" = "_proximal+distal")
        ident2 = switch(sub("\\..*","",pair), 
                        "P" = "_proximal",
                        "D" = "_distal", 
                        "T" = "_terminal")

        Lung_markers <- FindPairMarkers(sub_object,
                                        ident.1 = paste0(cell.type, ident1), 
                                        ident.2 = paste0(cell.type, ident2),
                                        logfc.threshold = 0.05, only.pos = F,
                                        return.thresh = 0.05,
                                        test.use = "MAST",
                                        save.files = FALSE)
}

write.csv(Lung_markers,paste0(path,"Lung_24_SCT_",p,"~.csv"))
Lung_markers = Lung_markers[Lung_markers$p_val_adj<0.05,]
write.csv(Lung_markers,paste0(path,"Lung_24_SCT_",p,".csv"))


#========= label correction switch cluster1 vs cluster 2=========
#DEGs_path <- "Yang/proximal_distal_terminal/Non-Integration/DEGs/Pairwise/"
#csv_files <- list.files(DEGs_path,pattern = "Lung_24") 
#length(csv_files)
#labels_pairs <- gsub(".csv","",csv_files) %>% gsub("Lung_24_","",.)
#which(!(labels_Pairs %in% labels_pairs))
#csv_files = csv_files[7:length(csv_files)]

#for (i in seq_along(csv_files)) {
#        deg <- read.csv(paste0(DEGs_path, csv_files[i]),row.names = 1)
#        temp <-  stringr::str_split(deg$cluster1.vs.cluster2," vs. ")
#        temp <- sapply(temp, function(x) paste(x[2],"/",x[1]))
#        deg$cluster1.vs.cluster2 = temp
#        write.csv(deg, file = paste0(DEGs_path, csv_files[i]), )
#        Progress(i,length(csv_files))
#        }
