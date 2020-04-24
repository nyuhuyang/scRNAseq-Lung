# ######################################################################
invisible(lapply(c("Seurat","dplyr","cowplot","openxlsx",
                   "magrittr"), function(x) {
                           suppressPackageStartupMessages(library(x,character.only = T))
                   }))
source("../R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)
#======1.2 load  Seurat =========================
object = readRDS(file = "data/Lung_28_Global_20200219.rds") 

#=======3.2 gene co-expression summary ================
cell_labels <- c("All_cells",sort(unique(object$cell.labels)))
Idents(object) = "cell.labels"
group_by = "orig.ident"
g = c("ACE2","TMPRSS2")
df_list <- list()
for(i in seq_along(cell_labels)){
        label = cell_labels[i]
        if(label == "All_cells") {
                single_object <- object
        } else single_object <- subset(object, idents = label)
        meta = single_object[[group_by]]
        meta[, g[1]] = as.vector(single_object@assays$SCT[g[1],]) >0 
        meta[, g[2]] = as.vector(single_object@assays$SCT[g[2],]) >0
        meta[, paste0(g[1],"_",g[2])] = meta[, g[1]] & meta[, g[2]]
        meta[, g[1]] %<>% as.integer()
        meta[, g[2]] %<>% as.integer()
        meta[, paste0(g[1],"_",g[2])] %<>% as.integer()
        df = aggregate(x = meta[,-1], by = list(meta[,group_by]), FUN = sum)
        df$sample = table(meta[,group_by]) %>% as.vector()

        fisher <- function(vector, g = g){
                output <- vector()
                #if(vector[g[1]] == 0 | vector[g[2]] == 0) return(c(1,1))
                vector = c(vector["sample"], vector[c(g[1], g[2],paste0(g[1],"_",g[2]))])
                vector %<>% unlist %>% matrix(nrow = 2)
                vector[1,2] = vector[1,2] - vector[2,2]
                vector[2,1] = vector[2,1] - vector[2,2]
                vector[1,1] = vector[1,1] - vector[2,1] - vector[2,1]
                FISH <- fisher.test(vector,workspace = 2000000)
                
                output["p_value"] = FISH$p.value
                output["p_val_adj"] = p.adjust(p = FISH$p.value, method = "BH", 
                                        n = nrow(df))
                return(output)
        }
        df_fisher = apply(df[, -1],1,fisher, g = g)
        df %<>% cbind(t(df_fisher))
        df_percentage = aggregate(x = meta[,-1], by = list(meta[,group_by]), 
                                  FUN = function(x) sum(x >0) / length(x) *100)
        colnames(df_percentage) %<>% paste0("_%")
        df %<>% cbind(df_percentage[,-1])
        
        df = df[,c("Group.1","sample", g[1],paste0(g[1],"_%"),g[2],paste0(g[2],"_%"),
                   paste0(g[1],"_",g[2]),paste0(g[1],"_",g[2],"_%"),"p_value","p_val_adj")]
        
        Idents(single_object) = "orig.ident"
        exp = AverageExpression(single_object, assays = "SCT", features = c("ACE2","TMPRSS2","ST14","FURIN"))
        exp = exp$SCT %>% t()
        
        df %<>% cbind(exp[df[,1],])
        colnames(df) = c("Sample ID","Total number of cells (A)",
                         "Number of ACE+ cells", "% ACE2+ cells (% of A)",
                         "Number of TMPRSS+ cells","% TMPRSS2+ cells (% of A)",
                         "Number of double+ cells","% double+ cells (% of A)",
                         "p_value","p_val_adj",
                         "ACE2","TMPRSS2","ST14","FURIN")
        df_list[[label]] = df
        Progress(i, length(cell_labels))
}
write.xlsx(df_list, file = "Yang/Covid19/CoV target gene summary.xlsx",
           colNames = TRUE, borders = "surrounding",colWidths = c(NA, "auto", "auto"))



#=======3.2 dotplot================
markers.to.plot <- c("ACE2","TMPRSS2","ST14","FURIN","ANPEP","BSG","CTSL","MUC16")
Idents(object) = "orig.ident"
Idents(object) %<>% factor(levels = rev(levels(object)))
jpeg(paste0(path,"dotplot_sample.jpeg"), units="in", width=6, height=7,res=600)
DotPlot(object, features = rev(markers.to.plot), assay = "SCT",
        cols = c("blue","red")) + RotatedAxis()
dev.off()

Idents(object) = "cell.labels"
Idents(object) %<>% factor(levels = rev(levels(object)))
jpeg(paste0(path,"dotplot_celltypes.jpeg"), units="in", width=6, height=10,res=600)
DotPlot(object, features = rev(markers.to.plot), assay = "SCT",
        cols = c("blue","red")) + RotatedAxis()
dev.off()

#=======3.2 dotplot================
(load("data/Lung_GTEx_20200307.Rda"))
markers.to.plot <- c("ACE2","TMPRSS2","ST14","FURIN","ANPEP","BSG","CTSL","MUC16")
Idents(object) = "Age"
object %<>% sortIdent()
Idents(object) %<>% factor(levels = rev(levels(object)))
jpeg(paste0(path,"dotplot_age.jpeg"), units="in", width=5, height=7,res=600)
DotPlot(object, features = rev(markers.to.plot), assay = "SCT",
        split.by = "Sex",
        cols = c("blue","red")) + RotatedAxis()
dev.off()

write.csv(table(object$Sex, object$Age), file = paste0(path, "bulk_sex.csv"))


object[["Age_Sex"]] = paste0(object$Age, "_", object$Sex)
Idents(object) = "Age_Sex"
exp = AverageExpression(object, assays = "SCT", 
                        features = c("ACE2","TMPRSS2","ST14","FURIN","ANPEP","BSG","CTSL","MUC16")
)
exp = exp$SCT
write.csv(t(exp), file = paste0(path, "cov_target_bulk_expression.csv"))

#=======3.3 scatterplot================
object = readRDS(file = "data/Lung_28_Global_20200219.rds") 
cell_labels <- sort(unique(object$cell.labels))
Idents(object) = "cell.labels"
group_by = "orig.ident"
genes = c("ACE2","TMPRSS2")
for(i in seq_along(cell_labels)){
        label = cell_labels[i]
        if(label == "All_cells") {
                single_object <- object
        } else single_object <- subset(object, idents = label)
        jpeg(paste0(path,"ScatterPlot_",paste(genes,collapse = "_"),"_",label,".jpeg"), 
             units="in", width=7, height=7,res=600)
        g <- FeatureScatter(single_object, 
                            feature1 = genes[1],
                            feature2 = genes[2],
                            pt.size = 2,#cols = "black",
                            slot = "data")+ NoLegend()
        print(g)
        dev.off()
        Progress(i, length(cell_labels))
}

#=======3.2 gene co-expression summary ================
Idents(object) = "cell.labels"
Alveolar_epithelium <- subset(object, idents = c("AT1", "AT2"))
write.csv(Alveolar_epithelium[["SCT"]]@data, file = paste0(path,"Alveolar_epithelium.csv"))

Airway_epithelium <- subset(object, idents = c("BC", "BC-p", "IC", "Sq", "S", "S-d", "C1",
                                               "C2", "C3", "C4", "p-C", "H", "Ion", "NEC",
                                               "SMG-Muc", "SMG-Ser", "MEC"))
write.csv(Airway_epithelium[["SCT"]]@data, file = paste0(path,"Airway_epithelium.csv"))

epithelium <- subset(object, idents = c("AT1", "AT2","BC", "BC-p", "IC", "Sq", "S", "S-d", "C1",
                                        "C2", "C3", "C4", "p-C", "H", "Ion", "NEC",
                                        "SMG-Muc", "SMG-Ser", "MEC"))
write.csv(epithelium[["SCT"]]@data, file = paste0(path,"epithelium.csv"))

Macrophages <- subset(object, idents = c("M0", "M1", "M2"))
write.csv(Macrophages[["SCT"]]@data, file = paste0(path,"Macrophages.csv"))

cell_labels <- sort(as.character(unique(Idents(object))))
for(i in seq_along(cell_labels)){
        sub_object <- subset(object, idents = cell_labels[i])
        write.csv(sub_object[["SCT"]]@data, file = paste0(path,cell_labels[i],".csv"))
        Progress(i, length(cell_labels))
}

