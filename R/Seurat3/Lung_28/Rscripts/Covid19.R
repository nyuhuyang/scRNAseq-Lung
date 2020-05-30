# ######################################################################
invisible(lapply(c("Seurat","dplyr","cowplot","openxlsx",
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
args <- as.numeric(slurm_arrayid)
print(paste0("slurm_arrayid=",args))

#======1.2 load  Seurat =========================
object = readRDS(file = "data/Lung_28_Global_20200219.rds") 

#=======3.2 gene co-expression summary ================
cell_labels <- c("All_cells",sort(unique(object$cell.labels)))
Idents(object) = "cell.labels"
group_by = "orig.ident"
gene_pairs <- list(c("ACE2","TMPRSS2"),
                   c("ACE2" , "CTSL"),
                   c("ANPEP" , "TMPRSS2"),
                   c("ANPEP" , "CTSL"),
                   c("BSG" , "TMPRSS2"),
                   c("BSG" , "CTSL"),
                   c("MUC16" , "TMPRSS2"),
                   c("MUC16" , "CTSL"))
(g <- gene_pairs[[args]])
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
                vector = c(vector["sample"], vector[c(g[1],g[2],paste0(g[1],"_",g[2]))])
                m  <- matrix(unlist(vector), nrow = 2)
                m[1,2] = m[1,2] - m[2,2]
                m[2,1] = m[2,1] - m[2,2]
                m[1,1] = m[1,1] - m[2,1] - m[1,2] - m[2,2]
                FISH <- fisher.test(m,workspace = 2000000)
                output <- vector()
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
write.xlsx(df_list, file = paste0(path,"CoV target genes- ",paste(g, collapse = " and "),".xlsx"),
           colNames = TRUE, borders = "surrounding",colWidths = c(NA, "auto", "auto"))
