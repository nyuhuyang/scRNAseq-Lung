#Identification of important genes
list.of.packages <- c("devtools","dplyr","pheatmap","VGAM", "irlba",
                      "matrixStats", "igraph", "combinat", "fastICA",
                      "grid", "reshape2", "plyr", "parallel", "methods")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
#check package
library("devtools")
library("DDRTree")
library("pheatmap")
library("M3Drop")
library("monocle")
library("reshape")
source("../R/Seurat_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
# 5.1 Importing data from Seurat object=================
(load(file="data/epithelial_Harmony_5_20190123.Rda"))
table(epithelial@ident)
cds <- importCDS(epithelial, import_all = TRUE)

# 5.1.1 Estimate size factors and dispersions
# estimateSizeFactors() and estimateDispersions() will only work,
# and are only needed, if you are working with a CellDataSet 
# with a negbinomial() or negbinomial.size() expression family.
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

# 5.1.2 Filtering low-quality cells
cds <- detectGenes(cds, min_expr = 0.1)
print(head(fData(cds)))
print(head(pData(cds)))

# 5.1.3 If you are using RPC values to measure expresion, 
# as we are in this vignette, it's also good to look at the distribution
# of mRNA totals across the cells:
pData(cds)$Total_mRNAs <- Matrix::colSums(exprs(cds))
upper_bound <- 10^(mean(log10(pData(cds)$Total_mRNAs)) + 
                     2*sd(log10(pData(cds)$Total_mRNAs)))
lower_bound <- 10^(mean(log10(pData(cds)$Total_mRNAs)) - 
                     2*sd(log10(pData(cds)$Total_mRNAs)))
qplot(Total_mRNAs, data = pData(cds), color = conditions, geom = "density") +
  geom_vline(xintercept = lower_bound) +
  geom_vline(xintercept = upper_bound)

# 5.2 Classifying and counting cells of different types
# 5.2.1 Classifying cells with manualHierarchy
epithelial <- SetAllIdent(epithelial, id = "manual")
all(row.names(pData(cds)) == names(epithelial@ident))
pData(cds)$manual <- epithelial@ident
table(pData(cds)$manual)
pie <- ggplot(pData(cds), aes(x = factor(1), fill = factor(manual))) +
  geom_bar(width = 1)
pie + coord_polar(theta = "y") +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())

# 5.3 Constructing Single Cell Trajectories
# 5.3.1 choosing genes that define progress
expressed_genes <- row.names(subset(fData(cds), num_cells_expressed >= 10))
length(expressed_genes)
#diff_test_res <- differentialGeneTest(cds[expressed_genes,],
#                                      fullModelFormulaStr = "~ manual",
#                                      cores = 4) #takes long time
#ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))

# get DE genes from Seurat
markers.Epi = read.csv(paste0("output/20190311/markers.Epi.csv"))

# a. All markers identified by Seurat for epithelial clusters (in epithelial-only analysis; let’s call it list A)
ordering_genes = unique(markers.Epi$gene)
#b. A subset of most significant differentially expressed markers in the list A: p val adj <10(-10)
ordering_genes = unique(markers.Epi[markers.Epi$p_val_adj<0.1^10,"gene"])
#c. Transcription factors (TFs) in the list A
TF = c("HMGA1", "TP63", "SOX15", "SOX7", "ELF3", "SPDEF", "HES4", "RARRES1",
       "XBP1", "SOX4", "FOXJ1", "RFX2", "TP73", "MYB", "SOX2", "HOPX", "NKX2-1", "ETV5", "FOXA2")
ordering_genes = FilterGenes(epithelial,TF)

length(ordering_genes)
cds <- setOrderingFilter(cds, ordering_genes)
plot_ordering_genes(cds)

# 5.3.2 reduce data dimensionality
#Now we're ready to try clustering the cells:.
plot_pc_variance_explained(cds, return_all = F) # norm_method = 'log',
cds <- reduceDimension(cds, max_components = 2,
                                  method = 'DDRTree')
#Trajectory step 3: order cells along the trajectory
cds <- orderCells(cds)

g1 <- plot_cell_trajectory(cds, color_by = "manual",cell_size = 3)
g2 <- plot_cell_trajectory(cds, color_by = "State",cell_size = 3)

jpeg(paste0(path,"a-trajectory_epithelial.jpeg"), units="in", width=10, height=7,res=600)
g1
dev.off()

jpeg(paste0(path,"a-trajectory_epithelial_state.jpeg"), units="in", width=10, height=7,res=600)
g2
dev.off()
