# Install devtools from CRAN
install.packages("devtools")
# Use devtools to install hdf5r and loomR from GitHub
install.packages("hdf5r")
devtools::install_github(repo = "mojaveazure/loomR", ref = "develop")

library(Seurat) #v2.3.4 main
library(magrittr)
source("../R/Seurat_functions.R")
source("R/util.R")
(load(file="data/epithelial_Harmony_5_20190123.Rda"))
object <- epithelial
remove(epithelial);GC()

object <- RenameCells.1(object, to.upper =TRUE)
object <- RunUMAP(object = object, dims = 1:30,reduction.use = "MNN")
object %<>% SetAllIdent("manual")
g_umap <- DimPlot.1(object = object, do.label = F, group.by = "ident", 
                    reduction.use = "umap",
                    do.return = TRUE, no.legend = F, 
                    cols.use = ExtractMetaColor(object),
                    pt.size = 1,label.size = 6 )+
        ggtitle("UMAP plot of all epithelial cell types")+
        theme(text = element_text(size=15),							
              plot.title = element_text(hjust = 0.5,size = 18)) 

jpeg(paste0(path,"UMAplot-epithelial.jpeg"), units="in", width=10, height=7,res=600)
print(g_umap)
dev.off()

pfile <- Convert(from = object, to = "anndata", filename = paste0("data/h5ad/epithelial.h5ad"), 
                 display.progress = F)
remove(pfile);GC()