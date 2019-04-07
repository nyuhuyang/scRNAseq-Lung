library(Seurat)
#downgrade R to 3.4.4 https://cloud.r-project.org/bin/macosx/el-capitan/base/R-3.4.4.pkg
library(URD)
source("../R/Seurat_functions.R")
path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
#====== 6.1 convert from seurat to URD  ==========================================
# load Seurat
(load(file="data/epithelial_Harmony_5_20190123.Rda"))
object <- epithelial
remove(epithelial);GC()
object %<>% SetAllIdent("res.0.8")
colors <- ExtractMetaColor(object)
object <- RunUMAP(object = object, dims = 1:30,reduction.use = "MNN",
                  reduction.name = "tsne", reduction.key = "tSNE_")
URD_object <- seuratToURD(object)
remove(object);GC()
pcSDPlot(URD_object)
jpeg(paste0(path,"UMAP_Stage_manual.jpeg"), units="in", width=10, height=7,
     res=600)
plotDim(URD_object,"manual",plot.title = "UMAP: Stage",label.clusters = T,
        discrete.colors=colors,legend = F,point.size = 2)+
        xlab("UMAP_1")+ylab("UMAP_2")+
        theme(plot.title = element_text(hjust = 0.5))
dev.off()
#======= 6.2 Calculate Diffusion Map =============
# In this case, knn=100 (larger than sqrt(n.cells)) works well because there are not many cell types.
# Sigma 16 is slightly smaller than the sigma auto-determined by using NULL parameter.
system.time(URD_object <- calcDM(URD_object, knn = 100, sigma=16))
plotDimArray(URD_object, reduction.use = "dm", dims.to.plot = 1:8, 
             outer.title = "Diffusion Map (Sigma 16, 100 NNs): Stage", 
             label="res.0.8", plot.title="", legend=F)
jpeg(paste0(path,"Developmental_stage_umap.jpeg"), units="in", width=10, height=7,
     res=600)
plotDim(URD_object, "manual", transitions.plot = 10000, 
        plot.title="Developmental stage (with transitions)",
        discrete.colors=colors,legend = F)+
        xlab("UMAP_1")+ylab("UMAP_2")+
        theme(plot.title = element_text(hjust = 0.5))
dev.off()
# ========6.3 Calculate pseudotime ===================
# Here we use all cells from the first stage as the root
root.cells <- cellsInCluster(URD_object, "manual", "Basal cells")

# Then we run 'flood' simulations
URD_floods <- floodPseudotime(URD_object, root.cells = root.cells, n=50, minimum.cells.flooded = 2, verbose=F)

# The we process the simulations into a pseudotime
URD_object <- floodPseudotimeProcess(URD_object, URD_floods, floods.name="pseudotime")

pseudotimePlotStabilityOverall(URD_object)

jpeg(paste0(path,"pseudotime_umap.jpeg"), units="in", width=10, height=7,
     res=600)
plotDim(URD_object, "pseudotime")+
        xlab("UMAP_1")+ylab("UMAP_2")+
        theme(plot.title = element_text(hjust = 0.5))
dev.off()

jpeg(paste0(path,"plotDists.jpeg"), units="in", width=10, height=7,res=600)
plotDists(URD_object, "pseudotime", "manual", plot.title="Pseudotime by stage")
dev.off()
# =======6.4 Find tips =======================
# Create a subsetted object of just those cells from the final stage
final_stage <- urdSubset(URD_object, cells.keep=cellsInCluster(URD_object, "manual",
                                                                 "Alveolar type II cells"))

# Use the variable genes that were calculated only on the final group of stages (which
# contain the last stage).
#final_stage@var.genes <- var.by.stage[[4]]

# Calculate PCA and tSNE
final_stage <- calcPCA(final_stage, mp.factor = 1.5)
pcSDPlot(final_stage)

set.seed(20)
final_stage <- calcTsne(final_stage)

# Calculate graph clustering of these cells
final_stage <- graphClustering(final_stage, num.nn = 50, do.jaccard=T, method="Louvain")
# By plotting the expression of marker genes, we can determine that cluster 1 is the notochord and cluster 2 is the prechordal plate.
plotDim(final_stage, "Louvain-50", plot.title = "Louvain (50 NN) graph clustering", point.size=3)
plotDim(final_stage, "SFTPC", plot.title="SFTPC")
plotDim(final_stage, "SFTPA1", plot.title="SFTPA1")


#Biased random walks

# Copy cluster identities from final_stage object to a new clustering ("tip.clusters") in the full URD_object object.
URD_object@group.ids[rownames(final_stage@group.ids), "tip.clusters"] <- final_stage@group.ids$`Louvain-50`

# Determine the parameters of the logistic used to bias the transition probabilities. The procedure
# is relatively robust to this parameter, but the cell numbers may need to be modified for larger
# or smaller data sets.
URD_object.ptlogistic <- pseudotimeDetermineLogistic(URD_object, "pseudotime", optimal.cells.forward=20, max.cells.back=40, do.plot = T)
# Bias the transition matrix acording to pseudotime
URD_object.biased.tm <- as.matrix(pseudotimeWeightTransitionMatrix(URD_object, "pseudotime", logistic.params=URD_object.ptlogistic))

# Simulate the biased random walks from each tip
URD_object.walks <- simulateRandomWalksFromTips(URD_object, tip.group.id="tip.clusters", root.cells=root.cells, transition.matrix = URD_object.biased.tm, n.per.tip = 25000, root.visits = 1, max.steps = 5000, verbose = F)

# Process the biased random walks into visitation frequencies
URD_object <- processRandomWalksFromTips(URD_object, URD_object.walks, verbose = F)


# We can then visualize the tips and the visitation of cells from each tip on the dataset.

plotDim(URD_object, "tip.clusters", plot.title="Cells in each tip")

plotDim(URD_object, "visitfreq.log.1", plot.title="Visitation frequency from tip 1 (log10)", transitions.plot=10000)
plotDim(URD_object, "visitfreq.log.2", plot.title="Visitation frequency from tip 2 (log10)", transitions.plot=10000)

#####################
# Build tree
#####################

# Load the cells used for each tip into the URD object
URD_object.tree <- loadTipCells(URD_object, "tip.clusters")

# Build the tree
URD_object.tree <- buildTree(URD_object.tree, pseudotime = "pseudotime", tips.use=1:2, divergence.method = "preference", cells.per.pseudotime.bin = 25, bins.per.pseudotime.window = 8, save.all.breakpoint.info = T, p.thresh=0.001)

# Name the segments based on our previous determination of the identity of tips 1 and 2.
URD_object.tree <- nameSegments(URD_object.tree, segments=c("1","2"), segment.names = c("Notochord", "Prechordal Plate"), short.names = c("Noto", "PCP"))

plotTree(URD_object.tree, "manual", title="Developmental Stage")

plotTree(URD_object.tree, "GSC", title="GSC (early prechordal plate marker)")

plotTree(URD_object.tree, "NOTO", title="NOTO (early notochord marker)")


plotTree(URD_object.tree, "HE1A", title="HE1A (prechordal plate differentiation marker")

plotTree(URD_object.tree, "COL8A1A", title="COL8A1A (notochord differentiation marker")


#Additionally, we can refer back to the tSNE representation to see where the branchpoint was found.

plotTree(URD_object.tree, "segment", title="URD tree segment")

plotDim(URD_object.tree, "segment", plot.title="URD tree segment")

# Force-directed layout
# Generate the force-directed layout
URD_object.tree <- treeForceDirectedLayout(URD_object.tree, num.nn=100, cut.unconnected.segments=2, verbose=T)

plotTreeForce(URD_object.tree, "SFTPA1", title = "SFTPA1", title.cex = 2, title.line=2.5)

plotTreeForce(URD_object.tree, "SFTPC", title = "SFTPC", title.cex=2, title.line=2.5)

plotTreeForce(URD_object.tree, "COL8A1A", title="COL8A1A", title.cex=2, title.line=2.5)