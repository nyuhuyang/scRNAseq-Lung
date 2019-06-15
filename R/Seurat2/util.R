require(tidyr)

DotPlot.1 <- function (object, genes.plot, cols.use = c("lightgrey", "blue"), 
          col.min = -2.5, col.max = 2.5, dot.min = 0, dot.scale = 6, 
          scale.by = "radius", scale.min = NA, scale.max = NA, group.by, 
          plot.legend = FALSE, do.return = FALSE, x.lab.rot = FALSE,size=5 ) 
{
        scale.func <- switch(EXPR = scale.by, size = scale_size, 
                             radius = scale_radius, stop("'scale.by' must be either 'size' or 'radius'"))
        if (!missing(x = group.by)) {
                object <- SetAllIdent(object = object, id = group.by)
        }
        data.to.plot <- data.frame(FetchData(object = object, vars.all = genes.plot))
        colnames(x = data.to.plot) <- genes.plot
        data.to.plot$cell <- rownames(x = data.to.plot)
        data.to.plot$id <- object@ident
        data.to.plot <- data.to.plot %>% gather(key = genes.plot, 
                                                value = expression, -c(cell, id))
        data.to.plot <- data.to.plot %>% group_by(id, genes.plot) %>% 
                summarize(avg.exp = mean(expm1(x = expression)), pct.exp = Seurat:::PercentAbove(x = expression, 
                                                                                        threshold = 0))
        data.to.plot <- data.to.plot %>% ungroup() %>% group_by(genes.plot) %>% 
                mutate(avg.exp.scale = scale(x = avg.exp)) %>% mutate(avg.exp.scale = MinMax(data = avg.exp.scale, 
                                                                                             max = col.max, min = col.min))
        data.to.plot$genes.plot <- factor(x = data.to.plot$genes.plot, 
                                          levels = rev(x = genes.plot))
        data.to.plot$pct.exp[data.to.plot$pct.exp < dot.min] <- NA
        data.to.plot$pct.exp <- data.to.plot$pct.exp * 100
        p <- ggplot(data = data.to.plot, mapping = aes(x = genes.plot, 
                                                       y = id)) + geom_point(mapping = aes(size = pct.exp, 
                                                                                           color = avg.exp.scale)) + scale.func(range = c(0, dot.scale), 
                                                                                                                                limits = c(scale.min, scale.max)) + theme(axis.title.x = element_blank(), 
                                                                                                                                                                          axis.title.y = element_blank())
        if (length(x = cols.use) == 1) {
                p <- p + scale_color_distiller(palette = cols.use)
        }
        else {
                p <- p + scale_color_gradient(low = cols.use[1], high = cols.use[2])
        }
        if (!plot.legend) {
                p <- p + theme(legend.position = "none")
        }
        if (x.lab.rot) {
                p <- p + theme(axis.text.x = element_text(angle = 90, 
                                                          vjust = 0.5,
                                                          size=size))
        }
        suppressWarnings(print(p))
        if (do.return) {
                return(p)
        }
}


#' find marker across by conditions
#' Modified FindMarkers.UMI function, compare the same ident across conditions
#' @param ident.1 dent.1 list
#' @param ident.2 dent.2 list
#' @param ... all paramethers are the same as Seurat::FindMarkers
#' @export gde.all data frame
#' @export save.path folder to save
#' @example FindPairMarkers(object, ident.1 = 1:8, ident.2 = c(5:8,1:4))
FindPairMarkers <- function(object, ident.1, ident.2 = NULL, genes.use = NULL,return.thresh = 0.05,
                            logfc.threshold = 0.05, test.use = "MAST", min.pct = 0.1,
                            min.diff.pct = -Inf, print.bar = TRUE, only.pos = FALSE,
                            max.cells.per.ident = Inf, random.seed = 1, latent.vars = "nUMI",
                            min.cells.gene = 3,min.cells.group=3, pseudocount.use = 1, 
                            assay.type = "RNA",save.path = NULL,save.files = TRUE,...){
        #prepare save folder name
        if(is.null(save.path)){
                path <- paste0("output/",gsub("-","",Sys.Date()),"/")
                save.path <- paste0(path,ident1,"_vs_",ident2,"/")
        }
        gde <- list()
        for(i in 1:length(ident.1)) {
                if(class(ident.1) =="list") {
                        ident.1vs2 <- paste(paste(ident.1[[i]],collapse = "_"),
                                            paste(ident.2[[i]], collapse = "_"),
                                            sep = " vs ")
                        ident1 <- ident.1[[i]]; ident2 <- ident.2[[i]]
                } else {
                        ident.1vs2 <- paste(ident.1[i], ident.2[i], sep = " vs ")
                        ident1 <- ident.1[i]; ident2 <- ident.2[i]
                        }
                print(ident.1vs2)
                gde[[i]] <- FindMarkers.UMI(object = object, 
                                            ident.1 = ident1,
                                            ident.2 = ident2, assay.type = assay.type, 
                                            genes.use = genes.use, 
                                            logfc.threshold = logfc.threshold, test.use = test.use, 
                                            min.pct = min.pct, min.diff.pct = min.diff.pct, 
                                            print.bar = print.bar, only.pos = only.pos, min.cells.gene = min.cells.gene, 
                                            min.cells.group = min.cells.group, latent.vars = latent.vars, 
                                            max.cells.per.ident = max.cells.per.ident)
                gde[[i]] <- gde[[i]][order(-gde[[i]]$avg_logFC,gde[[i]]$p_val),]
                gde[[i]] <- subset(x = gde[[i]], subset = p_val < return.thresh)
                gde[[i]]$cluster1.vs.cluster2 <- ident.1vs2
                gde[[i]]$gene <- rownames(x = gde[[i]])
                if(save.files){
                        if(!dir.exists(save.path)) dir.create(save.path, recursive = T)
                        write.csv( gde[[i]], paste0(save.path,ident.1vs2,".csv"))
                }
        }
        return(bind_rows(gde))
}

RenameCells.1 <-  function(object, add.cell.id = NULL, new.names = NULL, 
                           for.merge = FALSE, to.upper =TRUE) 
{
        if (!is.null(add.cell.id)) {
                new.cell.names <- paste(add.cell.id, object@cell.names, 
                                        sep = "_")
                new.rawdata.names <- paste(add.cell.id, colnames(object@raw.data), 
                                           sep = "_")
        }
        if (to.upper){
                new.cell.names <- toupper(object@cell.names)
                new.rawdata.names <- toupper(colnames(object@raw.data))
        }
        colnames(object@raw.data) <- new.rawdata.names
        rownames(object@meta.data) <- new.cell.names
        object@cell.names <- new.cell.names
        if (for.merge) {
                return(object)
        }
        colnames(object@data) <- new.cell.names
        if (!is.null(object@scale.data)) {
                colnames(object@scale.data) <- new.cell.names
        }
        names(object@ident) <- new.cell.names
        if (length(object@dr) > 0) {
                for (dr in names(object@dr)) {
                        rownames(object@dr[[dr]]@cell.embeddings) <- new.cell.names
                }
        }
        if (nrow(object@snn) == length(new.cell.names)) {
                colnames(object@snn) <- new.cell.names
                rownames(object@snn) <- new.cell.names
        }
        if (!is.null(object@kmeans)) {
                if (!is.null(object@kmeans@gene.kmeans.obj)) {
                        colnames(object@kmeans@gene.kmeans.obj$centers) <- new.cell.names
                }
                if (!is.null(object@kmeans@cell.kmeans.obj)) {
                        names(object@kmeans@cell.kmeans.obj$cluster) <- new.cell.names
                }
        }
        return(object)
}


DoHeatmap_exp <- function (object, data.use = NULL, use.scaled = TRUE, cells.use = NULL, 
                           genes.use = NULL, disp.min = -2.5, disp.max = 2.5, group.by = "ident", 
                           group.order = NULL, draw.line = TRUE, col.low = "#FF00FF", 
                           col.mid = "#000000", col.high = "#FFFF00", slim.col.label = FALSE, 
                           remove.key = FALSE, rotate.key = FALSE, title = NULL, cex.col = 10, 
                           cex.row = 10, group.label.loc = "bottom", group.label.rot = FALSE, 
                           group.cex = 15, group.spacing = 0.15, assay.type = "RNA", 
                           do.plot = TRUE) 
{
        if (is.null(x = data.use)) {
                if (use.scaled) {
                        data.use <- GetAssayData(object, assay.type = assay.type, 
                                                 slot = "scale.data")
                }
                else {
                        data.use <- GetAssayData(object, assay.type = assay.type, 
                                                 slot = "data")
                }
        }
        cells.use <- Seurat:::SetIfNull(x = cells.use, default = object@cell.names)
        cells.use <- intersect(x = cells.use, y = colnames(x = data.use))
        if (length(x = cells.use) == 0) {
                stop("No cells given to cells.use present in object")
        }
        genes.use <- Seurat:::SetIfNull(x = genes.use, default = rownames(x = data.use))
        genes.use <- intersect(x = genes.use, y = rownames(x = data.use))
        if (length(x = genes.use) == 0) {
                stop("No genes given to genes.use present in object")
        }
        if (is.null(x = group.by) || group.by == "ident") {
                cells.ident <- object@ident[cells.use]
        }
        else {
                cells.ident <- factor(x = FetchData(object = object, 
                                                    cells.use = cells.use, vars.all = group.by)[, 1])
                names(x = cells.ident) <- cells.use
        }
        cells.ident <- factor(x = cells.ident, labels = intersect(x = levels(x = cells.ident), 
                                                                  y = cells.ident))
        data.use <- data.use[genes.use, cells.use, drop = FALSE]
        if ((!use.scaled)) {
                data.use = as.matrix(x = data.use)
                if (disp.max == 2.5) 
                        disp.max = 10
        }
        data.use <- MinMax(data = data.use, min = disp.min, max = disp.max)
        data.use <- as.data.frame(x = t(x = data.use))
        data.use$cell <- rownames(x = data.use)
        colnames(x = data.use) <- make.unique(names = colnames(x = data.use))
        data.use <- data.use %>% melt(id.vars = "cell")
        names(x = data.use)[names(x = data.use) == "variable"] <- "gene"
        names(x = data.use)[names(x = data.use) == "value"] <- "expression"
        data.use$ident <- cells.ident[data.use$cell]
        if (!is.null(group.order)) {
                if (length(group.order) == length(levels(data.use$ident)) && 
                    all(group.order %in% levels(data.use$ident))) {
                        data.use$ident <- factor(data.use$ident, levels = group.order)
                }
                else {
                        stop("Invalid group.order")
                }
        }
        data.use$gene <- with(data = data.use, expr = factor(x = gene, 
                                                             levels = rev(x = unique(x = data.use$gene))))
        data.use$cell <- with(data = data.use, expr = factor(x = cell, 
                                                             levels = cells.use))
        spred_data.use <- data.use %>% spread(gene, expression)
        rownames(spred_data.use) = spred_data.use$cell
        spred_data.use = spred_data.use[,-1]
        colnames(spred_data.use)[1] = "cell_type"
        spred_data.use = spred_data.use[,ncol(spred_data.use):1]
        
        sortbyfactor <- function(DataFrame, key="cell_type", levels = group.order){
                DataFrame$order <- plyr::mapvalues(x = DataFrame[,key],
                                                   from = group.order,
                                                   to =   unique(DataFrame[,key]))
                DataFrame = DataFrame[order(DataFrame$order),]
                DataFrame = DataFrame[,-ncol(DataFrame)]
                return(DataFrame)
        }
        spred_data.use = sortbyfactor(spred_data.use,key="cell_type", levels = group.order)
        
        return(spred_data.use)
}


# make horizontal color bar==========
# make horizontal corlor bar for DoHeatmap
#' @param object Seurat object
#' @param group_by The criteria to group horizontal bar
#' @param split.by labels for each  horizontal bar
#' @param color color scheme
#' @export g vertical ggplot
#' @example MakeHCorlorBar(B_cells_MCL,group_by = "X5_clusters")

MakeHCorlorBar <- function(object, group_by = "X5_clusters", remove.legend = F,
                           file_name = "DE_clusters_top50",
                           split.by = "singler1sub",do.print = TRUE,do.return=FALSE){
        
        (group_index <- unique(object@meta.data[,group_by]) %>% sort)
        if(!paste0(split.by,".colors") %in% colnames(object@meta.data)){
                # add color to cluster
                gg_color_hue <- function(n) {
                        hues = seq(15, 375, length = n + 1)
                        grDevices::hcl(h = hues, l = 65, c = 100)[1:n]
                }
                color_num <- length(unique(object@meta.data[,split.by]))
                object <- AddMetaColor(object, label= split.by, 
                                       colors = gg_color_hue(color_num))
                
        }
        for(ind in group_index){
                b_color_bar <- object@meta.data[(object@meta.data[,group_by] %in% ind),]
                b_color_bar$x <- 1:nrow(b_color_bar)

                col_index <- which(colnames(b_color_bar) %in% c(split.by,paste0(split.by,".colors")))
                
                b_color_bar[,col_index] = sapply(b_color_bar[,col_index],as.character)
                b_color_bar$labels = b_color_bar[,split.by]
                g <- ggplot(data = b_color_bar, aes(x, split.by, fill = labels)) +
                        geom_tile()+
                        theme_bw() +
                        theme(panel.border = element_blank(),
                              panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank(),
                              axis.line = element_blank(),
                              axis.title.x=element_blank(),
                              axis.text.x=element_blank(),
                              axis.ticks.x=element_blank(),
                              axis.title.y=element_blank(),
                              axis.text.y=element_blank(),
                              axis.ticks.y=element_blank())
                if(remove.legend) g = g + theme(legend.position="none")
                
                # prepare color scheme
                dup <- duplicated(b_color_bar[,col_index[1]])
                colors_df <- b_color_bar[!dup,col_index]
                colors_df <- colors_df[order(colors_df[,1]),]
                colors_fill = colors_df[,paste0(split.by,".colors")]
                g = g + scale_fill_manual(values = colors_fill)
                
                path <- paste0("output/",gsub("-","",Sys.Date()),"/",file_name,"/")
                if(!dir.exists(path)) dir.create(path, recursive = T)
                if(do.print) {
                        jpeg(paste0(path,ind,"_","horizontal_bar.jpeg"), units="in",width=10, height=7,res=600)
                        print(g)
                        dev.off()
                } else if(do.return) {
                        return(g)
                }
        }
}
