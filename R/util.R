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
