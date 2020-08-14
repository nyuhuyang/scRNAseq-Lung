#' Phylogenetic Analysis of Identity Classes
#' @param method agglomeration method to be used. parameter pass to hclust function
BuildClusterTree <- function (object, assay = NULL, features = NULL, dims = NULL, 
          graph = NULL, slot = "data", reorder = FALSE, reorder.numeric = FALSE, 
          method = "complete", verbose = TRUE) 
{
    assay <- assay %||% DefaultAssay(object = object)
    if (!is.null(x = graph)) {
        idents <- levels(x = object)
        nclusters <- length(x = idents)
        data.dist <- matrix(data = numeric(length = 1L), nrow = nclusters, 
                            ncol = nclusters, dimnames = list(idents, idents))
        graph <- object[[graph]]
        cxi <- CellsByIdentities(object = object)
        cpairs <- na.omit(object = unique(x = t(x = apply(X = expand.grid(1:nclusters, 
                                                                          1:nclusters)[, c(2, 1)], MARGIN = 1, FUN = function(x) {
                                                                              if (length(x = x) == length(x = unique(x = x))) {
                                                                                  return(sort(x = x))
                                                                              }
                                                                              return(c(NA, NA))
                                                                          }))))
        if (verbose) {
            pb <- txtProgressBar(style = 3, file = stderr())
        }
        for (i in 1:nrow(x = cpairs)) {
            i1 <- cpairs[i, ][1]
            i2 <- cpairs[i, ][2]
            graph.sub <- graph[cxi[[idents[i1]]], cxi[[idents[i2]]]]
            d <- mean(x = graph.sub)
            if (is.na(x = d)) {
                d <- 0
            }
            data.dist[i1, i2] <- d
            if (verbose) {
                setTxtProgressBar(pb = pb, value = i/nrow(x = cpairs))
            }
        }
        if (verbose) {
            close(con = pb)
        }
        diag(x = data.dist) <- 1
        data.dist <- dist(x = data.dist)
    }
    else if (!is.null(x = dims)) {
        my.lapply <- ifelse(test = verbose, yes = pblapply, 
                            no = lapply)
        embeddings <- Embeddings(object = object, reduction = "pca")[, 
                                                                     dims]
        data.dims <- my.lapply(X = levels(x = object), FUN = function(x) {
            cells <- WhichCells(object = object, idents = x)
            if (length(x = cells) == 1) {
                cells <- c(cells, cells)
            }
            temp <- colMeans(x = embeddings[cells, ])
        })
        data.dims <- do.call(what = "cbind", args = data.dims)
        colnames(x = data.dims) <- levels(x = object)
        data.dist <- dist(x = t(x = data.dims))
    }
    else {
        features <- features %||% VariableFeatures(object = object)
        features <- intersect(x = features, y = rownames(x = object))
        data.avg <- AverageExpression(object = object, assays = assay, 
                                      features = features, slot = slot, verbose = verbose)[[1]]
        data.dist <- dist(x = t(x = data.avg[features, ]))
    }
    data.tree <- ape::as.phylo(x = hclust(d = data.dist, method = method))
    Tool(object = object) <- data.tree
    if (reorder) {
        if (verbose) {
            message("Reordering identity classes and rebuilding tree")
        }
        old.ident.order <- levels(x = object)
        data.tree <- Tool(object = object, slot = "BuildClusterTree")
        all.desc <- GetDescendants(tree = data.tree, node = (data.tree$Nnode + 
                                                                 2))
        all.desc <- old.ident.order[all.desc[all.desc <= (data.tree$Nnode + 
                                                              1)]]
        Idents(object = object) <- factor(x = Idents(object = object), 
                                          levels = all.desc, ordered = TRUE)
        if (reorder.numeric) {
            new.levels <- sort(x = unique(x = as.integer(x = Idents(object = object))))
            Idents(object = object) <- factor(x = as.integer(x = Idents(object = object)), 
                                              levels = new.levels)
            object[["tree.ident"]] <- as.integer(x = Idents(object = object))
        }
        object <- BuildClusterTree(object = object, assay = assay, 
                                   features = features, dims = dims, graph = graph, 
                                   slot = slot, reorder = FALSE, verbose = verbose)
    }
    return(object)
}
