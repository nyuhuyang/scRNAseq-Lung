xCellAnalysis.1 <- function (expr, signatures = NULL, genes = NULL, spill = NULL, 
          rnaseq = TRUE, file.name = NULL, scale = TRUE, alpha = 0.5, 
          save.raw = FALSE, parallel.sz = 4, parallel.type = "SOCK", 
          cell.types.use = names(signatures)) 
{
        if (is.null(signatures)) 
                signatures = xCell.data$signatures
        if (is.null(genes)) 
                genes = xCell.data$genes
        if (is.null(spill)) {
                if (rnaseq == TRUE) {
                        spill = xCell.data$spill
                }
                else {
                        spill = xCell.data$spill.array
                }
        }
        if (is.null(file.name) || save.raw == FALSE) {
                fn <- NULL
        } else {
                fn <- paste0(file.name, "_RAW.txt")
        }
        #if (!is.null(cell.types.use)) {
        #        A = intersect(cell.types.use, rownames(spill$K))
        #        if (length(A) < length(cell.types.use)) {
        #                return("ERROR - not all cell types listed are available")
        #        }
        #}
        scores <- rawEnrichmentAnalysis.1(expr, signatures, genes, 
                                        fn, parallel.sz = parallel.sz, parallel.type = "SOCK")
        #@scores.transformed <- transformScores(scores, spill$fv, 
         #                                     scale)
        #scores.transformed <- scores
        #if (is.null(file.name)) {
        #        fn <- NULL
        #} else {
        #        fn <- file.name
        #}
        #if (is.null(cell.types.use)) {
        #        scores.adjusted <- spillOver(scores.transformed, spill$K, 
        #                                     alpha, fn)
        #        scores.adjusted = microenvironmentScores(scores.adjusted)
        #} else {
        #        scores.adjusted <- spillOver(scores.transformed[cell.types.use, 
        #                                                        ], spill$K, alpha, fn)
        #}
        return(scores)
}


rawEnrichmentAnalysis.1 <- function (expr, signatures, genes, file.name = NULL, parallel.sz = 4, 
          parallel.type = "SOCK") 
{
        shared.genes <- intersect(rownames(expr), genes)
        print(paste("Num. of genes:", length(shared.genes)))
        expr <- expr[shared.genes, ]
        if (dim(expr)[1] < 5000) {
                print(paste("ERROR: not enough genes"))
                return - 1
        }
        expr <- apply(expr, 2, rank)
        if (packageVersion("GSVA") >= "1.36.0") {
                scores <- GSVA::gsva(expr, signatures, method = "ssgsea", 
                                     ssgsea.norm = FALSE, parallel.sz = parallel.sz)
        }
        else {
                scores <- GSVA::gsva(expr, signatures, method = "ssgsea", 
                                     ssgsea.norm = FALSE, parallel.sz = parallel.sz, 
                                     parallel.type = parallel.type)
        }
        scores = scores - apply(scores, 1, min)
        cell_types <- unlist(strsplit(rownames(scores), "%"))
        #cell_types <- cell_types[seq(1, length(cell_types), 3)]
        agg <- aggregate(scores ~ cell_types, FUN = mean)
        rownames(agg) <- agg[, 1]
        scores <- agg[, -1]
        if (!is.null(file.name)) {
                write.table(scores, file = file.name, sep = "\t", col.names = NA, 
                            quote = FALSE)
        }
        scores
}

