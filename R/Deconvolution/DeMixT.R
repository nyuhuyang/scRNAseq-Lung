source("http://bioconductor.org/biocLite.R")
biocLite("DSS", dependencies=TRUE)
install.packages("DeMix", repos="http://r-forge.r-project.org", dependencies=TRUE)
devtools::install_github("wwylab/DeMixT/DeMixT_0.1") 
devtools::install_github("wwylab/DeMixT/DeMixT_0.2")
library("DeMixT")
# Example 1: simulated two-component data 
data(test.data1)
res <- DeMixT(data.Y = test.data1$y, data.comp1 = test.data1$comp1, if.filter = FALSE, output.more.info = TRUE)
res$pi
head(res$decovExprT, 3)
head(res$decovExprN1, 3)
head(res$decovMu, 3)
head(res$decovSigma, 3)
res$pi.iter
res$gene.name


# Example 2: estimate proportions for simulated two-component data 
data(test.data1)
res <- DeMixT.S1(data.Y = test.data1$y, data.comp1 = test.data1$comp1, if.filter = FALSE)
res$pi
res$pi.iter
res$gene.name

# Example 3: two-component deconvolution given proportions 
data(test.data1)
givenpi <- c(t(as.matrix(test.data1$truth[-2,])))
res <- DeMixT.S2(data.Y = test.data1$y, data.comp1 = test.data1$comp1, givenpi = givenpi)
str(res)
