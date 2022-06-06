source("https://bioconductor.org/biocLite.R")
biocLite(c("EpiDISH","GEOquery"))
## ----download, eval=F, echo=T, message=FALSE, warning=FALSE----------------
require(GEOquery)
require(Biobase)
GSE80559 <- getGEO("GSE80559")
beta.m <- exprs(GSE80559[[1]])

## ----load, eval=TRUE, echo=T, message=FALSE, warning=FALSE-----------------
library(EpiDISH)
data(centDHSbloodDMC.m)
data(DummyBeta.m)

## ----infer, eval=TRUE, echo=T, message=FALSE, warning=FALSE----------------
ref.m <- centDHSbloodDMC.m[,1:6]
out.l <- epidish(beta.m, ref.m, method = "RPC") 

## ----check, eval=TRUE, echo=T, message=FALSE, warning=FALSE----------------
out.l$estF
dim(out.l$ref)
dim(out.l$dataREF)

## ----sessionInfo, echo=FALSE-----------------------------------------------
sessionInfo()

## CIBERSORT

LM22 <- read.delim2("./data/LM22.txt",row.names = 1)
out_LM22.l <- epidish(beta.m, LM22, method = "RPC") 
