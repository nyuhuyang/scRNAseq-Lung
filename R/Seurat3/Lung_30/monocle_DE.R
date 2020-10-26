########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################
invisible(lapply(c("Seurat","dplyr","cowplot",
                   "magrittr","data.table","future","ggplot2","tidyr"), function(x) {
                       suppressPackageStartupMessages(library(x,character.only = T))
                   }))
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat3_functions.R")

path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)
# load files
read.path = "Yang/Lung_30/Monocle2/20-ReadDE-distal.terminal.proximal.COPD-root=BC/"
cds = readRDS(paste0(read.path,basename(read.path),"_cds.rds"))

object <- as.Seurat(cds)
meta.data = object@meta.data

object[["SCT"]] = object[["RNA"]]

write.csv(as.data.frame.matrix(table(object$State, object$conditions)),
          file = paste0(path,"cellCounts_groups_state.csv"))
object$Pseudotime_group <- cut(object$Pseudotime,
                               breaks = c(0, 2, 13, 14, 17, Inf),
                               labels = c("Pseudotime < 2",
                                          "Pseudotime >= 2 & < 13",
                                          "Pseudotime >= 13 & < 14",
                                          "Pseudotime >= 14 & < 17",
                                          "Pseudotime >= 17"),
                               right  = FALSE)

write.csv(as.data.frame.matrix(table(object$State, object$Pseudotime_group)),
          file = paste0(path,"cellCounts_State_Pseudotime.csv"))