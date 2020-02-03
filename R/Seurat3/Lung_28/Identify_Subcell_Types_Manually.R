########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################
invisible(lapply(c("Seurat","dplyr","kableExtra","cowplot",
                   "magrittr","MAST","future","ggplot2","tidyr"), function(x) {
                           suppressPackageStartupMessages(library(x,character.only = T))
                   }))
source("../R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)

groups <- c("AT","B","D","En","F","Mon","SAE","SMG","SMP","T")sdd
#======== rename ident =================
# AT
args = 1
(g <- groups[args])
(load(file = paste0("data/Lung_28_",g,"_20200126.Rda")))
Idents(object) = "integrated_snn_res.0.1"
object %<>% RenameIdents("0" = "AT2-a",
                         "1" = "AT2-b",
                         "2" = "AT2-c",
                         "3" = "AT1")
object[["cell_types"]] = as.character(Idents(object))

#======== rename ident =================
# D
args = 3
(g <- groups[args])
(load(file = paste0("data/Lung_28_",g,"_20200126.Rda")))
Idents(object) = "integrated_snn_res.1.4"
object %<>% RenameIdents("0" = #"AT2-a",
                         "1" = #"AT2-b",
                         "2" = #"AT2-c",
                         "3" = #"AT1",
                         "4" = #"AT2-a",
                         "5" = "BC",
                         "6" = #"AT2-c",
                         "7" = #"AT1",
                         "8" = #"AT2-a",
                         "9" = #"AT2-b",
                         "10" =# "AT2-c",
                         "11" =# "AT1",
                         "12" =# "AT1",)
object[["cell_types"]] = as.character(Idents(object))

