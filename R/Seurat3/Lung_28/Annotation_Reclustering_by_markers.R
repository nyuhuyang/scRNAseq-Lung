########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################
invisible(lapply(c("Seurat","dplyr","cowplot","eulerr",
                   "magrittr","harmony"), function(x) {
                           suppressPackageStartupMessages(library(x,character.only = T))
                   }))
source("../R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)

object = readRDS(file = "data/Lung_28_Global_20200219.rds") 
data = object[["SCT"]]@data
remove(object);GC()
df_samples <- readxl::read_excel("doc/Annotations-4-30-20-RS.xlsx",
                                 sheet = "Sheet2")
Annotations <- list()
for(i in 1:51){
        Annotations[[i]] = readRDS(file = paste0(path, "Annotation-",i,"-",
                              df_samples$`Cell type`[i],".rds"))
}
names(Annotations) = df_samples$`Cell type`

pairs <- replicate(40,sample(1:49,8,replace = FALSE)) %>% 
        apply(2, as.integer) %>% split(rep(1:ncol(.), each = nrow(.)))
pairs <- list(c("BC","S-d","NEC1","MEC2","Sq","SMG-Ser","AT2","IC"),
              c("Nr","SM1","En-C","Cr","Gli","F3","F1","F4"),
              c("C1","H","S","Ion","Mon","M1","Neu-2","Neu-1"),
              c("Cr","MEC1","M0","Mon","M2","Neu-2","MEC2","S"),
              c("C1","p-C","AT1","C1","S-d","Sq","F1","Pr","SM1","Sq"))


for(k in 4:length(pairs)){
        uniform_intersections <- euler(Annotations[unique(pairs[[k]])])
        #jpeg(paste0(path,"Venn _",paste(df_samples$`Cell type`[pairs[[k]]],
        jpeg(paste0(path,"Venn_",paste(pairs[[k]],
                                        collapse = "-"),".jpeg"), 
             units="in", width=10, height=7,res=300)
        print(plot(uniform_intersections))
        dev.off()
        Progress(k, length(pairs))       
}


# Modify Annotations
#' @param Annotations annotation list with cell types as names and cell barcodes as elements
#' @param type1 cell type 1
#' @param type2 cell type 2
#' @param fun intersect, create a new cell type using intersect, 
#' @example Annotations %<>% Intersect(type1 = "F1", type2 = "F4")

Intersect <- function(Annotations, type1 = NULL, type2 = NULL){
        type3 = paste0(type1,"_",type2)
        Annotations[[type3]] = intersect(Annotations[[type1]], Annotations[[type2]])
        Annotations[[type1]] %<>% setdiff(Annotations[[type3]])
        Annotations[[type2]] %<>% setdiff(Annotations[[type3]])

        return(Annotations)
}
Annotations[["NEC1"]] %<>% union(Annotations[["NEC0"]])
Annotations = Annotations[-grep("NEC0", names(Annotations))]

Annotations[["S-d"]] %<>% union(Annotations[["S-d-1"]])
Annotations = Annotations[-grep("S-d-1", names(Annotations))]

Annotations[["S"]] %<>% setdiff(Annotations[["S-d"]])
Annotations[["S"]] %<>% setdiff(union(Annotations[["SMG-Muc"]], Annotations[["SMG-Ser"]]))
Annotations[["T-7s"]] %<>% setdiff(Annotations[["T-cn"]])
Annotations[["T-7s"]] %<>% setdiff(Annotations[["T-rm"]])
Annotations[["T-7s"]] %<>% setdiff(Annotations[["T-NK"]])
Annotations[["T-7s"]] %<>% setdiff(Annotations[["T-reg"]])
Annotations[["T-7s"]] %<>% setdiff(Annotations[["T-ifn"]])
Annotations[["T-7s"]] %<>% setdiff(Annotations[["B"]])
Annotations[["T-7s"]] %<>% setdiff(Annotations[["P-DC"]])
Annotations[["T-7s"]] %<>% setdiff(Annotations[["PC"]])

#  May 9, 2020

#AT2 - SMG-Ser
#if proximal, then SMG-Ser
#if distal or terminal, then AT2
(temp <- intersect(Annotations[["AT2"]], Annotations[["SMG-Ser"]]))

Annotations[["SMG-Ser"]] %<>% setdiff(temp[!grepl("-P_",temp)])
Annotations[["AT2"]] %<>% setdiff(temp[!grepl("-D_",temp)])

#AT2 - S-d 
#Priority: S-d
Annotations[["AT2"]] %<>% setdiff(Annotations[["S-d"]])

#BC - MEC2
#if ACTA2+ and proximal, then MEC2
#If not, then BC
(temp <- intersect(Annotations[["BC"]], Annotations[["MEC2"]]))
data["ACTA2",temp] >0
Annotations[["MEC2"]] %<>% setdiff(temp)
#BC - S-d
#Check-point: 
#If proximal and KRT5+ then BC
#if distal/terminal and SCGB3A2+ or distal/terminal and SFTPB+, then S-d
(temp <- intersect(Annotations[["BC"]], Annotations[["S-d"]]))
data["KRT5",temp] >0
Annotations[["S-d"]] %<>% setdiff(temp)

#BC - Sq
#Priority: Sq 
Annotations[["BC"]] %<>% setdiff(Annotations[["Sq"]])

#C1 - H
#If MUC5AC+, then H
(temp <- intersect(Annotations[["C1"]], Annotations[["H"]]))
data["MUC5AC",temp] >0
Annotations[["H"]] %<>% setdiff(temp[data["MUC5AC",temp] == 0])
Annotations[["C1"]] %<>% setdiff(temp[data["MUC5AC",temp] > 0])

#C1 - Ion
#Priority: Ion
Annotations[["C1"]] %<>% setdiff(Annotations[["Ion"]])

#Cr - En-C
#Priority - Cr
Annotations[["En-C"]] %<>% setdiff(Annotations[["Cr"]])

#Cr - Gli
#Priority - Gli
Annotations[["Cr"]] %<>% setdiff(Annotations[["Gli"]])

#Cr - SM1
#Priority - Cr
Annotations[["SM1"]] %<>% setdiff(Annotations[["Cr"]])

#En-C - Nr
#Priority - Nr
Annotations[["En-C"]] %<>% setdiff(Annotations[["Nr"]])

#F1 - F4
#Priority - F4
Annotations[["F1"]] %<>% setdiff(Annotations[["F4"]])

#F3 - Gli
#Priority - Gli
Annotations[["F3"]] %<>% setdiff(Annotations[["Gli"]])

#H - S
#If FOXJ1+ or CAPS+ or TPPP3+, then H
(temp <- intersect(Annotations[["H"]], Annotations[["S"]]))
data["FOXJ1",temp] >0 | data["CAPS",temp] >0 | data["TPPP3",temp] >0
Annotations[["S"]] %<>% setdiff(temp)

#IC - Sq
#Priority - Sq
Annotations[["IC"]] %<>% setdiff(Annotations[["Sq"]])

#IC - S-d
#If distal or terminal, and SCGB3A2+ or SFTPB+, and KRT5- and S100A2-, then S-d
#If proximal, then IC
(temp <- intersect(Annotations[["IC"]], Annotations[["S-d"]]))
(data["SCGB3A2",temp] >0 | data["SFTPB",temp] >0 ) & (data["KRT5",temp]  == 0 | data["S100A2",temp]  == 0)
Annotations[["S-d"]] %<>% setdiff(temp)

#IC - SMG-Ser
#Priority - SMG-Ser
Annotations[["IC"]] %<>% setdiff(Annotations[["SMG-Ser"]])

#Ion - S
#Priority - Ion
Annotations[["S"]] %<>% setdiff(Annotations[["Ion"]])

#M1 - Mon
#If EREG+, then Mon
#If MARCO+ or CXCL9+ or CXCL10+ or SPP1+, then M1
(temp <- intersect(Annotations[["M1"]], Annotations[["Mon"]]))
data["EREG",temp] >0
data["MARCO",temp] >0 | data["CXCL9",temp] >0 | data["CXCL10",temp] >0 |data["SPP1",temp] >0
Annotations[["M1"]] %<>% setdiff(temp[data["EREG",temp] >0])
(temp <- intersect(Annotations[["M1"]], Annotations[["Mon"]]))
data["EREG",temp] >0
data["MARCO",temp] >0 | data["CXCL9",temp] >0 | data["CXCL10",temp] >0 |data["SPP1",temp] >0
Annotations[["Mon"]] %<>% setdiff(temp)

#M1 Neu-1 or Neu-2
#If S100A8+ or S100A9+ or S100A12 or CSF3R+, then Neu-1 or Neu-2 
(temp <- intersect(Annotations[["M1"]], union(Annotations[["Neu-1"]],Annotations[["Neu-2"]])))
data["S100A8",temp] >0 | data["S100A9",temp] >0 | data["S100A12",temp] >0 | data["CSF3R",temp] >0
Annotations[["M1"]] %<>% setdiff(temp)

#NEC1 - S-d
#If distal or terminal, and SCGB3A2+ or SFTPB+, then S-d
#If CHGA+ or CHGB+ or GRP+, then NEC
(temp <- intersect(Annotations[["NEC1"]], Annotations[["S-d"]]))
data["SCGB3A2",temp] >0 | data["SFTPB",temp] >0 
Annotations[["NEC1"]] %<>% setdiff(temp)

#SMG-Ser - S-d
#If proximal, then SMG-Ser
#If distal or terminal, then S-d
(temp <- intersect(Annotations[["SMG-Ser"]], Annotations[["S-d"]]))
Annotations[["S-d"]] %<>% setdiff(temp[grepl("-P_",temp)])

Annotations = Annotations[sort(names(Annotations))]

# May 11, 2020, at 3:33 ==================
#AT1 S-d
#Checkpoint: SCGB3A2+ then S-d, AGER+ then AT1
(temp <- intersect(Annotations[["AT1"]], Annotations[["S-d"]]))
temp[data["SCGB3A2",temp] >0]
Annotations[["AT1"]] %<>% setdiff(temp)

#C1, p-C
#Checkpoint: HES6+ or CCNP+, then p-C
(temp <- intersect(Annotations[["C1"]], Annotations[["p-C"]]))
data["HES6",temp] >0 | data["CNTD2",temp] >0
Annotations[["C1"]] %<>% setdiff(temp[data["HES6",temp] > 0])
Annotations[["p-C"]] %<>% setdiff(temp[data["HES6",temp] == 0])

#C1 S-d
#Check-point: if distal or terminal, and SCGB3A2+ or SFTPB+, then S-d
(temp <- intersect(Annotations[["C1"]], Annotations[["S-d"]]))
data["SCGB3A2",temp] >0 | data["SFTPB",temp] >0
Annotations[["C1"]] %<>% setdiff(temp)

# Cr MEC1
# Checkpoint: if KRT5+ or KRT14+ or ACTA2+, then MEC1
(temp <- intersect(Annotations[["Cr"]], Annotations[["MEC1"]]))
data["KRT5",temp] >0 | data["KRT14",temp] >0 | data["ACTA2",temp] >0
Annotations[["MEC1"]] %<>% setdiff(temp)

#F1 SM1
#Checkpoint: if ACTA2 > DCN, then SM1
(temp <- intersect(Annotations[["F1"]], Annotations[["SM1"]]))
data["ACTA2",temp] > data["DCN",temp]
Annotations[["F1"]] %<>% setdiff(temp)

#M0 Mon
#Checkpoint: if EREG > MARCO, then Mon
(temp <- intersect(Annotations[["M0"]], Annotations[["Mon"]]))
data["EREG",temp] > data["MARCO",temp]
Annotations[["M0"]] %<>% setdiff(temp[data["EREG",temp] > data["MARCO",temp]])
Annotations[["Mon"]] %<>% setdiff(temp[data["EREG",temp] < data["MARCO",temp]])

#M2 Neu-2
#Checkpoint: if S100A8 > MARCO, then Neu-2
(temp <- intersect(Annotations[["M2"]], Annotations[["Neu-2"]]))
Annotations[["M2"]] %<>% setdiff(temp[data["S100A8",temp] > data["MARCO",temp]])
Annotations[["Neu-2"]] %<>% setdiff(temp[data["S100A8",temp] <= data["MARCO",temp]])

#MEC2 S
#Checkpoint: if ACTA2+ or KRT14+, then MEC2
(temp <- intersect(Annotations[["MEC2"]], Annotations[["S"]]))
Annotations[["S"]] %<>% setdiff(temp[data["ACTA2",temp] | data["KRT14",temp]])
Annotations[["MEC2"]] %<>% setdiff(temp[!(data["ACTA2",temp] | data["KRT14",temp])])

#NEC1 NEC2
# Keep NEC1 as priority
Annotations[["NEC2"]] %<>% setdiff(Annotations[["NEC1"]])

#Pr SM1
#Checkpoint: if LAMC3+ or KCNK3+ or CSPG4+, then Pr
(temp <- intersect(Annotations[["Pr"]], Annotations[["SM1"]]))
(rmove = temp[!temp %in% colnames(data)])
Annotations[["Pr"]] %<>% setdiff(rmove)
Annotations[["SM1"]] %<>% setdiff(rmove)
(temp <- intersect(Annotations[["Pr"]], Annotations[["SM1"]]))

table(Pr <- (data["LAMC3",temp[1:3]] >0 | data["KCNK3",temp]>0 | data["CSPG4",temp]>0))
Annotations[["SM1"]] %<>% setdiff(temp[Pr])
Annotations[["Pr"]] %<>% setdiff(temp[!Pr])

#S-d Sq
#Checkpoint: if proximal - Sq, if distal or terminal and SCGB3A2+ or SFTPB+, then S-d
(temp <- intersect(Annotations[["S-d"]], Annotations[["Sq"]]))
Annotations[["S-d"]] %<>% setdiff(temp[grepl("-P_",temp)])

Annotations = Annotations[sort(names(Annotations))]
# unlist and use.names without additional index in the names
Unlist <- function(List){
        for(i in seq_along(List)){
                names(List[[i]]) = rep(names(List)[i],length(List[[i]]))
        }
        flatten_List <- base::unlist(List,use.names=TRUE)
        names(flatten_List) %<>% gsub("\\..*","",.)
        return(flatten_List)
}
# fild overlap genes among list
#' @param barcodes list
#' @param return.pair.names TRUE/FALSE, return cell types overlap list
#' @param return.names TRUE/FALSE, return cell type names
#' @export pair.names pair names with overlap
#' @export names all names with overlap
Check_duplicate <- function(barcodes, return.pair.names = FALSE, return.names = FALSE){
        all_cells = Unlist(barcodes)
        
        dup_cells <- all_cells[duplicated(all_cells)]
        dup_cells <- dup_cells[sort(names(dup_cells))]
        dup_names <- unique(names(dup_cells))
        if(return.pair.names) {
                pair.names <- list()
                for(i in seq_along(dup_names)){
                        pair.names[[i]] = unique(names(all_cells)[all_cells %in% barcodes[[dup_names[i]]]])
                }
                return(pair.names)
        }
        
        if(return.names) {
                dup_cells = all_cells[all_cells %in% dup_cells]
                return(sort(unique(names(dup_cells))))
        }
}
Check_duplicate(barcodes = Annotations,return.pair.names = TRUE)
saveRDS(Annotations, file = "Yang/proximal_distal_terminal_COPD/Subset_Reclustering_by_markers/Annotations.rds")
Annotations = readRDS("Yang/proximal_distal_terminal_COPD/Subset_Reclustering_by_markers/Annotations.rds")
Annotations = Unlist(Annotations)
df = data.frame(annotations = names(Annotations), row.names = Annotations)
table(colnames(object) %in% Annotations)
unknown_cells <- colnames(object)[!(colnames(object) %in% Annotations)]
df_unknown = data.frame(annotations = rep("unknown",length(unknown_cells)), row.names = unknown_cells)
df %<>% rbind(df_unknown)
table(colnames(object) %in% rownames(df))
object[["annotations"]] = df
object %<>% AddMetaColor(label= "annotations", colors = Singler.colors)
saveRDS(object, "data/Lung_28_Global_20200511.rds")
object = readRDS(file = "data/Lung_28_Global_20200511.rds") 

#object@meta.data[unknown_cells,"annotations.colors"] = "#A9A9A9"
Idents(object) = "annotations"
lapply(c(T, F), function(label)
        UMAPPlot.1(object, group.by="annotations",pt.size = 0.5,label = label,
                   label.repel = T,alpha = 0.9,cols = Singler.colors,
                   no.legend = T,label.size = 4, repel = T, title = "Final annotation",
                   do.print = T, do.return = F))
sub_object <- subset(object,idents = "unknown", invert = T)

lapply(c(T, F), function(label)
        UMAPPlot.1(sub_object, group.by="annotations",pt.size = 0.5,label = label,
                    label.repel = T,alpha = 0.9,cols = Singler.colors,
                    no.legend = T,label.size = 4, repel = T, title = "Final annotation",
                    do.print = T, do.return = F))
rm(sub_object);GC()
meta.data = object@meta.data[,c("annotations","orig.ident")]
meta.data %<>% cbind(object@reductions$umap@cell.embeddings)
meta.data %<>% cbind(FetchData(object, vars = c("ACE2","TMPRSS2")))
colnames(meta.data)[1:2] = c("cell type","sample ID")
meta.data = meta.data[order(meta.data$TMPRSS2,decreasing = T),]
meta.data = meta.data[order(meta.data$ACE2,decreasing = T),]
write.csv(meta.data, file = paste0(path,"expression_ACE2_TMPRSS2.csv"))

table(object$orig.ident,object$annotations) %>% as.data.frame.matrix() %>%
        write.csv(file= paste0(path,"cell_type_per_sample.csv"))

meta.data = object@meta.data[,c("annotations","orig.ident")]
meta.data %<>% cbind(FetchData(object, vars = c("ACE2","TMPRSS2")))
percent <- function(Vector) sum(Vector >0)/length(Vector)
meta.data %<>% group_by(annotations,orig.ident) %>% 
        summarise(ACE2_mean_UMI = mean(ACE2),
                  TMPRSS2_mean_UMI = mean(TMPRSS2),
                  ACE2_Pct = paste(format(percent(ACE2)*100, digits=3), "%"),
                  TMPRSS2_Pct = paste(format(percent(TMPRSS2)*100,digits = 3), "%"))
colnames(meta.data)[1:2] = c("cell type","sample ID")
write.csv(meta.data, file = paste0(path,"Expression_summary_ACE2_TMPRSS2.csv"))
