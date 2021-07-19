#########################################################################################################################
# Simple example for creating input matrices for SPaRTAN analysis
#########################################################################################################################
dir = dirname(rstudioapi::getSourceEditorContext()$path)
setwd(dir)
rm(list=ls())
options(stringsAsFactors = FALSE)
library(Seurat)
data.dir <- "../data/"
##Load CITE-seq dataset
pbmc.data <- Read10X(data.dir = paste0(data.dir,"/cite-seq/pbmc5k_v3nextgem//filtered_feature_bc_matrix/"))
rownames(x = pbmc.data[["Antibody Capture"]]) <- gsub(pattern = "_[control_]*TotalSeqB", replacement = "", x = rownames(x = pbmc.data[["Antibody Capture"]]))
pbmc <- CreateSeuratObject(counts = pbmc.data[["Gene Expression"]], min.cells = ncol(pbmc.data[["Gene Expression"]])*0.03, min.features = 1000)
##For scRNA-seq
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA < 5000 & percent.mt < 30)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
##For scADT-seq
pbmc[["ADT"]] <- CreateAssayObject(pbmc.data[["Antibody Capture"]][setdiff(rownames(pbmc.data[["Antibody Capture"]]), c("IgG2a","IgG1","IgG2b")), colnames(x = pbmc)])
pbmc <- NormalizeData(pbmc, assay = "ADT", normalization.method = "CLR")
######################################################
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 5000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
rb.gene = grep('RPL|RPS|MRPS|MRPL', rownames(pbmc@assays$RNA@data), value = TRUE)
pbmc <- RunPCA(pbmc, features = setdiff(VariableFeatures(object = pbmc),rb.gene))
ElbowPlot(pbmc, ndims = 100)
pbmc <- FindNeighbors(pbmc, dims = 1:25)
pbmc <- FindClusters(pbmc, resolution = 0.7)
pbmc  <- RunTSNE(pbmc , dims = 1:25, method = "FIt-SNE")
DimPlot(pbmc , label = TRUE, reduction = "tsne", label.size = 4) + NoLegend()
###Marker genes for cell type assignment
# FeaturePlot(pbmc, features = c("CD8A","CD8B","adt_CD8a","adt_CD4","FCER1G"), label.size = 1)          #CD8+
# FeaturePlot(pbmc, features = c("IL7R", "CCR7","CD3E","adt_CD45RA","adt_CD4"))                         #Naive CD4+ T
# FeaturePlot(pbmc, features = c("IL7R", "S100A4","CD3E","adt_CD45RO","adt_CD4"))                       #Memory CD4+
# FeaturePlot(pbmc, features = c("GNLY", "NKG7", "adt_CD56"))                                           #NK
# FeaturePlot(pbmc, features = c("FCER1A", "CST3"))                                                     #DC
# FeaturePlot(pbmc, features = c("CD14", "LYZ", "adt_CD14"))                                            #CD14+ Mono
# FeaturePlot(pbmc, features = c("FCGR3A", "MS4A7", "adt_CD16"))                                        #FCGR3A+/CD16+ Mono
# FeaturePlot(pbmc, features = c("MS4A1", "adt_CD20", "adt_CD19"))                                      #B
new.cluster.ids <- c("NaiveCD4","CD14+Mono","MemoryCD4","CD14+Mono","NK","CD8","B","B","MemoryCD4", "DC", "CD16+Mono", "CD16+Mono","DC") ###5kn
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
pbmc[["classified"]] <- Idents(object = pbmc)
DimPlot(pbmc , label = TRUE, reduction = "tsne", label.size = 4) + NoLegend()
########Creating SPaRTAN input files
###Upload TF-target gene interaction matrix
D0 <- read.csv("../data/GenevsTF_SPaRTAN_sample.csv", row.names = 1)
common.genes <- intersect(rownames(pbmc),rownames(D0))
D <- D0[common.genes,]

P <- as.matrix(t(pbmc@assays$ADT@data)) ##10 8617
Y <- as.matrix(pbmc@assays$RNA@data[common.genes,]) ##20501 8617
P <- P[names(which(rowMeans(P)>0.5)),]
Y <- Y[,names(which(rowMeans(P)>0.5))]

celltypes <- as.character(pbmc@active.ident); names(celltypes) <- names(pbmc@active.ident)
B <- intersect(names(which(celltypes == "B")),colnames(Y))
cd14mono <- intersect(names(which(celltypes == "CD14+Mono")),colnames(Y))
cd16mono <- intersect(names(which(celltypes == "CD16+Mono")),colnames(Y))
cd4mem <- intersect(names(which(celltypes == "MemoryCD4")),colnames(Y))
cd4nav <- intersect(names(which(celltypes == "NaiveCD4")),colnames(Y))
cd8 <- intersect(names(which(celltypes == "CD8")),colnames(Y))
dc <- intersect(names(which(celltypes == "DC")),colnames(Y))
NK <- intersect(names(which(celltypes == "NK")),colnames(Y))
#########
input_dir <- ""
write.csv(Y[rownames(D),B], file=paste0(input_dir,"Y_B.csv")); write.csv(P[B,], file=paste0(input_dir,"P_B.csv"))
write.csv(Y[rownames(D),cd4nav], file=paste0(input_dir,"Y_CD4nav.csv")); write.csv(P[cd4nav,], file=paste0(input_dir,"P_CD4nav.csv"))
write.csv(Y[rownames(D),cd4mem], file=paste0(input_dir,"Y_CD4mem.csv")); write.csv(P[cd4mem,], file=paste0(input_dir,"P_CD4mem.csv"))
write.csv(Y[rownames(D),cd8], file=paste0(input_dir,"Y_CD8.csv")); write.csv(P[cd8,], file=paste0(input_dir,"P_CD8.csv"))
write.csv(Y[rownames(D),NK], file=paste0(input_dir,"Y_NK.csv")); write.csv(P[NK,], file=paste0(input_dir,"P_NK.csv"))
write.csv(Y[rownames(D),cd16mono], file=paste0(input_dir,"Y_CD16+MONO.csv")); write.csv(P[cd16mono,], file=paste0(input_dir,"P_CD16+MONO.csv"))
write.csv(Y[rownames(D),cd14mono], file=paste0(input_dir,"Y_CD14+MONO.csv")); write.csv(P[cd14mono,], file=paste0(input_dir,"P_CD14+MONO.csv"))
write.csv(Y[rownames(D),dc], file=paste0(input_dir,"Y_DC.csv")); write.csv(P[dc,], file=paste0(input_dir,"P_DC.csv"))
