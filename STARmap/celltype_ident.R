# cell type identification and validation 
library(Matrix)
library("reticulate")
library(Seurat)
library(dplyr)
library("stringr")


# read in gene expression matrix
setwd('~/Desktop/STARmap/')

exp.data <- Read10X(data.dir = "~/Desktop/STARmap/36hpi_250genes_sum/all_regions", gene.column = 1)
marker <- read.table("marker_long.txt", header = 1)
cell <- read.table("~/Desktop/STARmap/36hpi_250genes_sum/all_regions/barcodes.tsv", sep = "\t", row.names = 1)
colnames(cell) <- c( "x", "y", "z", "volume", "np", "np_nuc", "np_cyto", "pctnuc", "pctcyto", "region")


## pre-processing
exp <- CreateSeuratObject(counts = exp.data, project = "all")
FeatureScatter(exp, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
VlnPlot(exp, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2, pt.size = 0)
exp <- AddMetaData(exp, cell)
VlnPlot(exp, features = c("volume"))
exp <- subset(exp, subset = nFeature_RNA >= 3 & 
                volume > 5e3 & volume <= 100000 & 
                nCount_RNA > 5 & nCount_RNA <= 500)

# get number of reads for marker genes only
expression_matrix <- GetAssayData(object = exp, assay = "RNA")[marker$marker, ]
exp[["marker_reads"]] <- colSums(expression_matrix) 

exp_sub <- subset(exp, subset = marker_reads > 1) 

# normalization 
exp_sub <- SCTransform(exp_sub)

## dimension reduction - only run on marker genes 
exp_sub <- RunPCA(exp_sub, features = marker$marker)
ElbowPlot(exp_sub, ndims = 50)

# find neighbors
exp_sub <- FindNeighbors(exp_sub, dims = 1:30 )
exp_sub <- FindClusters(exp_sub, resolution = 0.5, algorithm = 4)

exp_sub <- RunUMAP(exp_sub, dims = 1:30)
exp_sub <- RunTSNE(exp_sub, dims = 1:30, check_duplicates =FALSE)

DimPlot(exp_sub, reduction = "umap", label=TRUE, repel = TRUE)

# find marker genes for each cluster
exp_sub.markers <- FindAllMarkers(exp_sub, features = marker$marker, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

top_marker <- exp_sub.markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)

# for each identified cluster, run FindSubCluster for more pure population
exp_sub = FindSubCluster(exp_sub, "3", graph.name = "SCT_snn", resolution = 0.2, subcluster.name = "lung", algorithm = 4)

exp_sub <- SetIdent(exp_sub, value = exp_sub@meta.data$lung)

exp_sub.markers <- FindAllMarkers(exp_sub, features = marker$marker, only.pos = TRUE, min.pct = 0.10, logfc.threshold = 0.25)

top <- exp.markers %>%
  group_by(cluster) %>%
  slice_max(n = 3, order_by = avg_log2FC)


## subcluster airway epithelial cell cluster
epi_star = subset(exp_label, idents = "AirwayEpi")

epi_star <- SCTransform(epi_star)
epi_star <- RunPCA(epi_star, features = epi_marker$marker)

# find clusters within airway epithelia
epi_star <- FindNeighbors(epi_star, dims = 1:7 )
epi_star <- FindClusters(epi_star, resolution = 0.3, algorithm = 4)
DoHeatmap(epi_star, epi_marker$V1)
epi_star <- RunUMAP(epi_star, dims = 1:7)
DimPlot(epi_star, reduction = "umap", label=TRUE, repel = TRUE)


# determine markers for each cluster
epi_star.markers <- FindAllMarkers(epi_star, features = marker$marker)
FeaturePlot(epi, features = epi_marker$marker)
epi = FindSubCluster(epi, "1", graph.name = "SCT_snn", resolution = 0.2, subcluster.name = "lung", algorithm = 4)
epi <- SetIdent(epi, value = epi@meta.data$lung)

# map cells back to the main obj 
target_cells <- WhichCells(exp_label, idents = 'AirwayEpi')
Idents(exp_label, cells = target_cells) = epi_ident


## label transfer
# find anchor (scRNA-Seq as ref)
sc = readRDS('lungmap.rds')

common_genes <- intersect(rownames(exp), rownames(sc))
anchors <- FindTransferAnchors(reference = sc, query = exp_label, 
                               dims = 1:30, reduction = 'cca', 
                               features = common_genes, 
                               k.anchor = 5, k.filter = 500, k.score = 50)


# label transfer
predictions_cluster <- TransferData(anchorset = anchors, 
                                    refdata = sc_adult$celltype_level3, 
                                    dims = 1:n_pcs)

exp_label@meta.data$predicted_cluster = predictions_cluster$predicted.id
exp_label@meta.data$predicted_cluster_score = predictions_cluster$prediction.score.max



## generate statistics of infection and IFN+ cells per region 

cell_type_counts <- table(Idents(exp_label), exp_label$region)
cell_type_percentages <- prop.table(cell_type_counts, 2) * 100
write.csv(as.data.frame.matrix(cell_type_percentages), 'AllRegions_AllCell.csv')

infect <- subset(exp_label, subset = NP.gRNA > 0, slot = 'counts')
infect_counts <- table(Idents(infect), infect$region)
infect_percentages <- prop.table(infect_counts, 2) * 100
write.csv(as.data.frame.matrix(infect_percentages), 'AllRegions_AllInfect.csv')

ifn <- subset(infect, subset = Ifnb1 >0 | Ifnl2 > 0 | Ifna2 > 0 | Ifna4 > 0, slot = 'counts')
ifn_counts <- table(Idents(ifn), ifn$region)
ifn_percentages <- prop.table(ifn_counts, 2) * 100
write.csv(as.data.frame.matrix(ifn_percentages), 'AllRegions_AllIFN.csv')






