# cell type identification and validation 
library(Matrix)
library("reticulate")
library(Seurat)
library(dplyr)
library("stringr")


#---- read in gene expression matrix
setwd('~/Desktop/STARmap/36hpi_250genes_sum/')

exp.data <- Read10X(data.dir = "~/Desktop/STARmap/36hpi_250genes_sum/all_regions_1", gene.column = 1)
marker <- read.table("marker_long.txt", header = 1)
epi_marker <- read.tbale("epi_marker.txt", header = 1)
cell <- read.table("~/Desktop/STARmap/36hpi_250genes_sum/all_regions/barcodes.tsv", sep = "\t", row.names = 1)
colnames(cell) <- c( "x", "y", "z", "volume", "np", "np_nuc", "np_cyto", "pctnuc", "pctcyto", "region")

#----pre-processing

exp <- CreateSeuratObject(counts = exp.data, project = "all")
FeatureScatter(exp, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
VlnPlot(exp, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2, pt.size = 0)
exp <- AddMetaData(exp, cell)
VlnPlot(exp, features = c("volume"))
exp <- subset(exp, subset = nFeature_RNA >= 3 & 
                volume > 5e3 & volume <= 100000 & 
                nCount_RNA > 2 & nCount_RNA <= 500)

#---- get number of reads for marker genes only
expression_matrix <- GetAssayData(object = exp, assay = "RNA")[marker$marker, ]
exp[["marker_reads"]] <- colSums(expression_matrix) 

exp_sub <- subset(exp, subset = marker_reads > 1) 


#---- normalization 
exp_sub <- SCTransform(exp_sub)


#---- dimension reduction - only run on marker genes 
exp_sub <- RunPCA(exp_sub, features = marker$marker)
ElbowPlot(exp_sub, ndims = 50)


#---- find neighbors
exp_sub <- FindNeighbors(exp_sub, dims = 1:30 )
exp_sub <- FindClusters(exp_sub, resolution = 0.5, algorithm = 4)

exp_sub <- RunUMAP(exp_sub, dims = 1:30)
exp_sub <- RunTSNE(exp_sub, dims = 1:30, check_duplicates =FALSE)

DimPlot(exp_sub, reduction = "umap", label=TRUE, repel = TRUE)


#---- find marker genes for each cluster for 1st level clustering 
exp_sub.markers <- FindAllMarkers(exp_sub, features = marker$marker, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top_marker <- exp_sub.markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)


## for ambiguous clusters, run FindSubCluster for more pure population
# cluster 1
exp_sub = FindSubCluster(exp_sub, "1", graph.name = "SCT_snn", resolution = 0.2, subcluster.name = "lung", algorithm = 4)
exp_sub <- SetIdent(exp_sub, value = exp_sub@meta.data$lung)
exp_sub.markers <- FindAllMarkers(exp_sub, features = marker$marker, only.pos = TRUE, min.pct = 0.10, logfc.threshold = 0.25)
top <- exp.markers %>%
  group_by(cluster) %>%
  slice_max(n = 3, order_by = avg_log2FC)

# cluster 2
exp_sub = FindSubCluster(exp_sub, "2", graph.name = "SCT_snn", resolution = 0.2, subcluster.name = "lung", algorithm = 4)
exp_sub <- SetIdent(exp_sub, value = exp_sub@meta.data$lung)
exp_sub.markers <- FindAllMarkers(exp_sub, features = marker$marker, only.pos = TRUE, min.pct = 0.10, logfc.threshold = 0.25)
top <- exp.markers %>%
  group_by(cluster) %>%
  slice_max(n = 3, order_by = avg_log2FC)

# cluster 3
exp_sub = FindSubCluster(exp_sub, "3", graph.name = "SCT_snn", resolution = 0.2, subcluster.name = "lung", algorithm = 4)
exp_sub <- SetIdent(exp_sub, value = exp_sub@meta.data$lung)
exp_sub.markers <- FindAllMarkers(exp_sub, features = marker$marker, only.pos = TRUE, min.pct = 0.10, logfc.threshold = 0.25)
top <- exp.markers %>%
  group_by(cluster) %>%
  slice_max(n = 3, order_by = avg_log2FC)

# cluster 4
exp_sub = FindSubCluster(exp_sub, "4", graph.name = "SCT_snn", resolution = 0.2, subcluster.name = "lung", algorithm = 4)
exp_sub <- SetIdent(exp_sub, value = exp_sub@meta.data$lung)
exp_sub.markers <- FindAllMarkers(exp_sub, features = marker$marker, only.pos = TRUE, min.pct = 0.10, logfc.threshold = 0.25)
top <- exp.markers %>%
  group_by(cluster) %>%
  slice_max(n = 3, order_by = avg_log2FC)

# cluster 5
exp_sub = FindSubCluster(exp_sub, "5", graph.name = "SCT_snn", resolution = 0.2, subcluster.name = "lung", algorithm = 4)
exp_sub <- SetIdent(exp_sub, value = exp_sub@meta.data$lung)
exp_sub.markers <- FindAllMarkers(exp_sub, features = marker$marker, only.pos = TRUE, min.pct = 0.10, logfc.threshold = 0.25)
top <- exp.markers %>%
  group_by(cluster) %>%
  slice_max(n = 3, order_by = avg_log2FC)

# cluster 7
exp_sub = FindSubCluster(exp_sub, "7", graph.name = "SCT_snn", resolution = 0.2, subcluster.name = "lung", algorithm = 4)
exp_sub <- SetIdent(exp_sub, value = exp_sub@meta.data$lung)
exp_sub.markers <- FindAllMarkers(exp_sub, features = marker$marker, only.pos = TRUE, min.pct = 0.10, logfc.threshold = 0.25)
top <- exp.markers %>%
  group_by(cluster) %>%
  slice_max(n = 3, order_by = avg_log2FC)

# cluster 8
exp_sub = FindSubCluster(exp_sub, "8", graph.name = "SCT_snn", resolution = 0.2, subcluster.name = "lung", algorithm = 4)
exp_sub <- SetIdent(exp_sub, value = exp_sub@meta.data$lung)
exp_sub.markers <- FindAllMarkers(exp_sub, features = marker$marker, only.pos = TRUE, min.pct = 0.10, logfc.threshold = 0.25)
top <- exp.markers %>%
  group_by(cluster) %>%
  slice_max(n = 3, order_by = avg_log2FC)

exp_label <- RenameIdents(object = exp_sub,
                           "1_1" = "AT2", 
                           "1_3" = "AT2", 
                           "1_2_1" = "Tcell (CD4)",
                           "1_2_3" = "Tcell (CD4)",
                           "1_4_3" = "AlveolarFib",
                           "1_4_1" = "AlveolarFib",
                           "1_5_2" = "AlveolarFib",
                           "1_5_1" = "LymphaticEndo",
                           "1_5_4" = "pDC",
                           "1_6_2" = "Neutrophil",
                           "2_1" = "Capillary",
                           "2_2" = "InflamMac",
                           "2_5" = "DC",
                           "2_6" = "pDC",
                           "2_4" = "Artery",
                           "2_3" = "Vein",
                           "3_4" = "AdventitialFib", 
                           "3_3" = "AdventitialFib",
                           "3_5" = "SMC",
                           "3_6" = "SCMF",
                           "3_1" = "SMC",
                           "3_2" = "SMC",
                           "5_1" = "Pericyte",
                           "5_3" = "Pericyte",
                           "7_2" = "PatrolMac",
                           "7_4" = "PatrolMac",
                           "7_3" = "Tcell (CD8)",
                           "8_2" = "InterstitialMac",
                           "8_3" = "InterstitialMac",
                           "8_1" = "AlveolarMac",
                           "4_2" = "AirwayEpi",
                           "4_1_1" = "AirwayEpi",
                           "4_1_2" = "AirwayEpi",
                           "4_1_3" = "AirwayEpi",
                           "9" = "Bcell",
                           "6" = "AT1")


#---- subcluster airway epithelial cell cluster
epi_star = subset(exp_label, idents = "AirwayEpi")
epi_star <- SCTransform(epi_star)
epi_star <- RunPCA(epi_star, features = epi_marker$marker)

#---- find clusters within airway epithelia
epi_star <- FindNeighbors(epi_star, dims = 1:7 )
epi_star <- FindClusters(epi_star, resolution = 0.3, algorithm = 4)
DoHeatmap(epi_star, epi_marker$V1)
epi_star <- RunUMAP(epi_star, dims = 1:7)
DimPlot(epi_star, reduction = "umap", label=TRUE, repel = TRUE)


#---- determine markers for each cluster
epi_star.markers <- FindAllMarkers(epi_star, features = marker$marker)
FeaturePlot(epi, features = epi_marker$marker)
epi = FindSubCluster(epi, "1", graph.name = "SCT_snn", resolution = 0.2, subcluster.name = "lung", algorithm = 4)
epi <- SetIdent(epi, value = epi@meta.data$lung)
epi_label <- RenameIdents(object = epi,
                          "1_1" = "CLubEpi", 
                          "1_3" = "ClubEpi", 
                          "1_2" = "CiliatedEpi",
                          "2" = "CiliatedEpi",
                          "3" = "Unidentified",
                          "4" = "Unidentified",
                          "5" = "ClubEpi")
#---- map cells back to the main obj 
target_cells <- WhichCells(exp_label, idents = 'AirwayEpi')
Idents(exp_label, cells = target_cells) = epi_label

my_levels <- c("ClubEpi",
               "CiliatedEpi",
               "AT2",
               "AT1",
               'AdventitialFib',
               'AlveolarFib',
               'SMC',
               'SCMF',
               "Pericyte",
               "Artery",
               "Vein",
               "Capillary",
               "LymphaticEndo",
               "AlveolarMac",
               "InterstitialMac",
               "InflamMac",
               "PatrolMac",
               "DC",
               "pDC",
               "Neutrophil",
               "Tcell (CD4)",
               "Tcell (CD8)",
               "Bcell")

# Relevel object@ident
exp_label@active.ident <- factor(x = exp_label@active.ident, levels = my_levels)


#----average cluster average 
sc = readRDS('231006_day0_idents.RDS') # Seurat object for scRNA-seq
sc.average <- AggregateExpression(sc)
write.table(sc.average$RNA,
            '251111_sc_cluster_average.txt', sep = '\t', row.names = T, col.names = T, quote = F)
cluster.averages <- AggregateExpression(star)
write.table(cluster.average$RNA,
            '251111_spatial_cluster_average.txt', sep = '\t', row.names = T, col.names = T, quote = F)

#---- label transfer
# find anchor (scRNA-Seq as ref)
sc = readRDS('231006_day0_idents.RDS') # Seurat object for scRNA-seq

common_genes <- intersect(rownames(exp_label), rownames(sc))
anchors <- FindTransferAnchors(reference = sc, query = exp_label, 
                               dims = 1:30, 
                               features = common_genes, 
                               k.anchor = 5, k.filter = 500, k.score = 50)


# label transfer
predictions_cluster <- TransferData(anchorset = anchors, 
                                    refdata = Idents(sc_adult), 
                                    dims = 1:30)
exp_label@meta.data$predicted_cluster = predictions_cluster$predicted.id
exp_label@meta.data$predicted_cluster_score = predictions_cluster$prediction.score.max








