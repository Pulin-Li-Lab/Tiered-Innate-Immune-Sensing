# cell type identification and validation 
library(Matrix)
library("reticulate")
library(Seurat)
library(dplyr)
library("stringr")

## 1 - Process STARmap data 
#---- read in gene scression matrix
setwd('~/Desktop/STARmap/36hpi_250genes_sum/')

sc.data <- Read10X(data.dir = "~/Desktop/STARmap/36hpi_250genes_sum/all_regions_1", gene.column = 1)
marker <- read.table("marker_long.txt", header = 1)
epi_marker <- read.table("epi_marker.txt", header = 1)
cell <- read.table("~/Desktop/STARmap/36hpi_250genes_sum/all_regions/barcodes.tsv", sep = "\t", row.names = 1)
colnames(cell) <- c( "x", "y", "z", "volume", "np", "np_nuc", "np_cyto", "pctnuc", "pctcyto", "region")

#----pre-processing

sc <- CreateSeuratObject(counts = sc.data, project = "all")
FeatureScatter(sc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
VlnPlot(sc, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2, pt.size = 0)
sc <- AddMetaData(sc, cell)
VlnPlot(sc, features = c("volume"))
sc <- subset(sc, subset = nFeature_RNA >= 3 & 
                volume > 5e3 & volume <= 100000 & 
                nCount_RNA > 2 & nCount_RNA <= 500)

#---- get number of reads for marker genes only
scression_matrix <- GetAssayData(object = sc, assay = "RNA")[marker$marker, ]
sc[["marker_reads"]] <- colSums(scression_matrix) 

sc_sub <- subset(sc, subset = marker_reads > 1) 


#---- normalization 
sc_sub <- SCTransform(sc_sub)


#---- dimension reduction - only run on marker genes 
sc_sub <- RunPCA(sc_sub, features = marker$marker)
ElbowPlot(sc_sub, ndims = 50)


#---- find neighbors
sc_sub <- FindNeighbors(sc_sub, dims = 1:30 )
sc_sub <- FindClusters(sc_sub, resolution = 0.5, algorithm = 4)

sc_sub <- RunUMAP(sc_sub, dims = 1:30)
sc_sub <- RunTSNE(sc_sub, dims = 1:30, check_duplicates =FALSE)

DimPlot(sc_sub, reduction = "umap", label=TRUE, repel = TRUE)


#---- find marker genes for each cluster for 1st level clustering 
sc_sub.markers <- FindAllMarkers(sc_sub, features = marker$marker, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top_marker <- sc_sub.markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)


## for ambiguous clusters, run FindSubCluster for more pure population
# cluster 1
sc_sub = FindSubCluster(sc_sub, "1", graph.name = "SCT_snn", resolution = 0.2, subcluster.name = "lung", algorithm = 4)
sc_sub <- SetIdent(sc_sub, value = sc_sub@meta.data$lung)
sc_sub.markers <- FindAllMarkers(sc_sub, features = marker$marker, only.pos = TRUE, min.pct = 0.10, logfc.threshold = 0.25)
top <- sc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 3, order_by = avg_log2FC)

# cluster 2
sc_sub = FindSubCluster(sc_sub, "2", graph.name = "SCT_snn", resolution = 0.2, subcluster.name = "lung", algorithm = 4)
sc_sub <- SetIdent(sc_sub, value = sc_sub@meta.data$lung)
sc_sub.markers <- FindAllMarkers(sc_sub, features = marker$marker, only.pos = TRUE, min.pct = 0.10, logfc.threshold = 0.25)
top <- sc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 3, order_by = avg_log2FC)

# cluster 3
sc_sub = FindSubCluster(sc_sub, "3", graph.name = "SCT_snn", resolution = 0.2, subcluster.name = "lung", algorithm = 4)
sc_sub <- SetIdent(sc_sub, value = sc_sub@meta.data$lung)
sc_sub.markers <- FindAllMarkers(sc_sub, features = marker$marker, only.pos = TRUE, min.pct = 0.10, logfc.threshold = 0.25)
top <- sc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 3, order_by = avg_log2FC)

# cluster 4
sc_sub = FindSubCluster(sc_sub, "4", graph.name = "SCT_snn", resolution = 0.2, subcluster.name = "lung", algorithm = 4)
sc_sub <- SetIdent(sc_sub, value = sc_sub@meta.data$lung)
sc_sub.markers <- FindAllMarkers(sc_sub, features = marker$marker, only.pos = TRUE, min.pct = 0.10, logfc.threshold = 0.25)
top <- sc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 3, order_by = avg_log2FC)

# cluster 5
sc_sub = FindSubCluster(sc_sub, "5", graph.name = "SCT_snn", resolution = 0.2, subcluster.name = "lung", algorithm = 4)
sc_sub <- SetIdent(sc_sub, value = sc_sub@meta.data$lung)
sc_sub.markers <- FindAllMarkers(sc_sub, features = marker$marker, only.pos = TRUE, min.pct = 0.10, logfc.threshold = 0.25)
top <- sc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 3, order_by = avg_log2FC)

# cluster 7
sc_sub = FindSubCluster(sc_sub, "7", graph.name = "SCT_snn", resolution = 0.2, subcluster.name = "lung", algorithm = 4)
sc_sub <- SetIdent(sc_sub, value = sc_sub@meta.data$lung)
sc_sub.markers <- FindAllMarkers(sc_sub, features = marker$marker, only.pos = TRUE, min.pct = 0.10, logfc.threshold = 0.25)
top <- sc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 3, order_by = avg_log2FC)

# cluster 8
sc_sub = FindSubCluster(sc_sub, "8", graph.name = "SCT_snn", resolution = 0.2, subcluster.name = "lung", algorithm = 4)
sc_sub <- SetIdent(sc_sub, value = sc_sub@meta.data$lung)
sc_sub.markers <- FindAllMarkers(sc_sub, features = marker$marker, only.pos = TRUE, min.pct = 0.10, logfc.threshold = 0.25)
top <- sc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 3, order_by = avg_log2FC)

sc_label <- RenameIdents(object = sc_sub,
                           "1_1" = "AT2", 
                           "1_3" = "AT2", 
                           "1_2_1" = "Tcell (CD4)",
                           "1_2_3" = "Tcell (CD4)",
                           "1_4_3" = "Col13Fib",
                           "1_4_1" = "Col13Fib",
                           "1_5_2" = "Col13Fib",
                           "1_5_1" = "LymphaticEndo",
                           "1_5_4" = "pDC",
                           "1_6_2" = "Neutrophil",
                           "2_1" = "Capillary",
                           "2_2" = "InflamMac",
                           "2_5" = "DC",
                           "2_6" = "pDC",
                           "2_4" = "Artery",
                           "2_3" = "Vein",
                           "3_4" = "Col14Fib", 
                           "3_3" = "Col14Fib",
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
epi_star = subset(sc_label, idents = "AirwayEpi")
epi_star <- SCTransform(epi_star)
epi_star <- RunPCA(epi_star, features = epi_marker$marker)

#---- find clusters within airway epithelia
epi_star <- FindNeighbors(epi_star, dims = 1:7 )
epi_star <- FindClusters(epi_star, resolution = 0.3, algorithm = 4)
DoHeatmap(epi_star, epi_marker$V1)
epi_star <- RunUMAP(epi_star, dims = 1:7)
DimPlot(epi_star, reduction = "umap", label=TRUE, repel = TRUE)


#---- determine markers for each cluster
epi_star.markers <- FindAllMarkers(epi_star, features = epi_marker$marker)
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
target_cells <- WhichCells(sc_label, idents = 'AirwayEpi')
Idents(sc_label, cells = target_cells) = epi_label

my_levels <- c("ClubEpi",
               "CiliatedEpi",
               "AT2",
               "AT1",
               'Col14Fib',
               'Col13Fib',
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
sc_label@active.ident <- factor(x = sc_label@active.ident, levels = my_levels)


## 2 - Process single-cell RNA-seq data 
###
sc.data <- Read10X(data.dir = "~/Desktop/STARmap/SingleCell_Lung/Cell", gene.column = 1)
# pre-processing
sc <- CreateSeuratObject(counts = sc.data, project = "all", min.cell = 5, min.features = 200)
FeatureScatter(sc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
sc[["percent.mt"]] <- PercentageFeatureSet(sc, pattern = "^mt-")
VlnPlot(sc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
sc <- subset(sc, subset = nFeature_RNA >= 200 & nFeature_RNA <= 2500 & percent.mt < 5)


##  Quality control
meta_object <- sc
meta_object$log10GenesPerUMI <- log10(meta_object$nFeature_RNA) / log10(meta_object$nCount_RNA)
meta_object$MitoRatio <- meta_object$percent.mt / 100
meta <- meta_object@meta.data
meta$cells <- colnames(meta_object)
meta <- meta %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)
meta$sample <- meta$seq_folder

library(ggplot2)
# plot number of cells
meta %>% 
  ggplot(aes(x=sample, fill=sample)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Number of Cells")
# UMI counts per cell -- should be above 1000 for good sequencing depth 
meta %>% 
  ggplot(aes(color=sample, x=nUMI, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 1000)

# number of genes per cell 
# Visualize the distribution of genes detected per cell via histogram
meta %>% 
  ggplot(aes(color=sample, x=nGene, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 300)

# Visualize the distribution of genes detected per cell via boxplot
meta %>% 
  ggplot(aes(x=sample, y=log10(nGene), fill=sample)) + 
  geom_boxplot() + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells vs NGenes")

# ratio of mitochondrial genes
meta %>% 
  ggplot(aes(color=sample, x=MitoRatio, fill=sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 0.2)

# complexity 
meta %>%
  ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8)


meta %>% 
  ggplot(aes(x=nUMI, y=nGene, color=MitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 1000) +
  geom_hline(yintercept = 500) +
  facet_wrap(~sample)

meta_object@meta.data <- meta
filtered_seurat <- subset(meta_object, 
                          subset= (nUMI >= 1000) & 
                            (nGene >= 500) & 
                            (log10GenesPerUMI > 0.80) & 
                            (MitoRatio < 0.20)) 


# normalization
sc <- filtered_seurat
# sc <- NormalizeData(sc, normalization.method = "LogNormalize", scale.factor = 1e4)
# sc <- FindVariableFeatures(sc, selection.method = "vst", nfeatures = 2000)
# all.genes <- rownames(sc)
# sc <- ScaleData(sc, features = all.genes)
sc <- SCTransform(sc)

# Dimension reduction
marker <- read.table("marker_long.txt")
sc <- RunPCA(sc)
DimPlot(sc, reduction = "pca")
DimHeatmap(sc, dims = 1:5, cells = 500, balanced = TRUE)
ElbowPlot(sc, ndims= 50)

# cluster cells
sc <- FindNeighbors(sc, dims = 1:35)
sc <- FindClusters(sc, resolution = 0.5)

# umap
sc <- RunUMAP(sc, dims = 1:35)
DimPlot(sc, reduction = "umap", label = TRUE)

#find markers
sc.markers <- FindAllMarkers(sc, features = VariableFeatures(sc), only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top_markers <- sc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)

sc <- FindSubCluster(sc, "21", graph.name = "RNA_snn", resolution = 0.2, subcluster.name = "lung")

sc <- SetIdent(sc, value = sc@meta.data$lung)

sc <- RenameIdents(object = sc,
                    "0" = "AlveolarMac",
                    "1" = "Artery",
                    "2" = "Capillary" ,
                    "3" = "Col13+FB",
                    "4" = "InflammatoryMac",
                    "5" = "T cell" ,
                    "6" = "PatrollingMac", 
                    "7" = "Pericyte",
                    "8" = "Vein",
                    "9" = "B cell",
                    "10" = "CD11b+DC", 
                    "11" = "AT2",
                    "12" = "Germinal B cell",
                    "13" = "CD103+ DC",
                    "14" = "Col14+FB", 
                    "15" = "Th2",
                    "16" = "Monocyte",
                    "17" = "Neutrophil",
                    "18" = "AT1",
                    "19" = "InterstitialMac",
                    "20" = "MyoFB" ,
                    "23" = "Mesothelium",
                    "24" = "pDC", 
                    "25" = "CiliatedEpi",
                    "26" = "Lymphatic", 
                    "27" = "Club") 
sc <- RenameIdents(object = sc,
                          "21_0" = "MyoFB", 
                          "21_1" = "MyoFB")
sc <- subset(sc, idents = c("21", "22"), invert = TRUE)

## 3----average cluster average 
sc.average <- AverageExpression(sc, assay = 'SCT', slot = 'scale.data')
write.table(sc.average$SCT,
            '251111_sc_cluster_average.txt', sep = '\t', row.names = T, col.names = T, quote = F)
cluster.averages <- AverageExpression(star, assay = 'SCT', slot = 'scale.data')
write.table(cluster.averages$SCT,
            '251111_spatial_cluster_average.txt', sep = '\t', row.names = T, col.names = T, quote = F)

## 4---- label transfer
common_genes <- intersect(rownames(star_label), rownames(sc))
anchors <- FindTransferAnchors(reference = sc, query = star_label, 
                               dims = 1:30, 
                               features = common_genes,
                               reference.reduction = "pca")

# label transfer
predictions_cluster <- TransferData(anchorset = anchors, 
                                    refdata = Idents(sc), 
                                    dims = 1:30)
pred_id = predictions_cluster$predicted.id
pred_max = predictions_cluster$prediction.score.max

# Keep only high-confidence predictions (>0.5)
pred_id_filtered <- ifelse(pred_max > 0.5, pred_id, "Unidentified")

# Store in metadata
star_label$predicted_cluster       <- pred_id_filtered
star_label$predicted_cluster_score <- pred_max








