options(java.parameters = "-Xmx8g")
library(VoltRon)

####
# Visium Analysis ####
####

####
## Import Data ####
####

# Dependencies
if(!requireNamespace("rhdf5"))
  BiocManager::install("rhdf5")
library(rhdf5)

# import Visium data
Ant_Sec1 <- importVisium("../data/Mouse Brain/Sagittal_Anterior/Section1/",
                         sample_name = "Anterior1")
Pos_Sec1 <- importVisium("../data/Mouse Brain/Sagittal_Posterior/Section1/",
                         sample_name = "Posterior1")

# merge datasets
MBrain_Sec <- merge(Ant_Sec1, Pos_Sec1)

####
## Omic Profile Clustering ####
####

####
### Processing and Filtering ####
####

# features
head(vrFeatures(MBrain_Sec))
length(vrFeatures(MBrain_Sec))

# normalize and select features
MBrain_Sec <- normalizeData(MBrain_Sec)
MBrain_Sec <- getFeatures(MBrain_Sec, n = 3000)

# selected features
head(vrFeatureData(MBrain_Sec))
selected_features <- getVariableFeatures(MBrain_Sec)
head(selected_features, 20)

####
### Dimensional Reduction ####
####

# embedding
MBrain_Sec <- getPCA(MBrain_Sec, features = selected_features, dims = 30)
MBrain_Sec <- getUMAP(MBrain_Sec, dims = 1:30)
vrEmbeddingNames(MBrain_Sec)

# embedding visualization
vrEmbeddingPlot(MBrain_Sec, embedding = "umap")

####
### Clustering ####
####

# graph for neighbors
MBrain_Sec <- getProfileNeighbors(MBrain_Sec, dims = 1:30, k = 10, method = "SNN")
vrGraphNames(MBrain_Sec)

# clustering
MBrain_Sec <- getClusters(MBrain_Sec, resolution = 0.5, label = "Clusters", graph = "SNN")

####
### Visualization ####
####

# embedding
vrEmbeddingPlot(MBrain_Sec, embedding = "umap", group.by = "Clusters")

####
## Niche Clustering ####
####

####
### Deconvolution ####
####

# libraries
if(!requireNamespace("Seurat"))
  install.packages("Seurat")
library(Seurat)

# import single cell data
allen_reference <- readRDS("../data/Mouse Brain/scRNA Mouse Brain/allen_cortex_analyzed_subset.rds")

# visualize
Idents(allen_reference) <- "subclass"
gsubclass <- DimPlot(allen_reference, reduction = "umap", label = T) + NoLegend()
Idents(allen_reference) <- "class"
gclass <- DimPlot(allen_reference, reduction = "umap", label = T) + NoLegend()
gsubclass | gclass

# Deconvolute anterior and posterior Visium spots
if(!requireNamespace("spacexr"))
  devtools::install_github("dmcable/spacexr")
library(spacexr)

# (OPTIONAL) deconvolute spots, but you can skip this as it is lenghty in time
MBrain_Sec <- getDeconvolution(MBrain_Sec, sc.object = allen_reference,
                               sc.cluster = "subclass", max_cores = 2)

# alternatively you can run the line below to get MBrain_Sec with deconvoluted spots
MBrain_Sec <- readRDS("../data/Mouse Brain/MBrain_Sec_decon.rds")

# Visualize
vrMainFeatureType(MBrain_Sec) <- "Decon"
vrFeatures(MBrain_Sec)
vrSpatialFeaturePlot(MBrain_Sec, features = c("L4", "L5 PT", "Oligo", "Vip"),
                     crop = TRUE, ncol = 2, alpha = 1, keep.scale = "all")

####
### Processing ####
####

vrMainFeatureType(MBrain_Sec) <- "Decon"
MBrain_Sec <- normalizeData(MBrain_Sec, method = "CLR")

# embedding visualization
MBrain_Sec <- getUMAP(MBrain_Sec, data.type = "norm", umap.key = "umap_niche")
vrEmbeddingPlot(MBrain_Sec, embedding = "umap_niche", group.by = "Sample")

####
### Niche Clustering ####
####

# clustering
MBrain_Sec <- getProfileNeighbors(MBrain_Sec, data.type = "norm", method = "SNN", graph.key = "SNN_niche")
MBrain_Sec <- getClusters(MBrain_Sec, resolution = 0.4, graph = "SNN_niche", label = "Niche_Clusters")

# visualize clustering
g1 <- vrEmbeddingPlot(MBrain_Sec, embedding = "umap", group.by = "Sample")
g2 <- vrEmbeddingPlot(MBrain_Sec, embedding = "umap", group.by = "Niche_Clusters", label = TRUE)
g1 | g2

####
### Visualization ####
####

# spatial clustering plot
vrSpatialPlot(MBrain_Sec, group.by = "Niche_Clusters", crop = TRUE, alpha = 1)

# heatmap plot
if(!requireNamespace("ComplexHeatmap"))
  BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)
vrHeatmapPlot(MBrain_Sec, features = vrFeatures(MBrain_Sec), group.by = "Niche_Clusters",
              show_row_names = T, show_heatmap_legend = T)

####
# Xenium Analysis ####
####

# Dependencies
if(!requireNamespace("rhdf5"))
  BiocManager::install("rhdf5")
if(!requireNamespace("RBioFormats"))
  BiocManager::install("RBioFormats")
library(rhdf5)

# import data
Xen_R1 <- importXenium("../data/BreastCancer/Xenium_R1/outs/", 
                       sample_name = "XeniumR1", 
                       overwrite_resolution = TRUE, 
                       resolution_level = 3)

####
## Filtering ####
####

# Metadata
Xen_R1@metadata
head(Metadata(Xen_R1))

# filter out counts
Xen_R1 <- subset(Xen_R1, Count > 5)

####
## Disk Representation ####
####

# Dependencies
if(!requireNamespace("rhdf5"))
  BiocManager::install("rhdf5")
if(!requireNamespace("HDF5Array"))
  BiocManager::install("HDF5Array")
if(!requireNamespace("ImageArray"))
  devtools::install_github("BIMSBbioinfo/ImageArray")
if(!requireNamespace("BPCells"))
  devtools::install_github("bnprks/BPCells/r")
library(rhdf5)
library(HDF5Array)
library(ImageArray)
library(BPCells)

# save voltron on disk
Xen_R1_ondisk <- saveVoltRon(Xen_R1, 
                             format = "HDF5VoltRon", 
                             output = "../data/ondisk/Xen_R1", 
                             replace = TRUE)

# load voltron from disk
Xen_R1_ondisk <- loadVoltRon("../data/ondisk/Xen_R1/")

####
## Omic Profile Clustering ####
####

# normalize
Xen_R1_ondisk <- normalizeData(Xen_R1_ondisk, sizefactor = 1000)

# PCA reduction
Xen_R1_ondisk <- getPCA(Xen_R1_ondisk, dims = 20, overwrite = TRUE)
Xen_R1_ondisk <- getUMAP(Xen_R1_ondisk, dims = 1:20)

# neighbors
Xen_R1_ondisk <- getProfileNeighbors(Xen_R1_ondisk, dims = 1:20, method = "SNN")
vrGraphNames(Xen_R1_ondisk)

# clustering
Xen_R1_ondisk <- getClusters(Xen_R1_ondisk, resolution = 1.3, label = "Clusters", graph = "SNN")

# visualization
vrEmbeddingPlot(Xen_R1_ondisk, group.by = "Clusters", embedding = "umap", 
                pt.size = 0.4, label = TRUE)

# spatial plot
vrSpatialPlot(Xen_R1_ondisk, group.by = "Clusters", pt.size = 0.18)

####
### Marker Analysis ####
####

library(Seurat)
Xen_R1$Clusters <- Xen_R1_ondisk$Clusters
Xen_R1_seu <- VoltRon::as.Seurat(Xen_R1, cell.assay = "Xenium", type = "image")
Idents(Xen_R1_seu) <- "Clusters"
Xen_R1_seu <- NormalizeData(Xen_R1_seu, scale.factor = 1000)
markers <- FindAllMarkers(Xen_R1_seu)

####
### Annotation ####
####

# get predefined annotations
annotations <- read.table("../data/BreastCancer/Xenium_R1/annotation.txt")[,1]

# annotate clusters
CellType <- factor(Xen_R1_ondisk$Clusters)
levels(CellType) <- annotations
Xen_R1_ondisk$CellType <- as.character(CellType)

# visualization
vrSpatialPlot(Xen_R1_ondisk, group.by = "CellType", pt.size = 0.18, alpha = 1)
vrEmbeddingPlot(Xen_R1_ondisk, group.by = "CellType", embedding = "umap",
                pt.size = 0.4, label = TRUE)

# save clusters and annotations to disk
Xen_R1_ondisk <- saveVoltRon(Xen_R1_ondisk)

####
## Spatially aware clustering ####
####

####
### Niche Clustering ####
####

# spatial neighbors
Xen_R1_ondisk <- getSpatialNeighbors(Xen_R1_ondisk, radius = 30, method = "radius")
vrGraphNames(Xen_R1_ondisk)

# get niche assay
Xen_R1_ondisk <- getNicheAssay(Xen_R1_ondisk, label = "CellType", graph.type = "radius")
Xen_R1_ondisk

# normalizing niche assay
vrMainFeatureType(Xen_R1_ondisk) <- "Niche"
Xen_R1_ondisk <- normalizeData(Xen_R1_ondisk, method = "CLR")

# clustering niches
Xen_R1_ondisk <- getClusters(Xen_R1_ondisk, nclus = 9, method = "kmeans", label = "niche_clusters")

# visualization
vrSpatialPlot(Xen_R1_ondisk, group.by = "niche_clusters", alpha = 1)
library(ComplexHeatmap)
vrHeatmapPlot(Xen_R1_ondisk, features = vrFeatures(Xen_R1_ondisk), group.by = "niche_clusters")

# visualization of specific cell type
vrSpatialPlot(Xen_R1_ondisk, group.by = "CellType", pt.size = 0.18, alpha = 1, group.ids = c("ACTA2_myoepithelial", "KRT15_myoepithelial"))
vrSpatialPlot(Xen_R1_ondisk, group.by = "CellType", pt.size = 1, alpha = 1, group.ids = c("CD4_TCells", "CD8_TCells", "BCells"), n.tile = 400)

####
### Hot Spot Analysis ####
####

# get spatial neighbor plot
Xen_R1_ondisk <- getSpatialNeighbors(Xen_R1_ondisk, method = "radius", radius = 15, graph.key = "radius_hot")

# visualize 
vrMainFeatureType(Xen_R1_ondisk) <- "RNA"
vrSpatialFeaturePlot(Xen_R1_ondisk, features = "PGR", alpha = 1, background.color = "black", n.tile = 300)

# analysis 
Xen_R1_ondisk <- getHotSpotAnalysis(Xen_R1_ondisk, features = "PGR", graph.type = "radius_hot", alpha.value = 0.001)

# visualize
vrSpatialFeaturePlot(Xen_R1_ondisk, features = "PGR_hotspot_stat", alpha = 1, background.color = "black", n.tile = 400)
vrSpatialPlot(Xen_R1_ondisk, group.by = "PGR_hotspot_flag", alpha = 1, background.color = "black", n.tile = 400)
