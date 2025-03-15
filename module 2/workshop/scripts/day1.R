options(java.parameters = "-Xmx8g")
library(VoltRon)

####
# VoltRon object ####
####

# import data
data("merged_object")

# object
merged_object

# sample metadata
SampleMetadata(merged_object)

# metadata
Metadata(merged_object)

# main assay
vrMainAssay(merged_object) <- "MolAssay"

# molecule metadata
Metadata(merged_object)

# ROI metadata
Metadata(merged_object, assay = "ROIAssay")

# data matrix
vrMainAssay(merged_object) <- "CellAssay"
datax <- vrData(merged_object)

# features 
vrFeatures(merged_object)

####
# Images ####
####

# import data
data("melc_data")

# sample metadata
SampleMetadata(melc_data)

# Features
vrFeatures(melc_data)

# images and channels
vrImageChannelNames(melc_data)

# get images
vrImages(melc_data)
vrImages(melc_data, scale.perc = 20)
vrImages(melc_data, channel = "DAPI", scale.perc = 20)
vrImages(melc_data, channel = "CD45", scale.perc = 20)

# combine image channels
melc_data <- combineChannels(melc_data, channels = c("DAPI", "CD45"), colors = c("blue", "yellow"), channel_key = "combined")

# check new channels
vrImageChannelNames(melc_data)
vrImages(melc_data, channel = "combined", scale.perc = 20)

####
# Visualization ####
####

# import data
data("xenium_data")

# sample metadata
SampleMetadata(xenium_data)

# visualization
vrSpatialPlot(xenium_data, group.by = "clusters")
vrSpatialPlot(xenium_data, group.by = "clusters", plot.segments = TRUE)
vrSpatialPlot(xenium_data, group.by = "clusters", plot.segments = TRUE, alpha = 0.5)

# interactive visualization
vrSpatialPlot(xenium_data, group.by = "clusters", plot.segments = TRUE, alpha = 0.5, interactive = TRUE)

# import data
data("merged_object")
merged_object

# visualization of molecules
vrSpatialPlot(merged_object, assay = "MolAssay", group.by = "gene")
vrMainAssay(merged_object) <- "MolAssay"
vrSpatialPlot(merged_object, group.by = "gene")

# visualization of one molecule
vrSpatialPlot(merged_object, assay = "MolAssay", group.by = "gene")
vrMainAssay(merged_object) <- "MolAssay"
vrSpatialPlot(merged_object, group.by = "gene")
vrSpatialPlot(merged_object, group.by = "gene", group.ids = "KRT15")
vrSpatialPlot(merged_object, group.by = "gene", group.ids = "KRT14", colors = list(KRT14 = "purple"))

# visualization of one assay
vrSpatialPlot(merged_object, assay = "Assay2", group.by = "gene")

####
# Multilayer Visualization ####
####

# import data
data("merged_object")

# single plot
vrSpatialPlot(merged_object, group.by = "clusters", plot.segments = TRUE)

# combine plots
vrSpatialPlot(merged_object, group.by = "clusters", plot.segments = TRUE) |> 
  addSpatialLayer(merged_object, assay = "Assay2", group.by = "gene")

# combine molecules and cells
vrSpatialPlot(merged_object, group.by = "clusters", plot.segments = TRUE) |> 
  addSpatialLayer(merged_object, assay = "Assay2", group.by = "gene", colors = list(KRT14 = "purple", KRT15 = "yellow"), pt.size = 1)

# combine cells and annotations
Metadata(merged_object, assay = "ROIAssay")

vrSpatialPlot(merged_object, group.by = "clusters") |> 
  addSpatialLayer(merged_object, assay = "Assay3", group.by = "annotation")

vrSpatialPlot(merged_object, group.by = "clusters", background.color = "black") |> 
  addSpatialLayer(merged_object, assay = "Assay3", group.by = "annotation", alpha = 0.5)
