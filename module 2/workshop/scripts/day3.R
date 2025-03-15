options(java.parameters = "-Xmx8g")
library(VoltRon)
library(dplyr)

####
# Spatial Alignment ####
####

####
## Import Data ####
####

# Import analyzed data
Xen_R1 <- loadVoltRon("../data/ondisk/Xen_R1/")

####
## Spatial Alignment (same section) ####
####

# (OPTIONAL) import H&E data, and save to disk
Xen_R1_image <- importImageData("../data/BreastCancer/Xenium_R1/Xenium_FFPE_Human_Breast_Cancer_Rep1_he_image_highres.tif",
                                sample_name = "HEimage",
                                channel_names = "H&E", tile.size = 100)
Xen_R1_image_disk <- saveVoltRon(Xen_R1_image, 
                                 format = "HDF5VoltRon", 
                                 output = "../data/ondisk/Xen_R1_image/", 
                                 replace = TRUE)

# load H&E data from disk
Xen_R1_image_disk <- loadVoltRon("../data/ondisk/Xen_R1_image/")

# automated alignment
xen_reg <- registerSpatialData(object_list = list(Xen_R1, Xen_R1_image_disk))

# transfer image
Xenium_reg <- xen_reg$registered_spat[[2]]
vrImages(Xen_R1[["Assay1"]], channel = "H&E") <- vrImages(Xenium_reg, name = "main_reg", channel = "H&E")
vrImageChannelNames(Xen_R1)

# visualize plots
vrSpatialPlot(Xen_R1, group.by = "Clusters", channel = "H&E")

# manual alignment
xen_reg <- registerSpatialData(object_list = list(Xen_R1, Xen_R1_image_disk))

####
# Multiomics ####
####

####
## Import Omics Data ####
####

# get Xenium data
vr2_merged_acute1 <- readRDS(file = "../data/COVID19/acutecase1_annotated.rds")
vr2_merged_acute1
SampleMetadata(vr2_merged_acute1)

####
## Visualize ####
####

# visualize cells
vrSpatialPlot(vr2_merged_acute1, assay = "Xenium", group.by = "CellType",
              plot.segments = TRUE)

# visualize molecules
vrSpatialPlot(vr2_merged_acute1, assay = "Xenium_mol", group.by = "gene", 
              group.ids = c("S2_N", "S2_orf1ab"), n.tile = 500)

####
## Automated H&E Registration ####
####

# get H&E data
imgdata <- importImageData("../data/COVID19/acutecase1_heimage.jpg", 
                           tile.size = 10, 
                           segments = "../data/COVID19/acutecase1_membrane.geojson", 
                           sample_name = "acute case 1 (HE)", 
                           channel_names = "H&E")
imgdata <- flipCoordinates(imgdata, assay = "ROIAnnotation")
imgdata

# modulate image
vr2_merged_acute1 <- modulateImage(vr2_merged_acute1, brightness = 300)

# registration
# xen_reg <- registerSpatialData(object_list = list(vr2_merged_acute1, imgdata))

# non interactive registration
xen_reg <- registerSpatialData(object_list = list(vr2_merged_acute1, imgdata), 
                               mapping_parameters = readRDS("../data/COVID19/acutecase1_align_parameters.rds"), 
                               interactive = FALSE)

# spatial coordinate systems and image channels
vrSpatialNames(vr2_merged_acute1)
vrImageChannelNames(vr2_merged_acute1)

# transfer images
imgdata_reg <- xen_reg$registered_spat[[2]]
vrSpatialNames(imgdata_reg)
vrImages(vr2_merged_acute1[["Assay7"]], name = "main", channel = "H&E") <- 
  vrImages(imgdata_reg, assay = "Assay1", name = "main_reg")
vrImages(vr2_merged_acute1[["Assay8"]], name = "main", channel = "H&E") <- 
  vrImages(imgdata_reg, assay = "Assay1", name = "main_reg")

# new image channel names
vrImageChannelNames(vr2_merged_acute1)

# add ROI assay
vr2_merged_acute1 <- addAssay(vr2_merged_acute1,
                              assay = imgdata_reg[["Assay2"]],
                              metadata = Metadata(imgdata_reg, assay = "ROIAnnotation"),
                              assay_name = "ROIAnnotation",
                              sample = "acute case 1", layer = "Section1")
vrMainSpatial(vr2_merged_acute1[["Assay9"]]) <- "main_reg"
vr2_merged_acute1

# call regions membrane
vrMainAssay(vr2_merged_acute1) <- "ROIAnnotation"
vr2_merged_acute1$Region <- "Hyaline Membrane"

# visualize
vrSpatialPlot(vr2_merged_acute1, assay = "Xenium_mol", group.by = "gene", 
              group.ids = c("S2_N", "S2_orf1ab"), n.tile = 500) |>
  addSpatialLayer(vr2_merged_acute1, assay = "ROIAnnotation", 
                  group.by = "Region", alpha = 0.3, spatial = "main_reg",
                  colors = list(`Hyaline Membrane` = "blue")) 

####
## Label Transfer ####
####

# set the spatial coordinate system of ROI Annotations assay
vrSpatialNames(vr2_merged_acute1, assay = "all")
vrMainSpatial(vr2_merged_acute1[["Assay9"]]) <- "main_reg"

# transfer ROI annotations to molecules
vr2_merged_acute1 <- transferData(object = vr2_merged_acute1, 
                                  from = "Assay9", 
                                  to = "Assay8", 
                                  features = "Region")

# Metadata of molecules
Metadata(vr2_merged_acute1, assay = "Xenium_mol")

# check the abundance of molecules of hyaline membranes
s2_summary_hyaline <- 
  Metadata(vr2_merged_acute1, assay = "Xenium_mol") %>%
  filter(gene %in% c("S2_N", "S2_orf1ab"), 
         Region == "Hyaline Membrane") %>% 
  summarise(S2_N = sum(gene == "S2_N"), 
            S2_orf1ab = sum(gene == "S2_orf1ab"), 
            ratio = sum(gene == "S2_N")/sum(gene == "S2_orf1ab")) %>% 
  as.matrix()
s2_summary_hyaline

####
## Interactive Visualization ####
####

# basilisk: DONT RUN using rolv
if(!requireNamespace("basilisk"))
  BiocManager::install("basilisk")
library(basilisk)

# Dependencies
if(!requireNamespace("vitessceR"))
  devtools::install_github("vitessce/vitessceR")
library(vitessceR)

# write zarr
as.AnnData(vr2_merged_acute1, assay = "Assay7", file = "../data/ondisk/COVID19/vr2_merged_acute1_annotated_registered.zarr", 
           flip_coordinates = TRUE, create.ometiff = TRUE, name = "main", channel = "H&E")

# visualize
vrSpatialPlot("../data/ondisk/COVID19/vr2_merged_acute1_annotated_registered.zarr", group.by = "CellType")

# visualize for rolv
port <- httpuv::randomPort(min = 8000, max = 9000, n = 1000)
vrSpatialPlot("../data/ondisk/COVID19/vr2_merged_acute1_annotated_registered.zarr", group.by = "CellType", 
              shiny.options = list(host = "http://localhost", port = port))

# visualize docker
vrSpatialPlot(vr2_merged_acute1, channel = "H&E", group.by = "CellType", interactive = TRUE)

####
## Interactive Annotation ####
####

# annotate region
vr2_merged_acute1 <- annotateSpatialData(vr2_merged_acute1, assay = "Xenium_mol", 
                                         label = "Region2", use.image = TRUE, 
                                         channel = "H&E", annotation_assay = "InfectedAnnotation")

# Metadata
Metadata(vr2_merged_acute1, assay = "Xenium_mol")

# visualize multilayer
if(!requireNamespace("ggnewscale"))
  install.packages("ggnewscale")
library(ggnewscale)
vrSpatialPlot(vr2_merged_acute1, assay = "Xenium_mol", group.by = "gene", 
              group.ids = c("S2_N", "S2_orf1ab"), n.tile = 500) |>
  addSpatialLayer(vr2_merged_acute1, assay = "ROIAnnotation", group.by = "Region", alpha = 0.3, spatial = "main_reg",
                  colors = list(`Hyaline Membrane` = "blue")) |>
  addSpatialLayer(vr2_merged_acute1, assay = "InfectedAnnotation", group.by = "Region2", alpha = 0.3, 
                  colors = list(`Infected` = "yellow"))

# check abundance of molecules
s2_summary_infected <- 
  Metadata(vr2_merged_acute1, assay = "Xenium_mol") %>%
  filter(gene %in% c("S2_N", "S2_orf1ab"), 
         Region2 == "Infected") %>% 
  summarise(S2_N = sum(gene == "S2_N"), 
            S2_orf1ab = sum(gene == "S2_orf1ab"), 
            ratio = sum(gene == "S2_N")/sum(gene == "S2_orf1ab")) %>% 
  as.matrix()
s2_summary_infected
rbind(s2_summary_infected, s2_summary_hyaline)

# check all abundances
S2_table <- matrix(c(s2_summary_hyaline[,1:2], 
                     s2_summary_infected[,1:2]), 
                   dimnames = list(Region = c("Hyaline", "Infected"), 
                                   S2 = c("N", "orf1ab")),
                   ncol = 2,  byrow = TRUE)
S2_table
fisher.test(S2_table, alternative = "two.sided")

####
## Hot Spot Analysis ####
####

# get graph
vr2_merged_acute1 <- getSpatialNeighbors(vr2_merged_acute1, assay = "Xenium_mol", 
                                         group.by = "gene", group.ids = c("S2_N", "S2_orf1ab"), 
                                         method = "radius", radius = 50)

# get hotspot analysis
vr2_merged_acute1 <- getHotSpotAnalysis(vr2_merged_acute1, assay = "Xenium_mol", 
                                        features = "gene", graph.type = "radius")

# visualize
vrSpatialPlot(vr2_merged_acute1, assay = "Xenium_mol", 
              group.by = "gene_hotspot_flag", group.ids = c("cold", "hot"), 
              alpha = 1, background.color = "white")