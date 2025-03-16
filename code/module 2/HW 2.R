options(java.parameters = "-Xmx8g")
library(VoltRon)

vis <- importVisium("module 2/workshop/data/BreastCancer/Visium/")
vrSpatialPlot(vis, crop = FALSE, alpha = 1)

Xen_R1 <- loadVoltRon("module 2/workshop/data/ondisk/Xen_R1/")
vrSpatialPlot(Xen_R1, pt.size = 0.18)

vis_regis<-registerSpatialData(object_list = list(vis, Xen_R1))
registerSpatialData(object_list = list(vis, Xen_R1))

vrSpatialPlot(vrSpatialPlot)
