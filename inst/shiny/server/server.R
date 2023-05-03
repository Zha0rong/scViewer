# Basic structure test:
## reative value: pre-determined name:

options(shiny.maxRequestSize=600000*1024^2)
library(SummarizedExperiment)
#source("../../R/import.R")




reactivevalue=reactiveValues(object_location=NULL,
                             Seurat_Object=NULL)

source('server/observer.R',local = T)








