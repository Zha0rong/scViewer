
options(shiny.maxRequestSize=600000*1024^2)
library(SummarizedExperiment)



if (!require("HGNChelper")) install.packages("HGNChelper")
if (!require("openxlsx")) install.packages("openxlsx")

reactivevalue=reactiveValues(object_location=NULL,
                             function_input_object_location=NULL,
                             Seurat_Object=NULL,
                             Loaded=F,
                             Experiment_Metadata=NULL)

reactivevalue$object_location=filelocation

output$object_location=renderText(ifelse(is.null(reactivevalue$object_location),
                                         yes='No File Input, Use Upload data function.',
                                  no='The Data is being loaded.'))



source('server/observer.R',local = T)








