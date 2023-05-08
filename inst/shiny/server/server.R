
options(shiny.maxRequestSize=600000*1024^2)
library(SummarizedExperiment)




reactivevalue=reactiveValues(object_location=NULL,
                             function_input_object_location=NULL,
                             Seurat_Object=NULL)

filelocation <- getShinyOption("filelocation")
print(filelocation)

source('server/observer.R',local = T)








