
options(shiny.maxRequestSize=600000*1024^2)
library(SummarizedExperiment)




reactivevalue=reactiveValues(object_location=NULL,
                             function_input_object_location=NULL,
                             Seurat_Object=NULL)

reactivevalue$object_location=filelocation

print(filelocation)
output$object_location=renderText(ifelse(is.null(reactivevalue$object_location),
                                         yes='No File Input, Use Upload data function.',
                                  no='File already uploaded, loading now.'))

observe({if (!is.null(reactivevalue$object_location)) {
  disable("Seurat_Object")
  disable("submit")










}})


source('server/observer.R',local = T)








