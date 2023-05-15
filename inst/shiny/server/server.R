
options(shiny.maxRequestSize=600000*1024^2)
library(SummarizedExperiment)




reactivevalue=reactiveValues(object_location=NULL,
                             function_input_object_location=NULL,
                             Seurat_Object=NULL)

reactivevalue$object_location=filelocation

output$object_location=renderText(ifelse(is.null(reactivevalue$object_location),
                                         yes='No File Input, Use Upload data function.',
                                  no='File already uploaded, loading now.'))

observe({if (!is.null(reactivevalue$object_location)) {
  disable("Seurat_Object")
  disable("submit")
  if (!is.null(reactivevalue$object_location)){
    withProgress(message = 'Load Data Object',value = 0, {

      n=8
      reactivevalue$SeuratObject=readRDS(reactivevalue$object_location)
      incProgress(1/n,detail = 'Start Loading')


      tryCatch({reactivevalue$SeuratObject=as.Seurat(reactivevalue$SeuratObject)},
               error=function(cond) {
                 message('This is not a SingleCellExperiment Object')
               }
      )

      incProgress(1/n,detail = 'Finish Loading')

      reactivevalue$SeuratObject$ShinyGroup='SCViewer'
      reactivevalue$Experiment_Metadata=reactivevalue$SeuratObject@meta.data
      DefaultAssay(reactivevalue$SeuratObject)='RNA'
      output$MainFigure=renderPlot(DimPlot(reactivevalue$SeuratObject))
      updateSelectizeInput(session = session,inputId = 'reduction',choices =names((reactivevalue$SeuratObject@reductions)),selected = NULL)
      updateSelectizeInput(session = session,inputId = 'GenesToInterrogate',choices =rownames(reactivevalue$SeuratObject),selected = NULL,server=T)
      incProgress(1/n,detail = 'Finish parsing Genes')

      updateSelectizeInput(session = session,inputId = 'variabletogroup',choices=colnames(reactivevalue$Experiment_Metadata)
                           [
                             !grepl("nCount",colnames(reactivevalue$Experiment_Metadata))&!grepl("nFeature",colnames(reactivevalue$Experiment_Metadata))&!grepl("^percent.",colnames(reactivevalue$Experiment_Metadata))
                             #&!grepl("_res.",colnames(reactivevalue$Experiment_Metadata))
                           ]
                           ,selected = 'ShinyGroup')

      incProgress(1/n,detail = 'Finish Parsing Columns step 1')

      updateSelectizeInput(session = session,inputId = 'variabletosplit',choices=colnames(reactivevalue$Experiment_Metadata)
                           [
                             !grepl("nCount",colnames(reactivevalue$Experiment_Metadata))&!grepl("nFeature",colnames(reactivevalue$Experiment_Metadata))&!grepl("^percent.",colnames(reactivevalue$Experiment_Metadata))
                             #  &!grepl("_res.",colnames(reactivevalue$Experiment_Metadata))
                           ]
                           ,selected = 'ShinyGroup')
      incProgress(1/n,detail = 'Finish Parsing Columns step 2')

      updateSelectizeInput(session = session,inputId = 'SampleColumn',choices=colnames(reactivevalue$Experiment_Metadata)
                           [
                             !grepl("nCount",colnames(reactivevalue$Experiment_Metadata))&!grepl("nFeature",colnames(reactivevalue$Experiment_Metadata))&!grepl("^percent.",colnames(reactivevalue$Experiment_Metadata))
                             #  &!grepl("_res.",colnames(reactivevalue$Experiment_Metadata))
                           ]
                           ,selected = NULL)
      incProgress(1/n,detail = 'Finish Parsing Columns step 3')

      updateSelectizeInput(session = session,inputId = 'PlotGroup',choices=colnames(reactivevalue$Experiment_Metadata)
                           [
                             !grepl("nCount",colnames(reactivevalue$Experiment_Metadata))&!grepl("nFeature",colnames(reactivevalue$Experiment_Metadata))&!grepl("^percent.",colnames(reactivevalue$Experiment_Metadata))
                             #  &!grepl("_res.",colnames(reactivevalue$Experiment_Metadata))
                           ]
                           ,selected = NULL)

      incProgress(1/n,detail = 'Finish Parsing Columns step 4')

      updateSelectizeInput(session = session,inputId = 'BarGraph1',choices=colnames(reactivevalue$Experiment_Metadata)
                           [
                             !grepl("nCount",colnames(reactivevalue$Experiment_Metadata))&!grepl("nFeature",colnames(reactivevalue$Experiment_Metadata))&!grepl("^percent.",colnames(reactivevalue$Experiment_Metadata))
                             #  &!grepl("_res.",colnames(reactivevalue$Experiment_Metadata))
                           ]
                           ,selected = NULL)

      updateSelectizeInput(session = session,inputId = 'BarGraph2',choices=colnames(reactivevalue$Experiment_Metadata)
                           [
                             !grepl("nCount",colnames(reactivevalue$Experiment_Metadata))&!grepl("nFeature",colnames(reactivevalue$Experiment_Metadata))&!grepl("^percent.",colnames(reactivevalue$Experiment_Metadata))
                             #  &!grepl("_res.",colnames(reactivevalue$Experiment_Metadata))
                           ]
                           ,selected = NULL)

      DGE_Group_Candidate=c()
      for (i in colnames(reactivevalue$Experiment_Metadata)[
        !grepl("nCount",colnames(reactivevalue$Experiment_Metadata))&!grepl("nFeature",colnames(reactivevalue$Experiment_Metadata))&!grepl("^percent.",colnames(reactivevalue$Experiment_Metadata))
        &!grepl("_res.",colnames(reactivevalue$Experiment_Metadata))
      ]) {
        DGE_Group_Candidate=c(DGE_Group_Candidate,paste0(i,'-',unique(as.character(reactivevalue$Experiment_Metadata[,i]))))
      }
      updateSelectizeInput(session = session,inputId = 'Assay',choices=names(reactivevalue$SeuratObject@assays),selected = NULL)
      updateSelectizeInput(session = session,inputId = 'DGEGroup1',choices=DGE_Group_Candidate,selected = NULL)
      updateSelectizeInput(session = session,inputId = 'DGEGroup2',choices=DGE_Group_Candidate,selected = NULL)
      incProgress(1/n,detail = 'Finish Parsing Columns step 5')



    })

  }









}})


source('server/observer.R',local = T)








