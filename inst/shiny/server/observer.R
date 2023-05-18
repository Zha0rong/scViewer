observe(if (!is.null(reactivevalue$object_location)&(!reactivevalue$Loaded)){
  disable(id='submit')
  disable(id='Seurat_Object')
  reactivevalue$Loaded=T

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
    output$MainFigure=renderPlot(DimPlot(reactivevalue$SeuratObject,raster = F))
    updateSelectizeInput(session = session,inputId = 'DRreduction',choices =names((reactivevalue$SeuratObject@reductions)),selected = NULL)
    updateSelectizeInput(session = session,inputId = 'GIreduction',choices =names((reactivevalue$SeuratObject@reductions)),selected = NULL)

    updateSelectizeInput(session = session,inputId = 'GenesToInterrogate',choices =rownames(reactivevalue$SeuratObject),selected = NULL,server=T)

    incProgress(1/n,detail = 'Finish parsing Genes')

    updateSelectizeInput(session = session,inputId = 'variabletogroup',choices=colnames(reactivevalue$Experiment_Metadata)
                         [
                           !grepl("nCount",colnames(reactivevalue$Experiment_Metadata))&!grepl("nFeature",colnames(reactivevalue$Experiment_Metadata))&!grepl("^percent.",colnames(reactivevalue$Experiment_Metadata))
                         ]
                         ,selected = 'ShinyGroup')

    incProgress(1/n,detail = 'Finish Parsing Columns step 1')

    updateSelectizeInput(session = session,inputId = 'variabletosplit',choices=colnames(reactivevalue$Experiment_Metadata)
                         [
                           !grepl("nCount",colnames(reactivevalue$Experiment_Metadata))&!grepl("nFeature",colnames(reactivevalue$Experiment_Metadata))&!grepl("^percent.",colnames(reactivevalue$Experiment_Metadata))
                         ]
                         ,selected = 'ShinyGroup')
    incProgress(1/n,detail = 'Finish Parsing Columns step 2')

    updateSelectizeInput(session = session,inputId = 'SampleColumn',choices=colnames(reactivevalue$Experiment_Metadata)
                         [
                           !grepl("nCount",colnames(reactivevalue$Experiment_Metadata))&!grepl("nFeature",colnames(reactivevalue$Experiment_Metadata))&!grepl("^percent.",colnames(reactivevalue$Experiment_Metadata))
                         ]
                         ,selected = NULL)
    incProgress(1/n,detail = 'Finish Parsing Columns step 3')

    updateSelectizeInput(session = session,inputId = 'PlotGroup',choices=colnames(reactivevalue$Experiment_Metadata)
                         [
                           !grepl("nCount",colnames(reactivevalue$Experiment_Metadata))&!grepl("nFeature",colnames(reactivevalue$Experiment_Metadata))&!grepl("^percent.",colnames(reactivevalue$Experiment_Metadata))
                         ])

    incProgress(1/n,detail = 'Finish Parsing Columns step 4')

    updateSelectizeInput(session = session,inputId = 'BarGraph1',choices=colnames(reactivevalue$Experiment_Metadata)
                         [
                           !grepl("nCount",colnames(reactivevalue$Experiment_Metadata))&!grepl("nFeature",colnames(reactivevalue$Experiment_Metadata))&!grepl("^percent.",colnames(reactivevalue$Experiment_Metadata))
                         ]
                         ,selected = NULL)

    updateSelectizeInput(session = session,inputId = 'BarGraph2',choices=colnames(reactivevalue$Experiment_Metadata)
                         [
                           !grepl("nCount",colnames(reactivevalue$Experiment_Metadata))&!grepl("nFeature",colnames(reactivevalue$Experiment_Metadata))&!grepl("^percent.",colnames(reactivevalue$Experiment_Metadata))
                         ]
                         ,selected = NULL)

    DGE_Group_Candidate=c()
    for (i in colnames(reactivevalue$SeuratObject@meta.data)[
      !grepl("nCount",colnames(reactivevalue$SeuratObject@meta.data))&!grepl("nFeature",colnames(reactivevalue$SeuratObject@meta.data))&!grepl("^percent.",colnames(reactivevalue$SeuratObject@meta.data))
      &!grepl("_res.",colnames(reactivevalue$SeuratObject@meta.data))
    ]) {
      DGE_Group_Candidate=c(DGE_Group_Candidate,paste0(i,'-',unique(as.character(reactivevalue$Experiment_Metadata[,i]))))
    }
    updateSelectizeInput(session = session,inputId = 'Assay',choices=names(reactivevalue$SeuratObject@assays),selected = NULL)
    updateSelectizeInput(session = session,inputId = 'DGEGroup1',choices=DGE_Group_Candidate,selected = NULL)
    updateSelectizeInput(session = session,inputId = 'DGEGroup2',choices=DGE_Group_Candidate,selected = NULL)
    incProgress(1/n,detail = 'Finish Parsing Columns step 5')



  })

})

observeEvent( input$Seurat_Object, {

  if (is.null(input$Seurat_Object)) return()
  reactivevalue$object_location=input$Seurat_Object$datapath
  output$object_location=renderText(reactivevalue$object_location)

  })


observeEvent( input$submit, {
  if (!is.null(reactivevalue$object_location)&(!reactivevalue$Loaded)){
    reactivevalue$Loaded=T
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
    output$MainFigure=renderPlot(DimPlot(reactivevalue$SeuratObject,raster = F))
    updateSelectizeInput(session = session,inputId = 'reduction',choices =names((reactivevalue$SeuratObject@reductions)),selected = NULL)
    updateSelectizeInput(session = session,inputId = 'GIreduction',choices =names((reactivevalue$SeuratObject@reductions)),selected = NULL)

    updateSelectizeInput(session = session,inputId = 'GenesToInterrogate',choices =rownames(reactivevalue$SeuratObject),selected = NULL,server=T)

    incProgress(1/n,detail = 'Finish parsing Genes')

    updateSelectizeInput(session = session,inputId = 'variabletogroup',choices=colnames(reactivevalue$Experiment_Metadata)
                         [
      !grepl("nCount",colnames(reactivevalue$Experiment_Metadata))&!grepl("nFeature",colnames(reactivevalue$Experiment_Metadata))&!grepl("^percent.",colnames(reactivevalue$Experiment_Metadata))
    ]
    ,selected = 'ShinyGroup')

    incProgress(1/n,detail = 'Finish Parsing Columns step 1')

    updateSelectizeInput(session = session,inputId = 'variabletosplit',choices=colnames(reactivevalue$Experiment_Metadata)
                         [
      !grepl("nCount",colnames(reactivevalue$Experiment_Metadata))&!grepl("nFeature",colnames(reactivevalue$Experiment_Metadata))&!grepl("^percent.",colnames(reactivevalue$Experiment_Metadata))
    ]
    ,selected = 'ShinyGroup')
    incProgress(1/n,detail = 'Finish Parsing Columns step 2')

    updateSelectizeInput(session = session,inputId = 'SampleColumn',choices=colnames(reactivevalue$Experiment_Metadata)
                         [
      !grepl("nCount",colnames(reactivevalue$Experiment_Metadata))&!grepl("nFeature",colnames(reactivevalue$Experiment_Metadata))&!grepl("^percent.",colnames(reactivevalue$Experiment_Metadata))
    ]
    ,selected = NULL)
    incProgress(1/n,detail = 'Finish Parsing Columns step 3')

    updateSelectizeInput(session = session,inputId = 'PlotGroup',choices=colnames(reactivevalue$Experiment_Metadata)
                         [
      !grepl("nCount",colnames(reactivevalue$Experiment_Metadata))&!grepl("nFeature",colnames(reactivevalue$Experiment_Metadata))&!grepl("^percent.",colnames(reactivevalue$Experiment_Metadata))
    ])

    incProgress(1/n,detail = 'Finish Parsing Columns step 4')

    updateSelectizeInput(session = session,inputId = 'BarGraph1',choices=colnames(reactivevalue$Experiment_Metadata)
                         [
                           !grepl("nCount",colnames(reactivevalue$Experiment_Metadata))&!grepl("nFeature",colnames(reactivevalue$Experiment_Metadata))&!grepl("^percent.",colnames(reactivevalue$Experiment_Metadata))
                         ]
                         ,selected = NULL)

    updateSelectizeInput(session = session,inputId = 'BarGraph2',choices=colnames(reactivevalue$Experiment_Metadata)
                         [
                           !grepl("nCount",colnames(reactivevalue$Experiment_Metadata))&!grepl("nFeature",colnames(reactivevalue$Experiment_Metadata))&!grepl("^percent.",colnames(reactivevalue$Experiment_Metadata))
                         ]
                         ,selected = NULL)

    DGE_Group_Candidate=c()
    for (i in colnames(reactivevalue$SeuratObject@meta.data)[
      !grepl("nCount",colnames(reactivevalue$SeuratObject@meta.data))&!grepl("nFeature",colnames(reactivevalue$SeuratObject@meta.data))&!grepl("^percent.",colnames(reactivevalue$SeuratObject@meta.data))
     &!grepl("_res.",colnames(reactivevalue$SeuratObject@meta.data))
    ]) {
      DGE_Group_Candidate=c(DGE_Group_Candidate,paste0(i,'-',unique(as.character(reactivevalue$Experiment_Metadata[,i]))))
    }
    updateSelectizeInput(session = session,inputId = 'Assay',choices=names(reactivevalue$SeuratObject@assays),selected = NULL)
    updateSelectizeInput(session = session,inputId = 'DGEGroup1',choices=DGE_Group_Candidate,selected = NULL)
    updateSelectizeInput(session = session,inputId = 'DGEGroup2',choices=DGE_Group_Candidate,selected = NULL)
    incProgress(1/n,detail = 'Finish Parsing Columns step 5')



    })

}
}
)

observeEvent( input$SampleColumn, {
  if ((input$SampleColumn%in%colnames(reactivevalue$Experiment_Metadata))&!is.null(reactivevalue$SeuratObject)){
    updateSelectInput(session = session,inputId = 'SampletoSubset',choices=unique(reactivevalue$Experiment_Metadata[,input$SampleColumn]),selected = NULL)
  }
})




BarGraphListener <- reactive({
  list(input$BarGraph1,input$BarGraph2)
})
observeEvent(BarGraphListener()
             ,{
               if ((input$BarGraph1%in%colnames(reactivevalue$Experiment_Metadata))&(input$BarGraph2%in%colnames(reactivevalue$Experiment_Metadata))&!is.null(reactivevalue$SeuratObject)) {
                 if (input$BarGraph1!=input$BarGraph2){
                   temp=data.frame(table(reactivevalue$Experiment_Metadata[,input$BarGraph1],reactivevalue$Experiment_Metadata[,input$BarGraph2]))
                   colnames(temp)=c("Variable1","Variable2",'CellNumber')#

                   output$BarPlot=renderPlot(ggplot(temp,aes(
                     x=CellNumber,y=Variable1,fill=Variable2
                   ))+geom_bar(stat = 'identity',position = 'fill')+xlab('Cell Number')+ylab(input$BarGraph1)+labs(fill=input$BarGraph2) )
                 }
               }

             }
)

ReductionListener <- reactive({
  list(input$DRreduction,input$variabletogroup,input$variabletosplit,input$SampletoSubset,input$DRLabel)
})
observeEvent(ReductionListener(),
             {
               if (is.null(input$DRreduction)) return()
               if (input$variabletosplit%in%colnames(reactivevalue$Experiment_Metadata)){
                 if (!is.null(input$SampletoSubset)) {
                   kept=rownames(reactivevalue$Experiment_Metadata)[!reactivevalue$Experiment_Metadata[,input$SampleColumn]%in%input$SampletoSubset]
                   if (length(kept)!=0) {
                     temp=subset(reactivevalue$SeuratObject,cells=kept)
                     output$tsne=renderPlot(DimPlot(temp,reduction = input$DRreduction,
                                                    group.by = input$variabletogroup,split.by = input$variabletosplit,
                                                    ncol=ifelse(length(unique(reactivevalue$Experiment_Metadata[,input$variabletosplit]))==1,yes = 1,no=2),raster = F,label = input$DRLabel))

                   } else {
                     output$tsne=renderPlot(DimPlot(reactivevalue$SeuratObject,reduction = input$DRreduction,
                                                    group.by = input$variabletogroup,split.by = input$variabletosplit,
                                                    ncol=ifelse(length(unique(reactivevalue$Experiment_Metadata[,input$variabletosplit]))==1,yes = 1,no=2),raster = F,label = input$DRLabel))

                   }

                 } else {
                   output$tsne=renderPlot(DimPlot(reactivevalue$SeuratObject,reduction = input$DRreduction,
                                                  group.by = input$variabletogroup,split.by = input$variabletosplit,
                                                  ncol=ifelse(length(unique(reactivevalue$Experiment_Metadata[,input$variabletosplit]))==1,yes = 1,no=2),
                                                  raster = F,label = input$DRLabel))

                 }
               }

             })

GenesToInterrogateListener <- reactive({
  list(input$GenesToInterrogate,input$PlotGroup,input$GIreduction,input$GILabel)
})

observeEvent( GenesToInterrogateListener(), {
  if (is.null(input$GenesToInterrogate)) return()
  reactivevalue$temp=reactivevalue$SeuratObject
  Idents(reactivevalue$temp)=input$PlotGroup
  output$FeaturePlot=renderPlot(FeaturePlot(reactivevalue$temp,features = input$GenesToInterrogate,order = T,label = input$GILabel,reduction = input$GIreduction))

  output$ViolinPlot=renderPlot(VlnPlot(reactivevalue$temp,assay = 'RNA',features = input$GenesToInterrogate,group.by = input$PlotGroup,
                                       pt.size = ifelse(ncol(reactivevalue$SeuratObject)>1000,yes=0,no=NULL)))
  output$RidgePlot=renderPlot(RidgePlot(reactivevalue$temp,assay = 'RNA',features = input$GenesToInterrogate,group.by = input$PlotGroup))

  if (length(input$GenesToInterrogate)==1) {
    globalstats=list()
    for (i in 1:length(unique(reactivevalue$SeuratObject@meta.data[[input$PlotGroup]]))) {
      temp=summary(reactivevalue$SeuratObject@assays$RNA@data[input$GenesToInterrogate,colnames(reactivevalue$SeuratObject@assays$RNA@data)%in%
                                                                rownames(reactivevalue$SeuratObject@meta.data)[
                                                                  reactivevalue$SeuratObject@meta.data[[input$PlotGroup]]==
                                                                    unique(reactivevalue$SeuratObject@meta.data[[input$PlotGroup]])[i]
                                                                ]])
      globalstat=data.frame(t(as.matrix(temp)))
      colnames(globalstat)=c('Minimum',
                              'Quantile 25th',
                              'Median',
                              'Mean',
                              'Quantile 75th',
                              'Maximum'
      )
      rownames(globalstat)=input$GenesToInterrogate
      globalstat$gene=input$GenesToInterrogate
      globalstat$Group=unique(reactivevalue$SeuratObject@meta.data[[input$PlotGroup]])[i]
      globalstats[[i]]=data.frame(globalstat)
    }
    globalstats=do.call(rbind,globalstats)
    globalstats=data.frame(globalstats)
    rownames(globalstats)=seq(1,nrow(globalstats))
    globalstats=globalstats[,c(ncol(globalstats),ncol(globalstats)-1,seq(1,ncol(globalstats)-2))]


  } else {
    globalstats=list()
    for (i in 1:length(input$GenesToInterrogate)) {
      globalstat=list()

      gene=input$GenesToInterrogate[i]

      for (j in 1:length(unique(reactivevalue$SeuratObject@meta.data[[input$PlotGroup]]))) {
        temp=summary(reactivevalue$SeuratObject@assays$RNA@data[gene,colnames(reactivevalue$SeuratObject@assays$RNA@data)%in%
                                                                  rownames(reactivevalue$SeuratObject@meta.data)[
                                                                    reactivevalue$SeuratObject@meta.data[[input$PlotGroup]]==
                                                                      unique(reactivevalue$SeuratObject@meta.data[[input$PlotGroup]])[j]
                                                                  ]])

        stat=data.frame(t(as.matrix(temp)))
        colnames(stat)=c('Minimum',
                               'Quantile 25th',
                               'Median',
                               'Mean',
                               'Quantile 75th',
                               'Maximum'
        )
        rownames(stat)=gene
        stat$gene=gene
        stat$Group=unique(reactivevalue$SeuratObject@meta.data[[input$PlotGroup]])[j]
        globalstat[[j]]=data.frame(stat)

      }
      globalstat=do.call(rbind,globalstat)
      globalstat=data.frame(globalstat)
      globalstat=globalstat[,c(ncol(globalstat),ncol(globalstat)-1,seq(1,ncol(globalstat)-2))]

      globalstats[[i]]=globalstat
    }
    globalstats=do.call(rbind,globalstats)
    globalstats=data.frame(globalstats)
    rownames(globalstats)=seq(1,nrow(globalstats))
    globalstats$gene=input$GenesToInterrogate
  }

  output$GlobalStats=DT::renderDataTable(DT::datatable(globalstats,editable = F, options = list(dom = 'Bfrtip'), filter = list(position = "top")),server = T)
  reactivevalue$GeneStats=globalstats
  output$Gene.Expression.Statistics.downloadData <- downloadHandler(
    filename = function() {
      # Use the selected dataset as the suggested file name
      paste0("Gene.Expression.Statistics", ".tsv")
    },
    content = function(file) {
      # Write the dataset to the `file` that will be downloaded
      write.table(reactivevalue$GeneStats, file,sep = '\t',quote=F)
    }
  )
})




DGEListener <- reactive({
  list(input$DGEGroup1,input$DGEGroup2)
})


observeEvent(DGEListener(),
             {
               if(!is.null(input$DGEGroup1)&&!is.null(input$DGEGroup2)){
                 CellRanch=reactivevalue$Experiment_Metadata[,
                                                             colnames(reactivevalue$Experiment_Metadata)[
                                                               !grepl("nCount",colnames(reactivevalue$Experiment_Metadata))&!grepl("nFeature",colnames(reactivevalue$Experiment_Metadata))&!grepl("^percent.",colnames(reactivevalue$Experiment_Metadata))
                                                               &!grepl("_res.",colnames(reactivevalue$Experiment_Metadata))
                                                             ]]
                 Group1Wrangler=list()
                 Group2Wrangler=list()
                 for (i in 1:length(input$DGEGroup1)) {
                   candidates=unlist(strsplit(input$DGEGroup1[i],'-'))[1]
                   tryCatch({
                     cells=rownames(CellRanch)[paste0((candidates),'-',as.character(CellRanch[,candidates]))==input$DGEGroup1[i]]
                     if (candidates%in%names(Group1Wrangler)) {
                       Group1Wrangler[[candidates]]=unique(c(Group1Wrangler[[candidates]],cells))

                     } else {
                       Group1Wrangler[[candidates]]=unique(cells)
                     }
                   },
                   error=function(err){
                     print('Error')
                   })


                 }

                 for (i in 1:length(input$DGEGroup2)) {
                   candidates=unlist(strsplit(input$DGEGroup2[i],'-'))[1]
                   tryCatch({
                     cells=rownames(CellRanch)[paste0((candidates),'-',as.character(CellRanch[,candidates]))==input$DGEGroup2[i]]
                     if (candidates%in%names(Group2Wrangler)) {
                       Group2Wrangler[[candidates]]=unique(c(Group2Wrangler[[candidates]],cells))

                     } else {
                       Group2Wrangler[[candidates]]=unique(cells)
                     }
                   },
                   error=function(err){
                     print('Error')
                   })
                 }


                 Group1Wrangled=Reduce(intersect,Group1Wrangler)
                 Group2Wrangled=Reduce(intersect,Group2Wrangler)

                 analysis=data.frame(Group=c('Group 1','Group 2'),
                                     CellNumber=c(
                                       length(Group1Wrangled[!Group1Wrangled%in%intersect(Group1Wrangled,Group2Wrangled)])-length(
                                         intersect(Group1Wrangled[!Group1Wrangled%in%intersect(Group1Wrangled,Group2Wrangled)],
                                                   Group2Wrangled[!Group2Wrangled%in%intersect(Group1Wrangled,Group2Wrangled)]
                                         )
                                       ),
                                       length(Group2Wrangled[!Group2Wrangled%in%intersect(Group1Wrangled,Group2Wrangled)])-length(
                                         intersect(Group1Wrangled[!Group1Wrangled%in%intersect(Group1Wrangled,Group2Wrangled)],
                                                   Group2Wrangled[!Group2Wrangled%in%intersect(Group1Wrangled,Group2Wrangled)]
                                         )
                                       )
                                     ))
                 reactivevalue$analysis=analysis
                 reactivevalue$Group1Wrangled=Group1Wrangled[!Group1Wrangled%in%intersect(Group1Wrangled[!Group1Wrangled%in%intersect(Group1Wrangled,Group2Wrangled)],
                                                                                          Group2Wrangled[!Group2Wrangled%in%intersect(Group1Wrangled,Group2Wrangled)]
                 )]
                 reactivevalue$Group2Wrangled=Group2Wrangled[!Group2Wrangled%in%intersect(Group1Wrangled[!Group1Wrangled%in%intersect(Group1Wrangled,Group2Wrangled)],
                                                                                          Group2Wrangled[!Group2Wrangled%in%intersect(Group1Wrangled,Group2Wrangled)]
                 )]

                 output$GroupNumber=DT::renderDataTable(DT::datatable((analysis), options = list(dom = 'Bfrtip'), filter = list(position = "top")),server = T
                  )
               }
             })


observeEvent(input$submitDGE, {
  if (all(reactivevalue$analysis[,'CellNumber']!=0)) {
    FindMarkerstemp=reactivevalue$SeuratObject
    DefaultAssay(FindMarkerstemp)=input$Assay
    withProgress(message = 'Start to do differential expression analysis',value = 0, {

      n=3
      FindMarkerstemp=subset(FindMarkerstemp,cells=unique(c(reactivevalue$Group1Wrangled,
                                                            reactivevalue$Group2Wrangled)))
      incProgress(1/n,detail = 'Start Data Normalization')

      FindMarkerstemp=NormalizeData(FindMarkerstemp)
      incProgress(1/n,detail = 'Start Differential Expression Analysis')

      FindMarkerstemp@meta.data$Group=ifelse(rownames(FindMarkerstemp@meta.data)%in%reactivevalue$Group1Wrangled,
                                             yes='Group1',no='Group2')
      Results=FindMarkers(FindMarkerstemp,ident.1='Group1',ident.2='Group2',group.by='Group',
                          assay=input$Assay)


      incProgress(1/n,detail = 'Finish Differential Expression Analysis')
    })
    Results$gene=rownames(Results)
    output$DifferentialExpressionAnalysisResults=DT::renderDataTable(DT::datatable(Results, options = list(dom = 'Bfrtip'), filter = list(position = "top")),server = T)


    reactivevalue$DifferentialExpressionAnalysisResults=Results
    output$downloadData <- downloadHandler(
      filename = function() {
        # Use the selected dataset as the suggested file name
        paste0("Differential.Expression.Analysis", ".tsv")
      },
      content = function(file) {
        # Write the dataset to the `file` that will be downloaded
        write.table(reactivevalue$DifferentialExpressionAnalysisResults, file,sep = '\t',quote=F)
      }
    )
  }

})











