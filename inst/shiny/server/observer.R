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
    reactivevalue$metadatacolumn=colnames(reactivevalue$Experiment_Metadata)[!grepl("nCount",colnames(reactivevalue$Experiment_Metadata))&!grepl("nFeature",colnames(reactivevalue$Experiment_Metadata))&!grepl("^percent.",colnames(reactivevalue$Experiment_Metadata))]
    
    reactivevalue$assays=names(reactivevalue$SeuratObject@assays)
    reactivevalue$reduction=names((reactivevalue$SeuratObject@reductions))
    reactivevalue$reduction=reactivevalue$reduction[!grepl('pca',reactivevalue$reduction,ignore.case = T)]
    reactivevalue$reduction=reactivevalue$reduction[!grepl('harmony',reactivevalue$reduction,ignore.case = T)]
    reactivevalue$genes=rownames(reactivevalue$SeuratObject)
    
    DefaultAssay(reactivevalue$SeuratObject)='RNA'
    output$MainFigure=renderPlot(DimPlot(reactivevalue$SeuratObject,raster = F))
    updateSelectizeInput(session = session,inputId = 'Reference_Column',choices =reactivevalue$metadatacolumn,selected = NULL,server=T)
    updateSelectizeInput(session = session,inputId = 'FindMarkersVariable',choices =reactivevalue$metadatacolumn,selected = NULL,server=T)
    updateSelectizeInput(session = session,inputId = 'Annotate_Group',choices =reactivevalue$metadatacolumn,server=T)
    
    
    updateSelectizeInput(session = session,inputId = 'DRreduction',choices =reactivevalue$reduction,selected = NULL,server=T)
    updateSelectizeInput(session = session,inputId = 'GIreduction',choices =reactivevalue$reduction,selected = NULL,server=T)

    updateSelectizeInput(session = session,inputId = 'GenesToInterrogate',choices =reactivevalue$genes,selected = NULL,server=T)

    incProgress(1/n,detail = 'Finish parsing Genes')

    updateSelectizeInput(session = session,inputId = 'variabletogroup',choices=reactivevalue$metadatacolumn
                         ,selected = 'ShinyGroup')

    incProgress(1/n,detail = 'Finish Parsing Columns step 1')

    updateSelectizeInput(session = session,inputId = 'variabletosplit',choices=reactivevalue$metadatacolumn
                         ,selected = 'ShinyGroup')
    incProgress(1/n,detail = 'Finish Parsing Columns step 2')

    updateSelectizeInput(session = session,inputId = 'SampleColumn',choices=reactivevalue$metadatacolumn
                         ,selected = NULL)
    incProgress(1/n,detail = 'Finish Parsing Columns step 3')

    updateSelectizeInput(session = session,inputId = 'PlotGroup',choices=reactivevalue$metadatacolumn)

    incProgress(1/n,detail = 'Finish Parsing Columns step 4')

    updateSelectizeInput(session = session,inputId = 'BarGraph1',choices=reactivevalue$metadatacolumn,selected = NULL)

    updateSelectizeInput(session = session,inputId = 'BarGraph2',choices=reactivevalue$metadatacolumn
                         ,selected = NULL)

    DGE_Group_Candidate=c()
    for (i in reactivevalue$metadatacolumn) {
      DGE_Group_Candidate=c(DGE_Group_Candidate,paste0(i,'-',unique(as.character(reactivevalue$Experiment_Metadata[,i]))))
    }
    updateSelectizeInput(session = session,inputId = 'Assay',choices=reactivevalue$assays,selected = NULL)
    updateSelectizeInput(session = session,inputId = 'GenesToInterrogateAssay',choices=reactivevalue$assays,selected = NULL)
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
  list(input$BarGraph1,input$BarGraph2,input$BarGraphPercentage)
})
observeEvent(BarGraphListener(),{
               if ((input$BarGraph1%in%colnames(reactivevalue$Experiment_Metadata))&(input$BarGraph2%in%colnames(reactivevalue$Experiment_Metadata))&!is.null(reactivevalue$SeuratObject)) {
                 if (input$BarGraph1!=input$BarGraph2){
                   temp=data.frame(table(reactivevalue$Experiment_Metadata[,input$BarGraph1],reactivevalue$Experiment_Metadata[,input$BarGraph2]))
                   colnames(temp)=c("Variable1","Variable2",'CellNumber')
                   if (input$BarGraphPercentage){
                   output$BarPlot=renderPlot(ggplot(temp,aes(
                     x=CellNumber,y=Variable1,fill=Variable2
                   ))+geom_bar(stat = 'identity',position = 'fill')+xlab('Cell Number')+ylab(input$BarGraph1)+labs(fill=input$BarGraph2) )} else {
                     output$BarPlot=renderPlot(ggplot(temp,aes(
                       x=CellNumber,y=Variable1,fill=Variable2
                     ))+geom_bar(stat = 'identity',position = 'stack')+xlab('Cell Number')+ylab(input$BarGraph1)+labs(fill=input$BarGraph2) )
                   }
                   output$BarPlotStats=DT::renderDataTable(DT::datatable(temp,editable = F,rownames= FALSE, filter = list(position = "top")),server = T)

                 }
               }

             })

ReductionListener <- reactive({
  list(input$DRreduction,input$variabletogroup,input$variabletosplit,input$SampletoSubset,input$DRLabel)
})
observeEvent(ReductionListener(),{
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
  list(input$GenesToInterrogate,input$PlotGroup,input$GIreduction,input$GILabel,input$GenesToInterrogateAssay,input$FeaturePlotOverlay)
})

observeEvent( GenesToInterrogateListener(), {
  if (is.null(input$GenesToInterrogate)|is.null(input$PlotGroup)) return()
  if ((input$PlotGroup!=''&input$GenesToInterrogateAssay!='')) {
  reactivevalue$temp=NULL
  reactivevalue$temp=reactivevalue$SeuratObject
  DefaultAssay(reactivevalue$temp)=input$GenesToInterrogateAssay
  Idents(reactivevalue$temp)=input$PlotGroup
  output$FeaturePlot=NULL
  if (input$FeaturePlotOverlay) {
    if (length(input$GenesToInterrogate)==1) {
      output$FeaturePlot=renderPlot(FeaturePlot(reactivevalue$temp,features = input$GenesToInterrogate,order = T,
                                                label = input$GILabel,reduction = input$GIreduction))
    }
    else if (length(input$GenesToInterrogate)==2) {
      plotinput=reactive({temp=FeaturePlot(reactivevalue$temp,features = input$GenesToInterrogate,order = T,
                                                   label = input$GILabel,reduction = input$GIreduction,blend = T,ncol = 2)
      
      temp$patches$layout$ncol=2
      temp=ggplotify::as.ggplot(patchwork::patchworkGrob(temp))
      temp
      })
      
      output$FeaturePlot=renderPlot({

        plotinput()
      })
    } else {
      reactivevalue$temp=AddModuleScore(reactivevalue$temp,features = list(input$GenesToInterrogate))
      reactivevalue$temp$Overlay_Expression=reactivevalue$temp$Cluster1
      output$FeaturePlot=renderPlot(FeaturePlot(reactivevalue$temp,features = 'Overlay_Expression',order = T,
                                                label = input$GILabel,reduction = input$GIreduction))
    }
  } else {
  output$FeaturePlot=renderPlot(FeaturePlot(reactivevalue$temp,features = input$GenesToInterrogate,order = T,
                                            label = input$GILabel,reduction = input$GIreduction))
  }

  output$ViolinPlot=renderPlot(VlnPlot(reactivevalue$temp,assay = input$GenesToInterrogateAssay,features = input$GenesToInterrogate,group.by = input$PlotGroup,
                                       pt.size = ifelse(ncol(reactivevalue$SeuratObject)>1000,yes=0,no=NULL)))
  output$RidgePlot=renderPlot(RidgePlot(reactivevalue$temp,assay = input$GenesToInterrogateAssay,features = input$GenesToInterrogate,group.by = input$PlotGroup))

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
    variables=unique(as.character(reactivevalue$SeuratObject@meta.data[[input$PlotGroup]]))
    genes=input$GenesToInterrogate
    count=1
    for (i in 1:length(variables)) {

      for (j in 1:length(genes)) {
        gene=genes[j]
        variable=variables[i]
        temp=summary(reactivevalue$SeuratObject@assays$RNA@data[gene,colnames(reactivevalue$SeuratObject@assays$RNA@data)%in%
                                                                  rownames(reactivevalue$SeuratObject@meta.data)[
                                                                    as.character(reactivevalue$SeuratObject@meta.data[[input$PlotGroup]])==variable
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
        stat$Group=variable
        stat=stat[,c(ncol(stat),ncol(stat)-1,seq(1,ncol(stat)-2))]
        
        globalstats[[count]]=data.frame(stat)
        count=count+1
      }
  }
    globalstats=do.call(rbind,globalstats)
    globalstats=data.frame(globalstats)
    rownames(globalstats)=seq(1,nrow(globalstats))
  }

  output$GlobalStats=DT::renderDataTable(DT::datatable(globalstats,editable = F, options = list(dom = 'Bfrtip'),rownames= FALSE, filter = list(position = "top")),server = T)
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
  genes=input$GenesToInterrogate
  if (length(genes)>1){
  temp=reactivevalue$SeuratObject@assays$RNA@data[genes,]
  fun=function(x) {return(all(x>0))}
  temp=colnames(temp)[apply(temp,2,fun)]
  } else {
    temp=reactivevalue$SeuratObject@assays$RNA@data[genes,]
    temp=colnames(reactivevalue$SeuratObject@assays$RNA@data)[temp>0]
  }
  if (length(temp)>0) {
    temp=subset(reactivevalue$SeuratObject,cells=temp)
    temp=data.frame(table(temp[[input$PlotGroup]]))
    output$CellStatsTableBasedonGeneFiltering=DT::renderDataTable(DT::datatable(temp,editable = F, options = list(dom = 'Bfrtip'),rownames= FALSE, filter = list(position = "top")),server = T)

  }
  
}
})




DGEListener <- reactive({
  list(input$DGEGroup1,input$DGEGroup2)
})


observeEvent(DGEListener(),{
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

                 output$GroupNumber=DT::renderDataTable(DT::datatable((analysis), options = list(dom = 'Bfrtip'),rownames= FALSE, filter = list(position = "top")),server = T
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
    data=Results[,c('gene','avg_log2FC','p_val','p_val_adj')]
    output$DifferentialExpressionAnalysisResults=DT::renderDataTable(DT::datatable(data, options = list(dom = 'Bfrtip'),rownames= FALSE, filter = list(position = "top")),server = T)


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

AnnotationListener=reactive({list(input$Reference_Column,input$AnnotationName)})
observeEvent(AnnotationListener(),{
  if (is.null(input$Reference_Column)|is.null(reactivevalue$Experiment_Metadata)) return()
  if (input$Reference_Column%in%colnames(reactivevalue$Experiment_Metadata)){
  reactivevalue$annotationdata=data.frame(referencecolumn=unique(as.character(reactivevalue$Experiment_Metadata[[input$Reference_Column]])),
                                          row.names = unique(as.character(reactivevalue$Experiment_Metadata[[input$Reference_Column]])))
  reactivevalue$annotationdata$NewAnnotation=''
  reactivevalue$annotationdata_record=reactivevalue$annotationdata
  output$AnnotationTable=DT::renderDataTable(DT::datatable(reactivevalue$annotationdata,editable = list(target = "cell", disable = list(columns = c(0))),rownames= FALSE))
}
  })


observeEvent((input$AnnotationTable_cell_edit), {
  reactivevalue$annotationdata_record[input$AnnotationTable_cell_edit$row,'NewAnnotation']=input$AnnotationTable_cell_edit$value
})

observeEvent(input$AddAnnotation,{
  previousIdents=Idents(reactivevalue$SeuratObject)
  NewIdents=reactivevalue$annotationdata_record$NewAnnotation

  names(NewIdents)=reactivevalue$annotationdata_record$referencecolumn

  Idents(reactivevalue$SeuratObject)=input$Reference_Column
  reactivevalue$SeuratObject=RenameIdents(reactivevalue$SeuratObject,NewIdents)
  if (is.null(input$AnnotationName)) {
    reactivevalue$SeuratObject@meta.data[[paste0('New Annotation Column',ncol(reactivevalue$SeuratObject@meta.data))]]=Idents(reactivevalue$SeuratObject)
    }else {
    reactivevalue$SeuratObject@meta.data[[input$AnnotationName]]=Idents(reactivevalue$SeuratObject)
      
    }
  reactivevalue$Experiment_Metadata=reactivevalue$SeuratObject@meta.data
  reactivevalue$metadatacolumn=colnames(reactivevalue$Experiment_Metadata)[!grepl("nCount",colnames(reactivevalue$Experiment_Metadata))&!grepl("nFeature",colnames(reactivevalue$Experiment_Metadata))&!grepl("^percent.",colnames(reactivevalue$Experiment_Metadata))]
  
  DGE_Group_Candidate=c()
  for (i in reactivevalue$metadatacolumn) {
    DGE_Group_Candidate=c(DGE_Group_Candidate,paste0(i,'-',unique(as.character(reactivevalue$Experiment_Metadata[,i]))))
  }
  
  
  updateSelectizeInput(session = session,inputId = 'variabletogroup',choices=reactivevalue$metadatacolumn
                       ,selected = 'ShinyGroup')
  

  updateSelectizeInput(session = session,inputId = 'variabletosplit',choices=reactivevalue$metadatacolumn
                       ,selected = 'ShinyGroup')

  updateSelectizeInput(session = session,inputId = 'SampleColumn',choices=reactivevalue$metadatacolumn
                       ,selected = NULL)

  updateSelectizeInput(session = session,inputId = 'PlotGroup',choices=reactivevalue$metadatacolumn)
  

  updateSelectizeInput(session = session,inputId = 'BarGraph1',choices=reactivevalue$metadatacolumn,selected = NULL)
  
  updateSelectizeInput(session = session,inputId = 'BarGraph2',choices=reactivevalue$metadatacolumn
                       ,selected = NULL)
  updateSelectizeInput(session = session,inputId = 'Reference_Column',choices =reactivevalue$metadatacolumn,selected = NULL,server=T)
  updateSelectizeInput(session = session,inputId = 'Assay',choices=reactivevalue$assays,selected = NULL)
  updateSelectizeInput(session = session,inputId = 'DGEGroup1',choices=DGE_Group_Candidate,selected = NULL)
  updateSelectizeInput(session = session,inputId = 'DGEGroup2',choices=DGE_Group_Candidate,selected = NULL)
  updateSelectizeInput(session = session,inputId = 'FindMarkersVariable',choices=reactivevalue$metadatacolumn
                       ,selected = NULL)
  updateSelectizeInput(session = session,inputId = 'Annotate_Group',choices =reactivevalue$metadatacolumn,server=T)
  
})

observeEvent(input$submitFindMarkers, {
  variable=input$FindMarkersVariable
  groups=unique(as.character(reactivevalue$SeuratObject@meta.data[[variable]]))
  cells=c()
  sampling.size=500
  withProgress(message = 'Doing Marker Finder',value = 0, {
  n=3
  incProgress(1/n,detail = 'Start Sampling')
  
  for (i in groups) {
    temp=rownames(reactivevalue$SeuratObject@meta.data)[reactivevalue$SeuratObject@meta.data[[variable]]==i]
    if (length(temp)<=sampling.size) {
      cells=c(cells,temp)
    } else {
      cells=c(cells,sample(temp,sampling.size,replace = F))
      
    }
  }
  incProgress(1/n,detail = 'Start Normalizing')
  
  temp=subset(reactivevalue$SeuratObject,cells=cells)
  DefaultAssay(temp)='RNA'
  temp=NormalizeData(temp)
  Idents(temp)=variable
  incProgress(1/n,detail = 'Start Marker Finding')
  
  temp=FindAllMarkers(temp,only.pos = T,logfc.threshold = 1,test.use = 'wilcox')
  temp=temp[,c('gene','avg_log2FC','p_val','p_val_adj')]
  output$FindMarkersResults=DT::renderDataTable(DT::datatable(temp, options = list(dom = 'Bfrtip'),rownames= FALSE, filter = list(position = "top")),server = T)
  
  })
})





observeEvent(input$Reference, {

  withProgress(message = 'Reference Downloading',value = 0, {
    n=2
    incProgress(1/n,detail = 'Start Downloading Reference.')
    if (input$Reference=='HumanPrimaryCellAtlasData') {
      reactivevalue$reference=celldex::HumanPrimaryCellAtlasData()
    }
    else if (input$Reference=='BlueprintEncodeData') {
      reactivevalue$reference=celldex::BlueprintEncodeData()
      
    }
    else if (input$Reference=='MouseRNAseqData') {
      reactivevalue$reference=celldex::MouseRNAseqData()
      
    }
    else if (input$Reference=='ImmGenData') {
      reactivevalue$reference=celldex::ImmGenData()
      
    }
    else if (input$Reference=='DatabaseImmuneCellExpressionData') {
      reactivevalue$reference=celldex::DatabaseImmuneCellExpressionData()
      
    }
    else if (input$Reference=='NovershternHematopoieticData') {
      reactivevalue$reference=celldex::NovershternHematopoieticData()
      
    }
    else if (input$Reference=='MonacoImmuneData') {
      reactivevalue$reference=celldex::MonacoImmuneData()
      
    }
    incProgress(2/n,detail = 'Finish Downloading Reference.')
    output$ReferenceInformation=DT::renderDataTable(DT::datatable(data.frame(colData(reactivevalue$reference))))
    
  })
})

observeEvent(input$submitAnnotation, {
  withProgress(message = 'Annotation',value = 0, {
  n=3
  reactivevalue$SeuratObject=NormalizeData(reactivevalue$SeuratObject)
  counts=reactivevalue$SeuratObject@assays$RNA@counts
  rownames(counts)=toupper(rownames(counts))
  Normalized=reactivevalue$SeuratObject@assays$RNA@data
  rownames(Normalized)=toupper(rownames(Normalized))
  incProgress(1/n,detail = 'Finished fixing single cell experiment.')
  Fixed.sce=SingleCellExperiment(list(counts=counts,logcounts=Normalized),colData=reactivevalue$SeuratObject@meta.data)
  
  
  
  ref.coldata=reactivevalue$reference@colData
  ref.logcounts=reactivevalue$reference@assays@data@listData$logcounts
  rownames(ref.logcounts)=toupper(rownames(ref.logcounts))
  
  Fixed.ref=SummarizedExperiment(list(logcounts=ref.logcounts),colData=ref.coldata)
  incProgress(2/n,detail = 'Start annotation')
  
  
  if (input$levelofannotation=='Main'){
    reactivevalue$annotationresults=SingleR(Fixed.sce[intersect(rownames(Fixed.ref),rownames(Fixed.sce)),], 
                                          Fixed.ref[intersect(rownames(Fixed.ref),rownames(Fixed.sce)),], 
                                          clusters=reactivevalue$SeuratObject@meta.data[,input$Annotate_Group], labels=Fixed.ref$label.main)
  } else {
    reactivevalue$annotationresults=SingleR(Fixed.sce[intersect(rownames(Fixed.ref),rownames(Fixed.sce)),], 
                                            Fixed.ref[intersect(rownames(Fixed.ref),rownames(Fixed.sce)),], 
                                            clusters=reactivevalue$SeuratObject@meta.data[,input$Annotate_Group], labels=Fixed.ref$label.fine)
    }
  reactivevalue$annotationresults.df=data.frame(reactivevalue$annotationresults)
  
  
  reactivevalue$annotation=data.frame(reactivevalue$annotationresults.df[,seq(ncol(reactivevalue$annotationresults.df)-4,ncol(reactivevalue$annotationresults.df))])
  for (i in 1:ncol(reactivevalue$annotation)) {
    reactivevalue$annotation[,i][is.na(reactivevalue$annotation[,i])]='Unknown'
  }
  
  
  output$AnnotationResults=DT::renderDataTable(DT::datatable(data.frame(reactivevalue$annotation)))
  heatmap.data=reactivevalue$annotationresults.df[,-seq(ncol(reactivevalue$annotationresults.df)-4,ncol(reactivevalue$annotationresults.df))]
  colnames(heatmap.data)=gsub('scores.','',colnames(heatmap.data))
  output$visualized.AnnotationResults=renderPlot(pheatmap(
    heatmap.data,
    cluster_rows = F
  ))
  
  Idents(reactivevalue$SeuratObject)=input$Annotate_Group
  annotated=reactivevalue$annotation$label
  names(annotated)=rownames(reactivevalue$annotation)
  reactivevalue$SeuratObject=RenameIdents(reactivevalue$SeuratObject,annotated)
  reactivevalue$SeuratObject$singleR.annotation=Idents(reactivevalue$SeuratObject)
  
  
  reactivevalue$Experiment_Metadata=reactivevalue$SeuratObject@meta.data
  reactivevalue$metadatacolumn=colnames(reactivevalue$Experiment_Metadata)[!grepl("nCount",colnames(reactivevalue$Experiment_Metadata))&!grepl("nFeature",colnames(reactivevalue$Experiment_Metadata))&!grepl("^percent.",colnames(reactivevalue$Experiment_Metadata))]
  
  DGE_Group_Candidate=c()
  for (i in reactivevalue$metadatacolumn) {
    DGE_Group_Candidate=c(DGE_Group_Candidate,paste0(i,'-',unique(as.character(reactivevalue$Experiment_Metadata[,i]))))
  }
  
  
  updateSelectizeInput(session = session,inputId = 'variabletogroup',choices=reactivevalue$metadatacolumn
                       ,selected = 'ShinyGroup')
  
  
  updateSelectizeInput(session = session,inputId = 'variabletosplit',choices=reactivevalue$metadatacolumn
                       ,selected = 'ShinyGroup')
  
  updateSelectizeInput(session = session,inputId = 'SampleColumn',choices=reactivevalue$metadatacolumn
                       ,selected = NULL)
  
  updateSelectizeInput(session = session,inputId = 'PlotGroup',choices=reactivevalue$metadatacolumn)
  
  
  updateSelectizeInput(session = session,inputId = 'BarGraph1',choices=reactivevalue$metadatacolumn,selected = NULL)
  
  updateSelectizeInput(session = session,inputId = 'BarGraph2',choices=reactivevalue$metadatacolumn
                       ,selected = NULL)
  updateSelectizeInput(session = session,inputId = 'Reference_Column',choices =reactivevalue$metadatacolumn,selected = NULL,server=T)
  updateSelectizeInput(session = session,inputId = 'Assay',choices=reactivevalue$assays,selected = NULL)
  updateSelectizeInput(session = session,inputId = 'DGEGroup1',choices=DGE_Group_Candidate,selected = NULL)
  updateSelectizeInput(session = session,inputId = 'DGEGroup2',choices=DGE_Group_Candidate,selected = NULL)
  updateSelectizeInput(session = session,inputId = 'FindMarkersVariable',choices=reactivevalue$metadatacolumn
                       ,selected = NULL)
  updateSelectizeInput(session = session,inputId = 'Annotate_Group',choices =reactivevalue$metadatacolumn,server=T)
  
  
  
  
  
  
  incProgress(3/n,detail = 'Finished Annotation')
  
  })
  
})






