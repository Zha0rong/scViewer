#### Obtain the Counts matrix and Count table location ####
observeEvent( input$Seurat_Object, {
  if (is.null(input$Seurat_Object)) return()
  #reactivevalue$counts = input$counts
  reactivevalue$object_location=input$Seurat_Object$datapath

  #output$Seurat_Object=readRDS(reactivevalue$object_location)
  print(reactivevalue$object_location)
})

observeEvent( input$submit, {
  if (!is.null(reactivevalue$object_location)){
    print('Start Loading')
    reactivevalue$SeuratObject=readRDS(reactivevalue$object_location)
    reactivevalue$SeuratObject$ShinyGroup='SCViewer'
    DefaultAssay(reactivevalue$SeuratObject)='RNA'
    reactivevalue$Experiment_Metadata=reactivevalue$SeuratObject@meta.data
    updateSelectizeInput(session = session,inputId = 'reduction',choices =names((reactivevalue$SeuratObject@reductions)),selected = NULL)
    updateSelectizeInput(session = session,inputId = 'GenesToInterrogate',choices =rownames(reactivevalue$SeuratObject),selected = NULL)

    updateSelectizeInput(session = session,inputId = 'variabletogroup',choices=colnames(reactivevalue$SeuratObject@meta.data)[
      !grepl("nCount",colnames(reactivevalue$SeuratObject@meta.data))&!grepl("nFeature",colnames(reactivevalue$SeuratObject@meta.data))&!grepl("^percent.",colnames(reactivevalue$SeuratObject@meta.data))
      &!grepl("_res.",colnames(reactivevalue$SeuratObject@meta.data))
    ]
                         ,selected = NULL)
    updateSelectizeInput(session = session,inputId = 'variabletosplit',choices=colnames(reactivevalue$SeuratObject@meta.data)[
      !grepl("nCount",colnames(reactivevalue$SeuratObject@meta.data))&!grepl("nFeature",colnames(reactivevalue$SeuratObject@meta.data))&!grepl("^percent.",colnames(reactivevalue$SeuratObject@meta.data))
      &!grepl("_res.",colnames(reactivevalue$SeuratObject@meta.data))
    ]
    ,selected = 'ShinyGroup')

    updateSelectizeInput(session = session,inputId = 'SampleColumn',choices=colnames(reactivevalue$SeuratObject@meta.data)[
      !grepl("nCount",colnames(reactivevalue$SeuratObject@meta.data))&!grepl("nFeature",colnames(reactivevalue$SeuratObject@meta.data))&!grepl("^percent.",colnames(reactivevalue$SeuratObject@meta.data))
      &!grepl("_res.",colnames(reactivevalue$SeuratObject@meta.data))
    ]
    ,selected = 'ShinyGroup')



    #updateSelectizeInput(session = session,inputId = 'Violingroup',choices=colnames(reactivevalue$SeuratObject@meta.data)[
    #  !grepl("nCount",colnames(reactivevalue$SeuratObject@meta.data))&!grepl("nFeature",colnames(reactivevalue$SeuratObject@meta.data))&!grepl("^percent.",colnames(reactivevalue$SeuratObject@meta.data))
    #  &!grepl("_res.",colnames(reactivevalue$SeuratObject@meta.data))
    #]
    #,selected = NULL)
    updateSelectizeInput(session = session,inputId = 'PlotGroup',choices=colnames(reactivevalue$SeuratObject@meta.data)[
      !grepl("nCount",colnames(reactivevalue$SeuratObject@meta.data))&!grepl("nFeature",colnames(reactivevalue$SeuratObject@meta.data))&!grepl("^percent.",colnames(reactivevalue$SeuratObject@meta.data))
      &!grepl("_res.",colnames(reactivevalue$SeuratObject@meta.data))
    ]
    ,selected = NULL)

    DGE_Group_Candidate=c()
    for (i in colnames(reactivevalue$SeuratObject@meta.data)[
      !grepl("nCount",colnames(reactivevalue$SeuratObject@meta.data))&!grepl("nFeature",colnames(reactivevalue$SeuratObject@meta.data))&!grepl("^percent.",colnames(reactivevalue$SeuratObject@meta.data))
      &!grepl("_res.",colnames(reactivevalue$SeuratObject@meta.data))
    ]) {
      DGE_Group_Candidate=c(DGE_Group_Candidate,paste0(i,'-',unique(as.character(reactivevalue$SeuratObject@meta.data[,i]))))
    }
    updateSelectizeInput(session = session,inputId = 'Assay',choices=names(reactivevalue$SeuratObject@assays),selected = NULL)
    updateSelectizeInput(session = session,inputId = 'DGEGroup1',choices=DGE_Group_Candidate,selected = NULL)
    updateSelectizeInput(session = session,inputId = 'DGEGroup2',choices=DGE_Group_Candidate,selected = NULL)
    print('Finish Loading')
}
})

observeEvent( input$SampleColumn, {
  if (!is.null(input$SampleColumn)&!is.null(reactivevalue$SeuratObject)){
  updateSelectInput(session = session,inputId = 'SampletoSubset',choices=unique(reactivevalue$SeuratObject@meta.data[,input$SampleColumn]),selected = NULL)

}
})


ReductionListener <- reactive({
  list(input$reduction,input$variabletogroup,input$variabletosplit,input$SampletoSubset)
})

observeEvent( ReductionListener(), {
  if (is.null(input$reduction)) return()
  if (!is.null(input$SampletoSubset)) {
    kept=rownames(reactivevalue$SeuratObject@meta.data)[!reactivevalue$SeuratObject@meta.data[,input$SampleColumn]%in%input$SampletoSubset]
    if (length(kept)!=0) {
      temp=subset(reactivevalue$SeuratObject,cells=kept)
      output$tsne=renderPlot(DimPlot(temp,reduction = input$reduction,group.by = input$variabletogroup,split.by = input$variabletosplit,ncol=ifelse(length(unique(reactivevalue$SeuratObject@meta.data[,input$variabletosplit]))==1,yes = 1,no=2)))

    } else {
      output$tsne=renderPlot(DimPlot(reactivevalue$SeuratObject,reduction = input$reduction,group.by = input$variabletogroup,split.by = input$variabletosplit,ncol=ifelse(length(unique(reactivevalue$SeuratObject@meta.data[,input$variabletosplit]))==1,yes = 1,no=2)))

    }
  } else {
    output$tsne=renderPlot(DimPlot(reactivevalue$SeuratObject,reduction = input$reduction,group.by = input$variabletogroup,split.by = input$variabletosplit,ncol=ifelse(length(unique(reactivevalue$SeuratObject@meta.data[,input$variabletosplit]))==1,yes = 1,no=2)))

  }


})


observeEvent( input$GenesToInterrogate, {
  if (is.null(input$GenesToInterrogate)) return()

  output$FeaturePlot=renderPlot(FeaturePlot(reactivevalue$SeuratObject,features = input$GenesToInterrogate,order = T))

  output$ViolinPlot=renderPlot(VlnPlot(reactivevalue$SeuratObject,assay = 'RNA',features = input$GenesToInterrogate,group.by = input$PlotGroup,
                                       pt.size = ifelse(ncol(reactivevalue$SeuratObject)>1000,yes=0,no=NULL)))
  output$RidgePlot=renderPlot(RidgePlot(reactivevalue$SeuratObject,assay = 'RNA',features = input$GenesToInterrogate,group.by = input$PlotGroup
                                       ))
})

DGEListener <- reactive({
  list(input$DGEGroup1,input$DGEGroup2)
})


observeEvent(DGEListener(), {
  if(!is.null(input$DGEGroup1)&&!is.null(input$DGEGroup2)){
    CellRanch=reactivevalue$SeuratObject@meta.data[,
      colnames(reactivevalue$SeuratObject@meta.data)[
        !grepl("nCount",colnames(reactivevalue$SeuratObject@meta.data))&!grepl("nFeature",colnames(reactivevalue$SeuratObject@meta.data))&!grepl("^percent.",colnames(reactivevalue$SeuratObject@meta.data))
        &!grepl("_res.",colnames(reactivevalue$SeuratObject@meta.data))
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

    output$GroupNumber=renderDataTable(analysis)





    }
})

observeEvent(input$submitDGE, {
  if (all(reactivevalue$analysis[,'CellNumber']!=0)) {
    FindMarkerstemp=reactivevalue$SeuratObject
    DefaultAssay(FindMarkerstemp)=input$Assay
    FindMarkerstemp=subset(FindMarkerstemp,cells=unique(c(reactivevalue$Group1Wrangled,
                                                          reactivevalue$Group2Wrangled)))
    FindMarkerstemp=NormalizeData(FindMarkerstemp)
    FindMarkerstemp@meta.data$Group=ifelse(rownames(FindMarkerstemp@meta.data)%in%reactivevalue$Group1Wrangled,
                                           yes='Group1',no='Group2')
    #for (i in 1:nrow(FindMarkerstemp@meta.data)) {
    #  if (rownames(FindMarkerstemp@meta.data)[i]%in%reactivevalue$Group1Wrangled) {
    #    FindMarkerstemp@meta.data$Group[i]='Group1'
    #  } else if (rownames(FindMarkerstemp@meta.data)[i]%in%reactivevalue$Group2Wrangled) {
    #    FindMarkerstemp@meta.data$Group[i]='Group2'
    #  }
    #}


    Results=FindMarkers(FindMarkerstemp,ident.1='Group1',ident.2='Group2',group.by='Group',
                        assay=input$Assay)
    Results$gene=rownames(Results)
    output$DifferentialExpressionAnalysisResults=renderDataTable(Results)
    reactivevalue$DifferentialExpressionAnalysisResults=Results
    Results$p_val_adj=ifelse(Results$p_val_adj==0,yes=min(Results$p_val_adj)/10,no=Results$p_val_adj)
    volcano.plot <- function(Results,fcthresh=0.25,pvalthres=0.05) {
      Results$color='grey'
      Results$color=ifelse((Results$avg_log2FC>fcthresh&Results$p_val_adj<pvalthres),
                           yes='Red',
                           no=Results$color)
      Results$color=ifelse((Results$avg_log2FC<(-1*fcthresh)&Results$p_val_adj<pvalthres),
                           yes='blue',
                           no=Results$color)
      if (length(table(Results$color))==3) {
        Results$color=factor(Results$color,levels = c('red','grey','blue'))
        volcano=ggplot(Results,aes(x=avg_log2FC,y=-log10(p_val_adj),coluor=color))+geom_point()+xlim(c(
          max(abs(Results$avg_log2FC))*(-1),
          max(abs(Results$avg_log2FC))*(1),
        ))+scale_color_manual(values = c('darkred','grey','blue'))

      } else if (length(table(Results$color))==2) {
        if ('red'%in%Results$color) {
        Results$color=factor(Results$color,levels = c('red','grey'))
        volcano=ggplot(Results,aes(x=avg_log2FC,y=-log10(p_val_adj),coluor=color))+geom_point()+xlim(c(
          max(abs(Results$avg_log2FC))*(-1),
          max(abs(Results$avg_log2FC))*(1),
        ))+scale_color_manual(values = c('darkred','grey'))

        }
        if ('blue'%in%Results$color) {
          Results$color=factor(Results$color,levels = c('blue','grey'))
          volcano=ggplot(Results,aes(x=avg_log2FC,y=-log10(p_val_adj),coluor=color))+geom_point()+xlim(c(
            max(abs(Results$avg_log2FC))*(-1),
            max(abs(Results$avg_log2FC))*(1),
          ))+scale_color_manual(values = c('blue','grey'))

        }
      } else if (length(table(Results$color))==1) {
        volcano=ggplot(Results,aes(x=avg_log2FC,y=-log10(p_val_adj),coluor=color))+geom_point()+xlim(c(
          max(abs(Results$avg_log2FC))*(-1),
          max(abs(Results$avg_log2FC))*(1),
        ))+scale_color_manual(values = c('grey'))
      }

    }
    return(volcano)
  }
  output$volcano=renderPlot((volcano.plot(Results)))

})









