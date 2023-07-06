accepted = c(".rds",
             ".rdata",
             ".h5seurat",
             ".csv")
tabPanel("Upload Data",useShinyjs(),
    tags$div(class = "jumbotron",tags$div(
            class = "container",
            fluidRow(column(7, h2("SingleViewer"))),
            tags$p("Viewer and Analyzer of Single Cell Data"),
            uiOutput("tab"))),
        mainPanel(
          tabsetPanel(
            tabPanel('Main Figure of the Dataset',
                     imageOutput('MainFigure')
                     
            ),
            tabPanel('Data Distribution',
                     selectizeInput('BarGraph1','First Variable to interrogate distribution',choices=NULL,selected=NULL),
                     selectizeInput('BarGraph2','Second Variable to interrogate distribution',choices=NULL,selected=NULL),
                     checkboxInput('BarGraphPercentage','Transfer the BarGraph into percentage plot'),
                     plotOutput('BarPlot'),
                     DT::dataTableOutput('BarPlotStats')
                     
            ),
            tabPanel('Annotate Data',
                     selectizeInput('Reference_Column','Reference Column',choices=NULL,selected=NULL),
                     textInput('AnnotationName',label = 'Name of new column',value = NULL),
                     DT::dataTableOutput('AnnotationTable'),
                     actionButton('AddAnnotation',label = 'Add Annotation')
                     
            )
          )
        )

)
