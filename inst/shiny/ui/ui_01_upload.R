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
    h2("Upload Data",id='UploadDataTitle'),
    #titlePanel("Upload Data")
    sidebarLayout(sidebarPanel(id='UploadBar',
            h3("Upload Seurat Object"),
            fileInput(
                inputId = "Seurat_Object",
                label = "Seurat Object",
                multiple = FALSE,
                accept = accepted
            ),
            #selectizeInput('objecttype','SingleCellExperiment object or Seurat object',choices=c('Seurat','SingleCellExperiment'),selected='SingleCellExperiment'),
            verbatimTextOutput('object_location'),
            actionButton(inputId = 'submit',label = 'Submit'),width = 3

        ),
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

#placeholcer
)

)
