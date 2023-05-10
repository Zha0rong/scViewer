accepted = c(".rds",
             ".rdata",
             ".h5seurat",
             ".csv")
tabPanel("Upload Data",
    useShinyjs(),
    tags$div(
        class = "jumbotron",
        tags$div(
            class = "container",
            fluidRow(column(7, h2("SingleViewer"))),
            tags$p("Viewer and Analyzer of Single Cell Data"),
            uiOutput("tab")

        )
    ),
    titlePanel("Upload Data"),

    sidebarLayout(
        sidebarPanel(
            h3("Upload Seurat Object"),
            fileInput(
                inputId = "Seurat_Object",
                label = "Seurat Object",
                multiple = FALSE,
                accept = accepted
            ),
            #selectizeInput('objecttype','SingleCellExperiment object or Seurat object',choices=c('Seurat','SingleCellExperiment'),selected='SingleCellExperiment'),
            verbatimTextOutput('object_location'),
            actionButton(inputId = 'submit',label = 'Submit')

        ),

        mainPanel(
            tabsetPanel(
                tabPanel('Main Figure of the Dataset',
                         imageOutput('MainFigure')

                         ),
                tabPanel('Data Distribution BarGraph',
                         selectizeInput('BarGraph1','First Variable to interrogate distribution',choices=NULL,selected=NULL),
                         selectizeInput('BarGraph2','Second Variable to interrogate distribution',choices=NULL,selected=NULL),
                         plotOutput('BarPlot')

                )
        )
    )
)
)
