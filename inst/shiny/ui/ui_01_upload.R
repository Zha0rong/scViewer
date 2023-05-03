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
            checkboxInput('LocalFile','Check if the file is already in the root directory and named as "Final.Analysis.rds"'),
            actionButton(inputId = 'submit',label = 'Submit')

        ),

        mainPanel(
            tabsetPanel(
                tabPanel('Main Figure of the Dataset',
                         imageOutput('MainFigure')

                         )
        )
    )
)
)
