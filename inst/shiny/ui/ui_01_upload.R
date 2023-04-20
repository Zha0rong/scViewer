accepted = c("rds",
             "rdata",
             "h5seurat",
             ".csv")
tabPanel("Upload Data",
    useShinyjs(),
    #tags$style(appCSS),
    tags$div(
        class = "jumbotron",
        tags$div(
            class = "container",
            fluidRow(column(7, h2("SingleViewer"))),
            tags$p("Viewer and Analyzer of Single Cell Data"),
            uiOutput("tab")

        )
    ),
    # Application title
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
