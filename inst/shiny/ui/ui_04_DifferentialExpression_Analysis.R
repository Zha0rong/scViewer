
tabPanel("Differential Expression Analysis",

         titlePanel("Differential Expression Analysis"),
         sidebarLayout(
           sidebarPanel(
             selectizeInput('Assay','Select the Assay to do differential expression analysis.',choices = NULL,selected = NULL),
             selectInput('DGEGroup1','Enter the First Group for Differential Expression',choices = NULL,selected = NULL,multiple = T),
             selectInput('DGEGroup2','Enter the Second Group for Differential Expression',choices = NULL,selected = NULL,multiple = T),
             actionButton(inputId = 'submitDGE',label = 'Submit')

           ),

           mainPanel(
             tabsetPanel(
               tabPanel("Differential Expression Analysis Results",
                        dataTableOutput('GroupNumber'),
                        dataTableOutput('DifferentialExpressionAnalysisResults'),
                        downloadButton("downloadData", "Download")
               )
             )

           )
         )
)
