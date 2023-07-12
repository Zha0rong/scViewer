
tabPanel("CellFindR",
         
         titlePanel("CellFindR parameter selection"),
         sidebarLayout(
           sidebarPanel(
             actionButton(inputId = 'CellFindR',label = 'Submit')
             
           ),
           
           mainPanel(
             tabsetPanel(
               tabPanel("Reference Information",
                        DT::dataTableOutput('scType_ReferenceInformation')
               ),
               tabPanel("Annotation Results",
                        DT::dataTableOutput('AnnotationResults'),
                        plotOutput('scType_visualized.AnnotationResults',height='800px')
               )
             )
             
           )
         )
)
