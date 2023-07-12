
tabPanel("scType Annotation",
         
         titlePanel("SingleR Annotation"),
         sidebarLayout(
           sidebarPanel(
             actionButton(inputId = 'submitReference',label = 'Submit')
             
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
