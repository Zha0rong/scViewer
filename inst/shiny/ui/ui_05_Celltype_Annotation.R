
tabPanel("SingleR Annotation",
         
         titlePanel("SingleR Annotation"),
         sidebarLayout(
           sidebarPanel(
             selectInput('Annotate_Group','Enter Which Group is used to annotate',choices = NULL,selected = NULL,multiple = F),
             selectInput('Reference','Choose the dataset used as reference',choices = c(
               'HumanPrimaryCellAtlasData',
               'BlueprintEncodeData',
               'MouseRNAseqData',
               'ImmGenData',
               'DatabaseImmuneCellExpressionData',
               'NovershternHematopoieticData',
               
               'MonacoImmuneData'
               
               
             ),selected = 'HumanPrimaryCellAtlasData',multiple = F),
             selectInput('levelofannotation','Select the level of annotatioon',choices = c(
               'Main',
               'Fine'
               
               
             ),selected = 'Main',multiple = F),
             
             actionButton(inputId = 'submitAnnotation',label = 'Submit')
             
           ),
           
           mainPanel(
             tabsetPanel(
               tabPanel("Reference Information",
                        DT::dataTableOutput('ReferenceInformation')
               ),
               tabPanel("Annotation Results",
                        DT::dataTableOutput('AnnotationResults'),
                        plotOutput('visualized.AnnotationResults',height='800px')
               )
             )
             
           )
         )
)
