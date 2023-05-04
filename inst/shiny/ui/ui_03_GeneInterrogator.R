
tabPanel("Gene Expression Interrogator",

         titlePanel("Gene Expression Interrogator"),
         sidebarLayout(
           sidebarPanel(
             selectInput('GenesToInterrogate','Enter the Gene of interest',choices = NULL,selected = NULL,multiple = T),
             selectizeInput('PlotGroup','Enter the grouping variable for Violin Plot',choices = NULL,selected = NULL)
           ),

           mainPanel(
             tabsetPanel(
               tabPanel("Feature Plot",
                         imageOutput('FeaturePlot')
),

tabPanel("Violin Plot",
         #selectizeInput('Violingroup','Enter the grouping variable for Violin Plot',choices = NULL,selected = NULL),

         imageOutput('ViolinPlot')
),
tabPanel("Ridge Plot",
         #selectizeInput('Ridgegroup','Enter the grouping variable for Violin Plot',choices = NULL,selected = NULL),

         imageOutput('RidgePlot')
)
#tabPanel("Differential Expression Analysis",
#         selectizeInput('DGEGroup','Enter the grouping variable for Differential Expression',choices = NULL,selected = NULL),
#         selectizeInput('DGEGroup1','Enter the First Group for Differential Expression',choices = NULL,selected = NULL),
#         selectizeInput('DGEGroup2','Enter the Second Group for Differential Expression',choices = NULL,selected = NULL),
#         actionButton(inputId = 'submitDGE',label = 'Submit'),
#
#         tableOutput('GeneInterrogationResults')
#),
                         )

           )
         )
)
