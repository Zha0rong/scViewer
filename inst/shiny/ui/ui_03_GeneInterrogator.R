
tabPanel("Gene Expression Interrogator",

         titlePanel("Gene Expression Interrogator"),
         sidebarLayout(
           sidebarPanel(
             selectInput('GenesToInterrogate','Enter the Gene of interest',choices = NULL,selected = NULL,multiple = T),
             selectizeInput('GIreduction','Dimension Reduction to use to generate figures',choices = NULL,selected = NULL),
             selectizeInput('PlotGroup','Enter the grouping variable for Violin Plot',choices = NULL,selected = NULL),
             checkboxInput('GILabel','Whether to Label your Data')

           ),

           mainPanel(
             tabsetPanel(
               tabPanel("Feature Plot",
                         imageOutput('FeaturePlot')
),

tabPanel("Violin Plot",

         imageOutput('ViolinPlot')
),
tabPanel("Ridge Plot",

         imageOutput('RidgePlot')
),
tabPanel('Global Expression Statistics',
         dataTableOutput('GlobalStats')
         )

                         )

           )
         )
)
