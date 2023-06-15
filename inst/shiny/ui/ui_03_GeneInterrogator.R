
tabPanel("Gene Expression Interrogator",

         titlePanel("Gene Expression Interrogator"),
         sidebarLayout(
           sidebarPanel(
             selectizeInput('GenesToInterrogateAssay','Select the Assay to explore expression.',choices = NULL,selected = NULL),
             selectInput('GenesToInterrogate','Enter the Gene of interest',choices = NULL,selected = NULL,multiple = T),
             selectizeInput('GIreduction','Dimension Reduction to use to generate figures',choices = NULL,selected = NULL),
             selectizeInput('PlotGroup','Enter the grouping variable for Violin Plot',choices = NULL,selected = NULL),
             checkboxInput('GILabel','Whether to Label your Data')

           ),

           mainPanel(
             tabsetPanel(
               tabPanel("Feature Plot",
                        checkboxInput('FeaturePlotOverlay','Generate Overlay Plot for multiple Genes'),
                         imageOutput('FeaturePlot')
),

tabPanel("Violin Plot",

         imageOutput('ViolinPlot')
),
tabPanel("Ridge Plot",

         imageOutput('RidgePlot')
),
tabPanel('Global Expression Statistics',
         DT::dataTableOutput('GlobalStats'),
         downloadButton('Gene.Expression.Statistics.downloadData','Download')
         ),
tabPanel('Cell Statistics based on gene filtering',
         h3('In this tab the number of cells with non-zero expressions of selected genes per group in the selected variable will be displayed.'),
         DT::dataTableOutput('CellStatsTableBasedonGeneFiltering')
)

                         )

           )
         )
)
