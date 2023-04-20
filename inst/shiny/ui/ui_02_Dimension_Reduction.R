
tabPanel("Dimension Reduction Viewer",

         titlePanel("Dimension Reduction Viewer"),

         sidebarLayout(
           sidebarPanel(
             selectizeInput('reduction','Dimension Reduction to use to generate figures',choices = NULL,selected = NULL),
             selectizeInput('variabletogroup','variable used to color the plot',choices = NULL,selected = NULL),
             selectizeInput('variabletosplit','variable used to split the plot',choices = NULL,selected = NULL)
           ),

           mainPanel(


               imageOutput('tsne')

           )
         )
)
