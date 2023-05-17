
tabPanel("Dimension Reduction Viewer",

         titlePanel("Dimension Reduction Viewer"),

         sidebarLayout(
           sidebarPanel(
             selectizeInput('DRreduction','Dimension Reduction to use to generate figures',choices = NULL,selected = NULL),
             selectizeInput('variabletogroup','variable used to color the plot',choices = NULL,selected = NULL),
             selectizeInput('variabletosplit','variable used to split the plot',choices = NULL,selected = NULL),
             selectizeInput('SampleColumn','Specify based on which variable do you want to subset the data.',choices = NULL,selected = NULL),
             selectInput('SampletoSubset','Specify which sample(s) do you want to remove.',choices = NULL,selected = NULL,multiple = T),
             checkboxInput('DRLabel','Whether to Label your Data')

           ),

           mainPanel(


               imageOutput('tsne')

           )
         )
)
