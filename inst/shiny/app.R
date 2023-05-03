require(shiny)
require(shinyjs)
require(shinythemes)
require(data.table)
require(Seurat)

ui <- navbarPage(

  title = "SingleViewer",
  id="SingleViewer",
  fluid=TRUE,
  theme = shinytheme("yeti"),
  source(file.path("ui", "ui_01_upload.R"),  local = TRUE)$value,
  source(file.path("ui", "ui_02_Dimension_Reduction.R"),  local = TRUE)$value,
  source(file.path("ui", "ui_03_GeneInterrogator.R"),  local = TRUE)$value,
  source(file.path("ui", "ui_04_DifferentialExpression_Analysis.R"),  local = TRUE)$value
  #source(file.path("ui", "ui_05_heatmaps.R"),  local = TRUE)$value,
  #source(file.path("ui", "ui_06_circular_dendogram.R"),  local = TRUE)$value,
  #source(file.path("ui", "ui_07_pca.R"),  local = TRUE)$value,
  #source(file.path("ui", "ui_08_shape.R"),  local = TRUE)$value
)

server <- function(input, output, session) {
  source(file.path("server/", "server.R"),  local = TRUE)$value
}

shinyApp(ui = ui, server = server)
