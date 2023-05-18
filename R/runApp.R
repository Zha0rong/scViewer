#' Run scViewer shiny app
#' @import shiny
#' @import shinyjs
#' @import shinythemes
#' @import data.table
#' @import Seurat
#' @import SingleCellExperiment
#' @import ggplot2
#' @import DT

#' @return The shiny app will open
#'
#' @param filelocation The directory of sce object or seurat object. If not specific the file can be uploaded inside the app.
#'
#' @export

scViewer <- function(filelocation=NULL) {

  appDir <- system.file("shiny", package="scViewer")
  if (appDir == "") {
    stop("Could not find scViewer Try re-installing `scViewer`.",
         call. = FALSE)
  }
  shinyOptions('filelocation'=filelocation)
 shiny::runApp(appDir, display.mode="normal")
}
