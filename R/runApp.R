#' Run scViewer shiny app
#'
#' @return The shiny app will open
#'
#' @param dev Run the applicaiton in developer mode
#'
#' @examples
#' \dontrun{
#' scViewer()
#' }
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
