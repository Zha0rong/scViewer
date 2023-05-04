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

scViewer <- function(dev=FALSE) {
  appDir <- system.file("shiny", package="scViewer")
  if (appDir == "") {
    stop("Could not find scViewer Try re-installing `scViewer`.",
         call. = FALSE)
  }
  if (dev) {
    options(shiny.autoreload=TRUE)
  }
  shiny::runApp(appDir, display.mode="normal")
  print(getwd())
}
