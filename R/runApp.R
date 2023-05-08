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
  wd=getwd()
  filelocation=list.files(wd,full.names = T,pattern = '.rds',)
  if (length(filelocation)==0) {
    print('Empty')
    filelocation=''
  }
  shinyOptions('filelocation'=filelocation)
 shiny::runApp(appDir, display.mode="normal")
}
