#' Launch Shiny App for DiffCoExpr
#'
#' A function that launches the Shiny app for DiffCoExpr. The purpose of this 
#' app is only to demonstrate how DiffCoExpr can be used to generate and 
#' compare coexpression networks. 
#' 
#' The code has been placed in \code{./inst/shiny-scripts}.
#'
#' @return No return value but open up a Shiny page.
#'
#' @examples
#' \dontrun{
#'
#' DiffCoExpr::runDiffCoExpr()
#' }
#'
#' @references
#' Grolemund, G. (2015). Learn Shiny - Video Tutorials. 
#' \href{https://shiny.rstudio.com/tutorial/}{Link}
#'
#' @export
#' @importFrom shiny runApp
runDiffCoExpr <- function() {
  
  appDir <- system.file("shiny-scripts",
                        package = "DiffCoExpr")
  actionShiny <- shiny::runApp(appDir, display.mode = "normal")
  
  return(actionShiny)
}

# [END]