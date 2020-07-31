
#' Draw Venn diagram with 3 sets
#'
#' @param sets
#' @param main
#' @param cat.pos
#'
#' @return
#' @export
#'
#' @examples
#' draw.3venn(sets = list(a = c(1,2,3,4), b = c(4,5,6,7), d = c(8,9,10)))
draw.3venn <- function(sets, main = NULL, cat.pos = c(-10, 10, 0)) {
  gridExtra::grid.arrange(
    grid::grid.grabExpr(
      grid::grid.draw(
        VennDiagram::venn.diagram(
          x          = sets,
          main       = main,
          filename   = NULL,
          fontface   = "plain",
          fontfamily = "sans",
          # alpha      = 0.2,
          cat.cex    = 1.5,
          cex        = 1.4,
          margin     = c(0,0,0,0),
          cat.fontface = "plain",
          cat.default.pos = "outer",
          cat.pos         = cat.pos,
          #cat.dist        = c(0.055, 0.055, 0.055),
          cat.fontfamily  = "sans",
          rotation        = 1,
          fill            = c("#BEBADA", "#FB8072", "#8DD3C7"),
          cat.col         = c("#56B353", "#E41A1C", "#377EB8"),
          col             = c("#56B353", "#E41A1C", "#377EB8"),
          log.to.file     = FALSE
        )
      )
    )
  )
}
