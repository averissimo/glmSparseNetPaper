# #' Title
# #'
# #' @param all.pairs
# #' @param fitted.risk
# #' @param ydata
# #' @param n.cores
# #'
# #' @return
# #' @export
# c.index.fun <- function(all.pairs, fitted.risk, ydata, n.cores = params$n.cores) {
#   unlist(mclapply(seq(ncol(all.pairs)), function(ix) {
#     ix.1 <- all.pairs[1,ix]
#     ix.2 <- all.pairs[2,ix]
#     return(my.c.index.cmp(fitted.risk[ix.1], fitted.risk[ix.2],
#                           ydata[ix.1,1], ydata[ix.2,1],
#                           ydata[ix.1,2], ydata[ix.2,2]))
#   }, mc.cores = n.cores))
# }
