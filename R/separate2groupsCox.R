separate2GroupsCox <- function(...) {
  glmSparseNet::separate2GroupsCox(...,
                                   risk.table = TRUE,
                                   cumcensor = TRUE,
                                   xlab = "Time in months",
                                   surv.median.line = "hv")
}
