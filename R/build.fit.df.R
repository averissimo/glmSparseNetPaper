build.fit.df <- function(obj, xdata.ix) {
  fitted.network <- data.frame()

  for (ix.name in names(obj)) {
    risk.a <- predict(obj[[ix.name]]$model, newx = xdata.train[,xdata.ix], s = obj[[ix.name]]$target, type = 'response')
    risk.b <- predict(obj[[ix.name]]$model, newx = xdata.test[,xdata.ix], s = obj[[ix.name]]$target, type = 'response')

    a <- data.frame(relative.risk = as.vector(risk.a), set = 'Train', type = ix.name, stringsAsFactors = FALSE)
    b <- data.frame(relative.risk = as.vector(risk.b), set = 'Test', type = ix.name, stringsAsFactors = FALSE)

    a$mean <- mean(a$relative.risk)
    b$mean <- mean(b$relative.risk)

    fitted.network <- rbind(fitted.network, a, b)
  }
  return(fitted.network)
}
