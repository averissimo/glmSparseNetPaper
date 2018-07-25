#' Test all survival variables that are continuous or binary
#'
#' @param clinical clinical dataset from tcga.data package
#' @param ydata table with time and status variables
#'
#' @return
#' @export
#' @example
#' test.all.clinical(clinical, ydata)
test.all.clinical <- function(clinical, ydata) {
  df        <- data.frame()
  plot.list <- list()
  plot.p    <- c()

  #
  # Iterate on all variables
  for (ix.name in colnames(clinical$all)) {

    # skip days_to_death
    if (ix.name %in% c('days_to_death', 'days_to_last_followup')) {
      next
    }

    #
    # If variables has 2 values then make it a factor
    u <- length(unique(clinical$all[[ix.name]]))
    if (is.numeric(clinical$all[[ix.name]]) || u == 2) {

      #
      # Test scenario

      test.me        <- clinical$all[strtrim(rownames(ydata), 12), ix.name]
      names(test.me) <- rownames(ydata)

      if (!is.factor(test.me) && all(sapply(test.me, is.character))) {
        test.me <- factor(test.me)
      }

      #
      # filter out NA and NaN

      test.me <- test.me[!(is.nan(test.me) | is.na(test.me))]

      if (length(test.me) == 0) { next }

      test.me.df     <- data.frame(time    = ydata[names(test.me), 'time'],
                                   status  = ydata[names(test.me), 'status'],
                                   test.me = test.me)

      colnames(test.me.df)[3] <- ix.name
      # build model with a single variable
      test.me.model <- coxph(Surv(time, status) ~ ., data = test.me.df, control = coxph.control(iter.max = 500))
      p.value       <- summary(test.me.model)$coefficients[1,'Pr(>|z|)']
      #
      # skip if model is not right
      if (is.na(p.value)) { next }
      desc.values <- c(0,0)
      if (is.factor(test.me)) {
        desc.values <- paste(levels(test.me), collapse = ', ')
        my.names        <- names(test.me)
        levels(test.me) <- c(0,1)
        test.me         <- as.numeric(as.character(test.me))
        names(test.me)  <- my.names
      } else {
        desc.values <- sprintf('min: %g -- max %g -- median %g', min(test.me), max(test.me), median(test.me))
      }
      #
      df <- rbind(df, data.frame(name    = ix.name,
                                 numeric = is.numeric(clinical$all[[ix.name]]),
                                 unique  = u,
                                 p.value = p.value,
                                 desc    = desc.values))
      #
      my.coef        <- coef(test.me.model)
      names(my.coef) <- 'test.me'
      ydata.me       <- data.frame(time = ydata[names(test.me),1], status = ydata[names(test.me),2] * 1)
      test.me.km     <- draw.kaplan(as.vector(my.coef), as.vector(test.me), ydata.me, plot.title = ix.name, legend.outside = F)
      # build list to multiplot
      #  but skip if not significant
      if (p.value > 0.05 && test.me.km$pvalue > 0.05) { next }
      plot.list[[length(plot.list) + 1]] <- test.me.km$plot
      plot.p <- c(plot.p, test.me.km$pvalue)
    }
  }
  ncol       <- 3
  p.value.ix <- sort(plot.p, index.return = T)$ix # used to sort layout by p.values
  layout     <- matrix(c((1:length(plot.list))[p.value.ix], rep(NA, (ncol * ceiling(length(plot.list) / ncol)) - (length(plot.list)))),
                       byrow = T, ncol = ncol)
  return(list(plot.list = plot.list, layout = layout, df = df))
}
