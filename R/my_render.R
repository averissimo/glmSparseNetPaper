my.render <- function(input, output_file = NULL, ...) {
  tryCatch({
    if (!is.null(output_file)) {
      out_dir <- format(Sys.time(), 'Reports %Y%m%d')
      dir.create(out_dir, showWarnings = FALSE)
      output_file <- file.path(out_dir, output_file)
    }
    rmarkdown::render(input,
                      output_file   = output_file,
                      output_format = rmarkdown::html_document(dev='svg'),
                      ...)
  }, error = function(error) {
    futile.logger::flog.error('Something wrong happened while rendering %s',
                              error)
  }
  )
}
