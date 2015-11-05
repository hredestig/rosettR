makeReport <- function(path, report, browse=TRUE) {
  allReports <- list.files(system.file("reports", package=PKG), full.names=TRUE)
  names(allReports) <- gsub(".Rmd", "", basename(allReports))
  chosen <- match.arg(report, names(allReports))
  reportDir <- file.path(path, "Output", chosen)
  if(!file.exists(reportDir))
    dir.create(reportDir)
  file.copy(allReports[chosen], reportDir, overwrite=TRUE)
  cwd <- getwd()
  on.exit(setwd(cwd))
  setwd(reportDir)
  out <- knit(basename(allReports[chosen]), envir=new.env())
  markdownToHTML(sub(".Rmd", ".md", basename(allReports[chosen])),
                 sub(".Rmd", ".html", basename(allReports[chosen])))
  if(browse)
    browseURL(file.path(reportDir, paste0(chosen, ".html")))
}

