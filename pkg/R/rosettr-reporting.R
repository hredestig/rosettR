#' Compile a template report for a plate phenotyping experiment
#'
#' This package comes with a set of pre-defined reports for facilitating
#' preparing and evaluating your in-vitro plate phenotyping experiments. This
#' function allows to conveniently compile reports for created experiments.
#' @param path path on the file system to where the experiment directory is located
#' @param report the report to compile, see \code{\link{plate-reports}}
#' @param browse if true, open the compiled report in a web-browser
#' @return nothing, used for its side effect
#' @export
#' @examples
#' meta <- metaTemplate(letters[1:4], treatments=c("control", "osmotic"))
#' newExperiment(file.path(tempdir(), "testExperiment"), meta)
#' makeReport(file.path(tempdir(), "testExperiment"), "sowing")
#' @author Henning Redestig
makeReport <- function(path, report, browse=interactive()) {
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


#' Make a thumbnail gallery of the input images 
#'
#' Creates a html thumbnail gallery of the input images where the
#' images are arranged after how the will be interpreted by the
#' renaming facilities and which day they were taken. The images are
#' clickable and links to the original input images.
#' @param path a resolvable path to the experiment.
#' @param what either 'qc' or 'raw' for galleries of qc images or raw
#' images respectively.
#' @return nothing, used for its side effect.
#' @export
#' @examples
#' \dontrun{
#' dummyExperimentId <- createDummyPlateExperiment(242)
#' plateGallery(dummyExperimentId, "raw")
#' }
#' @author Henning Redestig
plateGallery <- function(path, what=c("raw", "qc")) {
  what <- match.arg(what)
  switch(what, qc={
    pda <- readPhenodata(path)
    cols <- c("plate", "timepoint", "image", "qc_picture")
    if(!all(cols %in% names(pda))) {
      message("missing qc pictures or images in phenodata, skipping thumbnails")
      return(NULL)
    }
    df <- unique(pda[,cols])
    df$orig <- file.path('..', df$image)
    df$thumb <- file.path('..', df$qc_picture)
    plateMakeTng(df, 'qc-thumbnails.html', path)
  }, raw={
    version <- system("convert --version", intern=TRUE)
    if(!any(grepl("ImageMagick", version)))
      stop("ImageMagick's 'convert' utility not available")
    mf <- readManifest(path)
    rnm_df <- ddply(mf, "timepoint", function(dd) {
      daydir <- unique(dirname(as.character(dd$image)))
      dd <- dd[with(dd, order(BLOCK, position)),]
      first_region <- unique(dd$germplasm_region)[1]
      expected_pics <-
        basename(as.character(subset(dd, dd$germplasm_region ==
                                       first_region)$image)) 
      renaming_df(cleanPath(file.path(path, daydir), mustWork=TRUE),
                  expected_pics)
    })
    for(s in unique(rnm_df$subdir)) {
      pi <- file.path(path, s)
      pt <- file.path(path, "Output", "thumbs", s)
      if(!file.exists(pt)) dir.create(pt, recursive=TRUE)
      system.time(res <- llply(rnm_df$image[rnm_df$subdir == s], function(im) {
        system2("convert", c("-define jpeg:size=1000x500 -thumbnail 400x200",
                             sprintf('"%s/%s"', pi, im),
                             sprintf('"%s/%s"', pt, im)))
      }, .progress=ifelse(interactive(), "text", "none")))
    }
    df <- data.frame(orig=file.path("..", rnm_df$subdir, rnm_df$image),
                     thumb=file.path("thumbs", rnm_df$subdir, rnm_df$image),
                     plate=rnm_df$newname,
                     timepoint=as.numeric(gsub(".*_*D(\\d+)", "\\1",
                       rnm_df$subdir)))
    plateMakeTng(df, "raw-thumbnails.html", path)
  })
}

plateMakeTng <- function(df, output, path) {
  df <- with(df, df[order(plate, orig),])
  df$link <- paste('<a href="', file.path(df$orig), '">',
                   '<img src="', file.path(df$thumb), '", rel="lightbox"></a>',
                   sep="")
  cdf <- reshape2::dcast(df, plate ~ timepoint, value.var="link")
  rownames(df) <- NULL
  if(nrow(df) == 0)
    return(NULL)
  html_out <- file.path(path, "Output", output)
  cat( "<!DOCTYPE html>"
      ,"<html> <head>"
      ,'<meta http-equiv="Content-Type" content="text/html; charset=utf-8"/>'
      ,"<title>Plates for experiment", basename(path), "</title></head>"
      ,"<body>"
      ,"<h1>Plates for experiment", basename(path), "</h1>", file=html_out)
  tab <- print(xtable::xtable(cdf), "html", sanitize.text.function=identity,
               file=html_out, append=TRUE)
  cat("</body></html>", file=html_out, append=TRUE)
}

