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
#' Creates a html thumbnailDir gallery of the input images where the
#' images are arranged after how the will be interpreted by the
#' renaming facilities and which day they were taken. The images are
#' clickable and links to the original input images.
#' @param path a resolvable path to the experiment.
#' @param what either 'qc' or 'raw' for galleries of qc images or raw
#' images respectively.
#' @param parallel use registered parallel backend for faster
#' generation of image thumbnails, see \code{\link{llply}}
#' @return nothing, used for its side effect.
#' @export
#' @seealso reports
#' @examples
#' \dontrun{
#' plateGallery(pathToExperiment, "raw")
#' }
#' @author Henning Redestig
plateGallery <- function(path, what=c("raw", "qc"), parallel=FALSE) {
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
    mf <- readManifest(path)
    rnmDf <- ddply(mf, "timepoint", function(dd) {
      daydir <- unique(dirname(as.character(dd$image)))
      dd <- dd[with(dd, order(BLOCK, position)),]
      first_region <- unique(dd$germplasm_region)[1]
      expectedPics <-
        basename(as.character(subset(dd, dd$germplasm_region ==
                                       first_region)$image)) 
      renamingDf(cleanPath(file.path(path, daydir), mustWork=TRUE),
                  expectedPics)
    })
    for(dayDir in unique(rnmDf$subdir)) {
      imSubDir <- file.path(path, dayDir)
      thumbnailDir <- file.path(path, "Output", "thumbs", dayDir)
      if(!file.exists(thumbnailDir)) dir.create(thumbnailDir, recursive=TRUE)
      res <- llply(rnmDf$image[rnmDf$subdir == dayDir], function(imFile) {
        if(!file.exists(file.path(thumbnailDir, imFile))) {
          image <- readImage(file.path(path, dayDir, imFile))
          thumb <- resize(image, w=300)
          writeImage(thumb, file=file.path(thumbnailDir, imFile))
        }
      }, .parallel=parallel)
    }
    df <- data.frame(orig=file.path("../..", rnmDf$subdir, rnmDf$image),
                     thumb=file.path("thumbs", rnmDf$subdir, rnmDf$image),
                     plate=rnmDf$newname,
                     timepoint=as.numeric(gsub(".*_*D(\\d+)", "\\1",
                       rnmDf$subdir)))
    plateMakeTng(df)
  })
}

plateMakeTng <- function(df) {
  df <- with(df, df[order(plate, orig),])
  df$link <- paste('<a href="', file.path(df$orig), '">',
                   '<img src="', "../",
                   file.path(df$thumb), '", rel="lightbox"></a>',
                   sep="")
  cdf <- reshape2::dcast(df, plate ~ timepoint, value.var="link")
  rownames(df) <- NULL
  if(nrow(df) == 0)
    return(NULL)
  print(xtable::xtable(cdf), "html", sanitize.text.function=identity)
}

