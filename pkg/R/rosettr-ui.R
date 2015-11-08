#' Start a new plate phenotyping experiment
#'
#' The experiment designs implemented in this package are managed in
#' file system directories with pre-defined structure. Use this
#' function to create a new experiment with sowing instructions and
#' directories in which to place pictures.
#' @param path the desired path for the experiment
#' @param meta a meta data object as created by \code{\link{metaTemplate}}
#' @return nothing, used for the side effect of create a new experiment directory
#' @export
#' @examples
#' meta <- metaTemplate(letters[1:4], c("control", "osmotic_stress"))
#' newExperiment(file.path(tempdir(), "test"), meta)
#' list.dir(tempdir())
newExperiment <- function(path, meta) {
  if(!file.exists(path))
    dir.create(path, recursive=TRUE)
  else if(length(list.files(path)) != 0)
    stop("directory exists but is not empty")
  dir.create(file.path(path, "Output"))
  lapply(meta$timepoints, function(day)
    dir.create(file.path(path, sprintf("D%02d", day))))
  writeManifest(expandManifest(meta), path)
  writeMeta(meta, path)
}

#' Create a dummy plate experiment
#'
#' To be used for testing purposed in documentation and unit-tests. 
#' @param exdir where to create the test experiment
#' @return nothing, used for side effect
#' @export
#' @examples
#' makeTestExperiment(tempdir())
#' @author Henning Redestig
makeTestExperiment <- function(exdir=".") {
  path <- file.path(exdir, "rosettrTest")
  unzip(system.file("examples/rosettrTest.zip", package=PKG))
  meta <- metaTemplate(letters[1:4], treatments=c("control", "osmotic"),
                       timepoints=c(11, 14, 16, 18), pixelsmm=7.538, nrepeats=3)
  writeManifest(expandManifest(meta), path)
  writeMeta(meta, path)
}

#' Create a meta data template
#'
#' Each experiment is annotated with a set of parameters stored in a 'meta' data
#' object. Parameters to use may differ between experiment setup and design. In
#' all cases, the basic experiment setup is to have plates with seedlings that
#' are taken from the growth chamber for imaging on a pre-defined set of
#' days-after-sowing. 
#' @details This package is delivered with three different experiment
#' designs. These are:
#'
#' \itemize{
#' \item{6x6.abcd}{Uses a 6x6 plate grid with germplasms (lines) allocated to
#' regions a, b, c and d that are placed clockwise around the plate. Total
#' number of germplasms must thus be divisible by four.}
#' 
#' \item{6x6.ab}{Uses a 6x6 plate grid with germplasms (lines) allocated to
#' regions a and b that positioned at the upper and lower half of the plate.}
#'
#' \item{6x6.one}{Uses a 6x6 plate grid with a single germplasm.}
#' }
#' @param name the name of the experiment type to use. See details.
#' @param germplasms the names of the germplasms (genotypes) used in the
#' experiment.
#' @param timepoints the timepoints (days) when pictures are taken.
#' @param nrepeats the number of replicates
#' @param pixelsmm the number of pixels per mm on a picture (must be the same
#' for all pictures)
#' @param plateRadius the radius of the plate (max radius of the lid)
#' @param r the number of wells in each row and column of a plate. Layout is
#' assumed to be square but with the corner wells left-out.
#' @param d the width of a well in millimeter
#' @return a list with meta data for the experiment
#' @export
#' @examples
#' metaTemplate(letters[1:4], c("control", "osmotic_stress"))
#' @author Henning Redestig
metaTemplate <- function(germplasms,
                         treatments="control",
                         timepoints=c(11, 14, 16, 18),
                         nrepeats=10,
                         pixelsmm=16,
                         name=c("6x6.abcd", "6x6.ab", "6x6.one"),
                         plateRadius=76, nBoxGrid=6, boxWidth=20) {
  name <- match.arg(name)
  template <- list(germplasms=germplasms, treatments=treatments,
                   nBoxGrid=nBoxGrid, boxWidth=boxWidth,
                   timepoints=timepoints, nrepeats=nrepeats,
                   pixelsmm=pixelsmm, name=name, plateRadius=plateRadius)
  
  reps <- expand.grid(RANGE=1:6, ROW=1:6)
  griddf <- subset(data.frame(RANGE=reps$RANGE, ROW=reps$ROW), 
                   !(ROW %in% c(1, 6) & RANGE %in% c(1, 6)))
  rownames(griddf) <- 1:32
  griddf$box_num <- 1:nrow(griddf)
  griddf$removed <- FALSE
  griddf$too_small <- FALSE

  if(name == "6x6.abcd") {
    griddf$germplasm_region <-
        c(     "c", "c", "d", "d", 
          "c", "c", "c", "d", "d", "d", 
          "c", "c", "c", "d", "d", "d", 
          "b", "b", "b", "a", "a", "a", 
          "b", "b", "b", "a", "a", "a", 
               "b", "b", "a", "a"      )
    griddf$Sample_ID <-
      c(   1, 2, 1, 2, 
        3, 4, 5, 3, 4, 5, 
        6, 7, 8, 6, 7, 8, 
        1, 2, 3, 1, 2, 3, 
        4, 5, 6, 4, 5, 6, 
           7, 8, 7, 8   )
  }
  
  if(name == "6x6.ab") {
    griddf$germplasm_region <-
      c(     "b", "b", "a", "a", 
        "b", "b", "b", "a", "a", "a", 
        "b", "b", "b", "a", "a", "a", 
        "b", "b", "b", "a", "a", "a", 
        "b", "b", "b", "a", "a", "a", 
             "b", "b", "a", "a"      )
    griddf$Sample_ID <- c(1:16, 1:16)
  }

  if(name == "6x6.one") {
    griddf$Sample_ID <- 32:1
    griddf$germplasm_region <- "a"
  }
  template$griddf <- griddf
  template
}
