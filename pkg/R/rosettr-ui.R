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
#' meta <- metaTemplate(letters[1:4], c("control", "osmotic_stress"), reference="a")
#' newExperiment(file.path(tempdir(), "test"), meta)
#' list.files(tempdir())
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
#' @return invisibly, the path to the created test experiment
#' @export
#' @examples
#' makeTestExperiment(tempdir())
#' @author Henning Redestig
makeTestExperiment <- function(exdir=".") {
  path <- file.path(exdir, "rosettrTest")
  unzip(system.file("examples/rosettrTest.zip", package=PKG),
        exdir=exdir)
  meta <- metaTemplate(c("foo", "bar", "baz", "qux"),
                       treatments=c("control", "osmotic"),
                       timepoints=c(11, 14, 16, 18),
                       pixelsmm=7.538, nblocks=3,
                       reference="foo")
  writeMeta(meta, path)
  invisible(path)
}

#' Create a meta data template
#'
#' Each experiment is annotated with a set of parameters stored in a
#' 'meta' data object. Parameters to use may differ between experiment
#' setup and design. In all cases, the basic experiment setup is to
#' have plates with seedlings that are taken from the growth chamber
#' for imaging on a pre-defined set of days-after-sowing.
#' @details This package is delivered with three different experiment
#' designs. These are:
#'
#' \itemize{ \item{6x6.abcd}{Uses a 6x6 plate grid with genotypes
#' (lines) allocated to regions a, b, c and d that are placed
#' clockwise around the plate. Total number of genotypes must thus be
#' divisible by four.}
#' 
#' \item{6x6.ab}{Uses a 6x6 plate grid with genotypes (lines)
#' allocated to regions a and b that positioned at the upper and lower
#' half of the plate.}
#' }
#' @param genotypes the names of the genotypes used in the experiment.
#' @param treatments a character vector listing the applied treatments
#' @param timepoints the timepoints (days) when pictures are taken.
#' @param nblocks the number of replicates
#' @param description a short description of the experiment to
#' annotate template reports
#' @param pixelsmm the number of pixels per mm on a picture (must be
#' the same for all pictures)
#' @param name the name of the experiment type to use. See details.
#' @param plateRadius the radius of the plate (max radius of the lid)
#' @param nBoxGrid the number of wells in each row and column of a
#' plate. Layout is assumed to be square but with the corner wells
#' left-out.
#' @param boxWidth the width of a well in millimeter
#' @param reference a list of genotypes that is meant to be used as the
#' reference. Necessary for the 'compare areas' report but can be left
#' \code{NULL} if not applicable.
#' @return a list with meta data for the experiment
#' @export
#' @examples
#' metaTemplate(c("foo", "bar", "baz", "qux"),
#'              c("control", "osmotic_stress"), reference="foo")
metaTemplate <- function(genotypes,
                         treatments="control",
                         timepoints=c(11, 14, 16, 18),
                         nblocks=10,
                         description="",
                         pixelsmm=16,
                         name=c("6x6.abcd", "6x6.ab"),
                         plateRadius=76, nBoxGrid=6, boxWidth=20,
                         reference) {
  name <- match.arg(name)
  if(!is.null(reference))
    stopifnot(reference %in% genotypes)
  template <- list(genotypes=genotypes, treatments=treatments,
                   nBoxGrid=nBoxGrid, boxWidth=boxWidth,
                   timepoints=timepoints, nblocks=nblocks,
                   pixelsmm=pixelsmm, name=name, plateRadius=plateRadius,
                   reference=reference)
  
  reps <- expand.grid(RANGE=1:6, ROW=1:6)
  preGriddf <- data.frame(RANGE=reps$RANGE, ROW=reps$ROW)
  griddf <- preGriddf[!(preGriddf$ROW %in% c(1, 6) &
                          preGriddf$RANGE %in% c(1, 6)),]
  rownames(griddf) <- 1:32
  griddf$box_num <- 1:nrow(griddf)
  griddf$removed <- FALSE
  griddf$too_small <- FALSE

  if(name == "6x6.abcd") {
    griddf$genotype_region <-
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
    griddf$genotype_region <-
      c(     "b", "b", "a", "a", 
        "b", "b", "b", "a", "a", "a", 
        "b", "b", "b", "a", "a", "a", 
        "b", "b", "b", "a", "a", "a", 
        "b", "b", "b", "a", "a", "a", 
             "b", "b", "a", "a"      )
    griddf$Sample_ID <- c(1:16, 1:16)
  }

  ## if(name == "6x6.one") {
  ##   griddf$Sample_ID <- 32:1
  ##   griddf$genotype_region <- "a"
  ## }
  template$griddf <- griddf
  template
}
