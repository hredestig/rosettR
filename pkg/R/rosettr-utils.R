#' Update the meta object
#'
#' Meta object is a list of parameters to use as defaults for a given
#' experiment. One can tweak these defaults for e.g. specific image by passing
#' parametesr as named arguments. This function then updates the meta object
#' with these tweaks
#' @param ... named arguments
#' @param meta the meta object
#' @return the updated meta object
#' @author Henning Redestig
updateMeta <- function(..., meta) {
  newArgs <- list(...)
  for(item in names(newArgs))
    meta[[item]] <- newArgs[[item]]
  meta
}

readPhenodata <- function(path)
  readRDS(file.path(path, "Output", "phenodata.rda"))

writePhenodata <- function(path, df) {
  saveRDS(df, file.path(path, "Output", "phenodata.rda"))
  write.csv(df, file=file.path(path, "Output", "data.csv"))
}

cleanPath <- function (path, mustWork=NA) {
  path <- gsub("/+$", "", normalizePath(path, winslash="/", 
                                        mustWork=mustWork))
  gsub("\\\\+$", "", path)
}

#' Quality control statistics for plate images
#'
#' Creates a data frame with QC statistics from the result of the
#' image processing pipeline
#'
#' The following components are given \describe{
#' 
#' \item{rotation}{the clockwise rotation in degrees that was applied
#' to the plate. Missing value indicates that no rotation was
#' done. An absolute value close to the maximum tested rotation may
#' indicate poor quality}
#' 
#' \item{eccentricity}{an estimate for how far away the plate
#' deviates from the center in millimeter. A large value may indicate
#' poor quality but the plate detection algorithm is unstable an
#' often over-estimates this value.}
#' 
#' \item{ambiguous\_boxes}{the number of boxes on the plate that had
#' features that could not be resolved to a specific boxes and hence
#' treated by hard-splitting of the box and summing features. A
#' number greater than 0 indicates poor quality.}
#' 
#' \item{max_features}{the maximum number of features in a box. The
#' larger the number the worse quality.}
#'
#' }
#' @param df a data frame with all data obtained from a plate
#' experiment (a \code{phenodata} dataset).
#' @param html should the dataframe be formatted to HTML or not
#' @return A data frame with quality control statistics ready for
#' HTML export (if so requested). Links are provided relatively
#' assuming that the output file is two levels below the output
#' directory, e.g., in Output/reports/myanalysis/plate-qc-df.html.
#' @export 
#' @author Henning Redestig
createPlateQcDf <- function(df, html=TRUE) {
  mymax <- function(x) {
    if(all(is.na(x))) return(NA)
    return(max(x, na.rm=TRUE))
  }
  dg <- ddply(df, c("timepoint", "qc_picture"), function(dd) {
    removed <- "None"
    if(any(dd$removed)) 
      removed <- with(subset(dd, removed),
                      paste(ROW, RANGE, collapse=";", sep=":"))
    data.frame(image=dd$image[1],
               rotation=dd$rotation[1],
               eccentricity=mymax(c(dd$deltax[1],dd$deltay[1])),
               ambiguous_boxes=sum(dd$ambig_box),
               max_features=mymax(dd$nfeats),
               removed=removed)
  })
  if(html) {
    dg$image <- file.path("../../..", dg$image)
    dg$image <-
      with(dg, paste('<a href=\"', image, '\">',
                     basename(as.character(image)), "</a>", sep=""))
    dg$qc_picture <- file.path("../../..", dg$qc_picture)
    dg$qc_picture <-
      with(dg, paste('<a href=\"', qc_picture, '\">',
                     basename(as.character(qc_picture)), "</a>", sep=""))
  }
  dg[with(dg, order(-ambiguous_boxes, -max_features, -rotation, -eccentricity)),]
}

#' The date an image was taken
#'
#' Read create date from the exif tags using 
#' @param file an image file (jpg)
#' @return the date
#' @export
#' @references
#' Uses Mayank Lahiri's exif parsing library \url{https://github.com/mayanklahiri/easyexif}
#' @author Henning Redestig
dateTaken <- function(file) {
  exif_date <- .Call("phenotyping_date_original", file, PACKAGE=PKG)
  if(is.null(exif_date))
    return(NA)
  as.POSIXct(strptime(exif_date,"%Y:%m:%d %H:%M:%S"))
}

#' Read and write the meta-data to a file
#'
#' rosettR saves meta data about an each experiment in a special file
#' in the experiment directory. Use these function to read/write meta data
#' to that file in json format.
#' @param meta the meta data object
#' @param path path to the experiment directory
#' @return nothing, used for side effect
#' @export
#' @name rwMeta
writeMeta <- function(meta, path)
  cat(toJSON(meta), file=file.path(path, "meta"))

#' @name rwMeta
#' @export
readMeta <- function(path)
  fromJSON(file.path(path, "meta"))

#' Write manifest to an experiment
#'
#' Not intended for common use. Write manifest to a given
#' experiment. Sorts the manifest after block and position prior to
#' writing it.
#' @param meta the meta data object of class \code{Oni} or
#' \code{lmetadata}
#' @param manifest the manifest data.frame
#' @param path path to the experiment
#' @return nothing
#' @author Henning Redestig
writeManifest <- function(manifest, path)
  write.table(manifest, file=file.path(path, "manifest.txt"),
              row.names=FALSE, sep=",")

#' Read the manifest from an experiment
#'
#' Each image in an experiment is described in the manifest file. Use
#' this function to retrieve the manifest as a data frame.
#' @param path the path to the manifest file itself or to the
#' experiment where a manifest can be found.
#' @return the manifest \code{data.frame}
#' @export
#' @author Henning Redestig
readManifest <- function(path) {
  mf <- read.table(file.path(path, "manifest.txt"), header=TRUE, sep=",")
  mf$image <- as.character(mf$image)
  mf
}

#' Expand experiment design
#'
#' Given the parameters defined in the metadata, expand the balanced
#' design matrix with each germplasm and each treatment represented
#' in each repeat of the experiment. 
#' @param meta the meta data object that defines the experiment design; the
#' number of repeats, timepoints, germplasms and treatments
#' @param platename canonical file name of the plates to use with
#' \code{\link{sprintf}}.
#' @param plateoffset the number that the plate numbering should be
#' offset to. The number of the first plate is 1 + plateoffset.
#' @param plate_order An pre-defined order of the plate after they
#' have been randomized. The plate order \emph{must} sort to correct
#' order when using \code{\link{rank}} but can be either numerical or
#' character.
#' @param ... not used
#' @return a data frame representing the suggested design of the
#' experiment
#' @export
#' @examples
#' data(exampleMetadata)
#' expandManifest(exampleMetadata)
#' @author Henning Redestig
expandManifest <- function(meta,
                           platename="plate%03d.jpg", plateoffset=0,
                           plate_order, ...) {
  if(length(meta$germplasm) == 0)
    stop("no germplasms defined")
  if(length(meta$treatments) == 0)
    stop("no treatments defined")
  if(length(meta$timepoints) == 0)
    stop("no timepoints defined")
  if(meta$nrepeats == 0) {
    warning("number of repeats not defined, assuming 1 repeat")
    meta$nrepeats <- 1
  }
  ntre <- length(meta$treatments)
  nger <- length(meta$germplasm)
  nrep <- meta$nrepeats
  germplasm_region <- unique(meta$griddf$germplasm_region)
  nreg <- length(germplasm_region)
  chun <- (nger / nreg) * ntre
  
  df <-
    data.frame(plate=rep(sprintf(platename,
                   (plateoffset + 1):(plateoffset + (nger / nreg) * ntre * nrep)),
                   each=nreg),
               chunk=rep(1:((nger / nreg) * ntre), each=nrep * nreg),
               treatment=rep(meta$treatments, each=nger  * nrep),
               GERMPLASM=as.vector(apply(matrix(meta$germplasm, nrow=nreg),
                 2, rep, nrep)),
               germplasm_region=germplasm_region,
               BLOCK=as.vector(replicate(chun, rep(sample(nrep), each=nreg))),
               stringsAsFactors=FALSE)
  df$o <- rep(1:(chun * nrep), each=nreg)
  df <- df[order(df$BLOCK),]
  df$position <-
    as.vector(replicate(nrep, rep(sample(chun), each=nreg)))
  df <- df[order(df$o),]
  if(!missing(plate_order)) {
    if(!(length(plate_order) * nreg) == nrow(df))
      stop("plate order does not match experiment dimensions")
    repdf <- data.frame( plt      =rep(rank(plate_order) , each=nreg)
                        ,BLOCK    =rep(gl(nrep, chun)    , each=nreg)
                        ,position =rep(1:chun            , each=nreg)
                        )
    mm <- match(repdf$plt, df$o) + seq(0, (nreg - 1))
    df[mm,c("BLOCK", "position")] <-
      repdf[,c("BLOCK", "position")]
  }
  df <- merge(df, data.frame(timepoint=meta$timepoints))
  if(nrep < length(LETTERS))
    df$label <- with(df, paste(LETTERS[BLOCK], position, sep=""))
  else {
    warning("cannot create labels for more than 26 repeats")
    df$label <- "no_label"
  }
  df$image <-
    file.path(".", sprintf("D%02d", df$timepoint), df$plate)
  df <- df[,-match(c("o", "chunk"), colnames(df))]
  df
}

#' Remove boxes
#'
#' Flag observations in an experiment dataset as removed to inform
#' downstream analysis that the corresponding data should be ignored.
#' @param path path to the experiment (i.e. the top-level folder in
#' which the Output directory is found)
#' @param pattern a character vector listing pattern on the form
#' "[day]/[plateXXX]/[row:column]". All units can be regular
#' expressions, e.g. use .*/plate11./3:.* to remove all boxes in the
#' third row of plate 110 to 119 from all days.
#' @param platere a regular expression that matches the plate
#' name. Can be left missing to use a default of plate[optional
#' spaces][numbers].jpg.
#' @param remove remove if true, un-remove if false.
#' @return nothing, writes in the \code{data.csv} file and saves the
#' \code{phenodata.rda} file for the corresponding experiment
#' @export
#' @author Henning Redestig
#' @examples
#' \dontrun{
#' ## remove the plant in 3:3 on plate120 from all days
#' removeBoxes('.*/plate120/3:3', '..path..')
#' ## list the removed plants
#' load('..path../Output/phenodata.rda')
#' subset(dataset(phenodata), removed)
#' }
removeBoxes <- function(path, pattern, platere="plate\\s*\\d+", remove=TRUE) {
  if(any(!grepl(paste(".+/", platere, "/.+:.+", sep=""), pattern)))
    stop(paste("malformatted plate indication. should be [day]/[plateXXX]/[row]:[col],",
               "for example 18/plate120/2:3"))
  df <- readPhenodata(path)
  day_plate <-
    with(df, paste(timepoint,
                   gsub(".jpg", "", basename(as.character(image))), sep="/"))
  dayregex <- strsplit(strsplit(pattern,"/")[[1]][1], ":")[[1]][1]
  pltregex <- strsplit(strsplit(pattern,"/")[[1]][2], ":")[[1]][1]
  rowregex <- strsplit(strsplit(pattern,"/")[[1]][3], ":")[[1]][1]
  colregex <- strsplit(strsplit(pattern,"/")[[1]][3], ":")[[1]][2]

  selection_to_remove <-
      grepl(dayregex, as.character(df$timepoint)) &
        grepl(pltregex, as.character(df$plate)) &
          grepl(rowregex, as.character(df$ROW)) &
            grepl(colregex, as.character(df$RANGE))

  df$removed[selection_to_remove] <- remove
  writePhenodata(path, df)
}

#' Process a plate experiment
#'
#' After pictures have been placed in their appropriate places in the
#' experiment folder, this function renames the files to
#' plate<number>.jpg and updates the manifest file to indicate the
#' old file names. After this, the experiment is 'sealed' to avoid
#' renaming the files again. Images are sorted after the day they
#' were taken as taken from the EXIF tag in the image file.
#'
#' After processing, images can be analyzed and data analysis reports
#' compiled.
#'
#' This function is intended for interactive use. You probably do not
#' want to use it in programs directly, instead you may want to look
#' at \code{\link{analyzeImage}} and \code{\link{processPlateImages}}.
#' @param path the directory where to find data for this experiment.
#' @param meta the metadata object for this experiment. Read from the
#' experiment's manifest file if missing.
#' @param rename if true, rename image files otherwise assume that
#' images are renamed externally.
#' @param analyze analyze the images or not
#' @param ... passed on to \code{process_plate_images}. Note in
#' particular the argument \code{wincores} that enables the use of
#' more than one CPU.
#' @return The resulting phenodata object if \code{analyze} was true,
#' otherwise \code{NULL}, both invisibly.
#' ## unpack an example experiment to working directory.
#' @examples
#' \dontrun{
#' dummyPath <- createDummyPlateExperiment(242)
#' ## You may want to inspect the folder structure
#' dir(dummyPath)
#' processPlateExperiment(dummyPath, analyze=FALSE)
#' readManifest(dummyExperimentId)
#' ## optionally analyze the pictures as well
#' processPlateExperiment(dummyExperimentId)
#' trashExperiment(as(dummyExperimentId, "Oni"))
#' }
#' @export
#' @author Henning Redestig
processPlateExperiment <- function(path, meta=readMeta(path), rename=TRUE,
                                   analyze=TRUE, ...) {
  outDir <- file.path(path, "Output")
  qcDir <- file.path(outDir, "qc")
  if(file.exists(qcDir))
    stop("previous QC directory already exists ", qcDir)
  mf <- readManifest(path)
  if(file.exists(file.path(outDir, "sealed")) & rename) {
    message("not renaming files since the experiment has already been sealed")
  } else {
    rnmDf <- ddply(mf, "timepoint", function(dd) {
      daydir <- unique(dirname(as.character(dd$image)))
      dd <- dd[with(dd, order(BLOCK, position)),]
      first_region <- unique(dd$germplasm_region)[1]
      expected_pics <-
        basename(as.character(subset(dd, dd$germplasm_region ==
                                       first_region)$image))
      renamingDf(cleanPath(file.path(path, daydir), mustWork=TRUE),
                  expected_pics)
    })
    renameImages(path, rnmDf, dry=!rename, verbose=FALSE)
    rdf <- with(rnmDf, data.frame(image=file.path(".", subdir, newname),
                                   oldname=image, date=date))
    mf <- merge(mf, rdf, by="image", all=TRUE)
    missing_images <- Filter(Negate(file.exists), unique(mf$image))
    if(length(missing_images > 0)) {
      warning(sprintf("%d missing images for %s, dropping them from manifest",
                      length(missing_images), path))
      mf <- mf[file.exists(mf$image),]
    }
    file.create(file.path(path, "Output", "sealed"))
    writeManifest(mf, path)
  }
  res <- NULL
  if(analyze) {
    message("processing images..")
    res <- processPlateImages(path, ...)
    message("making plate gallery..")
    plateGallery(path, "qc")
  }
  invisible(res)
}

#' Rename plate images
#'
#' Images uploaded from a digital camera may have file names that
#' make it difficult to track which image correspond to which
#' plate. This function can be used to order the files after when
#' they were taken (as defined in the EXIF tag) and to rename them to
#' desired file names. This function is destructive, use with care
#' and always do a dry run first to make sure that the new file names
#' are the desired ones.
#'
#' \emph{For typical plate image analysis, you should not use this
#' function}, but rather let \code{\link{processPlateExperiment}}
#' take care of all file renaming. This function is only exported in
#' order to facilitate transition from legacy software.
#' @param path the path to the experiment or path to a directory with
#' images
#' @param newnames the desired new names of the files after they have
#' been ordered by the time they were taken. The newnames should have
#' the correct jpg suffix. If they do not, a warning is given and the
#' jpg suffix is added.
#' @param dry if true, do simulation only, no files are renamed. If
#' false, files will be renamed destructively.
#' @param verbose print messages about the file renaming action to be
#' performed
#' @param subdirs character vector indicating subset of directories to process
#' (see example)
#' @return a data frame that indictes the old and the new name. It is
#' good to save this data frame to be able to revert file renaming
#' @export
#' @examples
#' \dontrun{
#' dummyExperimentId <- createDummyPlateExperiment(242)
#' newnames <- c('plate003.jpg', 'plate002.jpg', 'plate004.jpg', 'plate001.jpg')
#' ## simulate renaming all files
#' renamePlateImages(dummyExperimentId, newnames)
#' ## simulate for a single directory
#' renamePlateImages(dummyExperimentId, newnames, subdirs="D11")
#' ## to actually make changes
#' renaming <- renamePlateImages(dummyExperimentId, newnames, dry=FALSE)
#' list.files(file.path(expdir(dummyExperimentId), 'D11'))
#' ## undo the renaming
#' with(renaming, mapply(file.rename, file.path(expdir(dummyExperimentId), subdir, newname),
#' file.path(expdir(dummyExperimentId), subdir, intermediate)))
#' with(renaming, mapply(file.rename, file.path(expdir(dummyExperimentId), subdir, intermediate),
#' file.path(expdir(dummyExperimentId), subdir, image)))
#' list.files(file.path(expdir(dummyExperimentId), 'D11'))
#' trashExperiment(as(dummyExperimentId, "Oni"))
#' }
#' @author Henning Redestig
renamePlateImages <- function(path, newnames, dry=TRUE, verbose=TRUE,
                                subdirs=NULL) {
  rnm <- renamingDf(path, newnames, subdirs=subdirs)
  renameImages(path, rnm, dry, verbose)
}

renamingDf <- function(path, newnames, pattern="^(EXP.*_*)*D\\d+",
                       subdirs=NULL) {
  dirs <- file.path(path, list.files(path, pattern=pattern))
  if(!is.null(subdirs))
    dirs <- Filter(function(x) basename(x) %in% subdirs, dirs)
  if(length(dirs) == 0) {
    dirs <- path
    path <- dirname(path)
  }
  jpeg_suffix <- grepl("[jJ][pP][eE]*[gG]$", newnames)
  if(!all(jpeg_suffix)) {
    warning("some target newnames do not have the right (jpg) suffix, adjusting this")
    newnames[!jpeg_suffix] <- paste(newnames[!jpeg_suffix], ".jpg", sep="")
  }
  ldply(dirs, renamingDfSingle, newnames)
}

renamingDfSingle <- function(d, newnames) {
  null_df <-
    data.frame(image=character(), date=character(), newname=character(),
               subdir=character(),
               intermediate=character(), stringsAsFactors=FALSE)
  images <- list.files(d, pattern="[jJ][pP][eE]*[gG]$")
  if(length(images) == 0) {
    warning("no jpeg files in advertised directory ", d)
    return(null_df)
  }
  if(!length(images) == length(newnames))
    stop("the number of images in ", d, " (", length(images),
         ") does not match the length of the new names ",
         "(", length(newnames),")")
  date <- vapply(images, function(x) {
    as.character(dateTaken(file.path(d, x)))
  }, character(1), USE.NAMES=FALSE)
  if(any(is.na(date)))
    stop("failed reading date from files ", paste(dQuote(images[is.na(date)])))
  orderTaken <- order(as.POSIXct(date))
  if(any(duplicated(date))) {
    dateIsDuplicated <- date == date[duplicated(date)]
    message(paste(dQuote(images[dateIsDuplicated]), collapse=", "),
            " do not have unique exif timestamps; assuming they sort alphabetically.",
            " Confirm these images were renamed correctly.")
    orderTaken[dateIsDuplicated] <-
      orderTaken[dateIsDuplicated][order(images[dateIsDuplicated])]
  }
  tmpfiles <- replicate(length(images), {
    basename(tempfile(pattern="rnm_", tmpdir=file.path(d), fileext=".jpg"))
  })
  data.frame(image=images[orderTaken], date=date[orderTaken], newname=newnames,
             subdir=basename(d),
             intermediate=tmpfiles, stringsAsFactors=FALSE)
}

renameImages <- function(path, df, dry, verbose) {
  sdf <- df[df$image != df$newname,]
  if(nrow(sdf) > 0) {
    if( dry  &  verbose) { fn <- function(x, y) { message(x, " -> ", y) }
                         } 
    if(!dry  & !verbose) { fn <- file.rename
                         } 
    if(!dry  &  verbose) { fn <- function(x, y) { message(x, " -> ", y);
                                                  file.rename(x, y); }}
    if( dry  & !verbose) { fn <- function(x, y) {}
                         } 
    with(sdf, mapply(fn, file.path(path, subdir, image),
                     file.path(path, subdir, intermediate)))
    with(sdf, mapply(fn, file.path(path, subdir, intermediate),
                     file.path(path, subdir, newname)))
  }
  df
}
