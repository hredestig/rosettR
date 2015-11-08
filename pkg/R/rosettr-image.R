markPlate <- function(mat, radius, pixelsmm) {
  oarea <-
      ((row(mat) - (nrow(mat) / 2))^2 +
       (col(mat) - (ncol(mat) / 2))^2 < (radius * pixelsmm + 1)^2)
  iarea <-
      ((row(mat) - (nrow(mat) / 2))^2 +
       (col(mat) - (ncol(mat) / 2))^2 < (radius * pixelsmm - 1)^2)
  rim <- oarea & !iarea
  mat[rim] <- 1
  mat
}

colorWellFrames <- function(mat, df, boxWidthP) {
  for(i in which(!df$ambig_box)) {
    bi <- box_index(df, i, boxWidthP)
    mat[bi$rows, bi$cols] <-
        frame_box(mat[bi$rows, bi$cols], 1)
  }
  for(i in which(df$ambig_box)) {
    bi <- box_index(df, i, boxWidthP)
    mat[bi$rows, bi$cols] <-
        frame_box(mat[bi$rows, bi$cols], 2)
  }
  mat
}

wellOccupance <- function(largeLabeled, df, feats, boxWidthP) {
  nfeat <- max(largeLabeled)
  res <- matrix(0, nrow=nfeat, ncol=nrow(df))
  rownames(res) <- as.character(1:nfeat)
  for(i in 1:nrow(df)) {
    dd <- df[i,]
    rows <- floor(dd$centerx - (boxWidthP / 2)):floor(dd$centerx + (boxWidthP / 2))
    cols <- floor(dd$centery - (boxWidthP / 2)):floor(dd$centery + (boxWidthP / 2))
    box <- as.vector(imageData(largeLabeled[rows, cols]))
    if(any(box != 0)) {
      tab <- table(box[box > 0])
      res[names(tab),i] <- res[names(tab),i] + tab
    }
  }
  res
}

addWellCenters <- function(griddf, im, pxlsmm, boxWidth, nBoxGrid, plot=FALSE) {
  if(nBoxGrid %% 2 > 0) stop("grid must be symmetric")
  boxWidthP <- pxlsmm * boxWidth
  centerX <- nrow(im) / 2
  centerY <- ncol(im) / 2
  vertical <- sapply(-(nBoxGrid/2):(nBoxGrid/2), function(i)
    centerX + i * boxWidthP)
  horizontal <- sapply(-(nBoxGrid/2):(nBoxGrid/2), function(i)
    centerY + i * boxWidthP)
  if(any(is.na(c(vertical, horizontal))))
      stop("lost coordinates, perhaps wrong pixels / mm")
  if(any(is.na(c(vertical, horizontal))))
      stop("got negative coordinates, perhaps wrong pixels / mm")
  if(plot) {
    plotimage(im)
    abline(v=vertical, col="red")
    abline(h=horizontal, col="blue")
    abline(v=centerX, col="green")
    abline(h=centerY, col="green")
  }
  ddply(griddf, c("ROW", "RANGE"), function(dd) {
    j <- as.integer(dd$ROW)
    i <- dd$RANGE
    dd$centerx <- vertical[i] + (vertical[i + 1] - vertical[i]) / 2
    dd$centery <- horizontal[j] + (horizontal[j + 1] - horizontal[j]) / 2
    if(plot)
        points(dd$centerx, dd$centery, col="red", cex=2, pch=4)
    dd
  })
}

#' Estimate the optimum rotatory correction to an image with objects
#' in straight rows and columns
#'
#' Applies equidistant location-fixed mixed model to the summed
#' horizontal and vertical margin pixel sums to a series of rotated
#' images with the goal to minimize a the common \eqn{\sigma}.
#' @param im an image object
#' @param pxlsmm the number of pixels per millimeter
#' @param boxWidth the number of rows and columns in the expected symmetric grid
#' @param nBoxGrid the width of a box
#' @param angles min and max values for rotation correction
#' @param ... not used
#' @return the suggested rotation correction in degrees
#' @author Henning Redestig
optimAngle <- function(im, pxlsmm, boxWidth, nBoxGrid, angles=c(-3,3), ...) {
  par <- list(0)
  opt <- optim(par, rotationSigmaFun,
               im=im, pxlsmm=pxlsmm, boxWidth=boxWidth, nBoxGrid=nBoxGrid,
               method="Brent", lower=min(angles),
               upper=max(angles))
  opt$par[[1]]
}

rotationSigmaFun <- function(par, im, pxlsmm, boxWidth, nBoxGrid) {
  boxWidthP <- pxlsmm * boxWidth
  rim <- rotate(im, par[[1]])
  rim <- roundImage(rim)
  centerX <- floor(nrow(rim) / 2)
  centerY <- floor(ncol(rim) / 2)
  dataX <- row(rim)[rim > 0]
  dataY <- col(rim)[rim > 0]
  sigmax <- rotationSigma(dataX, centerX, boxWidthP, nBoxGrid)
  sigmay <- rotationSigma(dataY, centerY, boxWidthP, nBoxGrid)
  sigmax + sigmay
}

findCombinedCurve <- function(mu1=NULL, mu2=NULL, sigma=NULL,
                              lambda=NULL, x=NULL) {
  suppressWarnings(res <- lambda * dnorm(x, mean=mu1, sd=sigma) +
                   (1 - lambda) * dnorm(x, mean=mu2, sd=sigma))
  res[is.nan(res)] <- 0
  res
}

diffFun <- function(par, fitX=NULL, fitY=NULL){
  names(par) <- c("mu1", "mu2", "sigma", "lambda")
  myY <-
    findCombinedCurve(mu1=par["mu1"], mu2=par["mu2"],  sigma=par["sigma"],
                        lambda=par["lambda"], x=fitX)
  sum( (myY - fitY)^2 )
}

makeBinary <- function(im, pixelsmm, thresh=NULL, radius=85,
                        fraction=0.25, minNonEmpty=1e-3, max_prop=10,
                        step=0.01, verbose=FALSE) {
  plate_area <-
      (row(im) - (nrow(im) / 2))^2 + (col(im) - (ncol(im) / 2))^2 <
          (pixelsmm * radius)^2
  platevec <- as.vector(imageData(im)[plate_area])
  if(is.null(thresh)) {
    smpl <- sample(platevec, length(platevec) * fraction)
    breaks <- seq(from=min(smpl) - step / 2,
                  to=max(smpl) + 3 * step / 2, by=step)
    max_bin_size <- length(platevec) / length(breaks) / max_prop
    smpl <- unlist(tapply(smpl, cut(smpl, breaks),
                          function(x) sample(x, min(length(x), max_bin_size))))
    fitX <- hist(smpl, breaks, plot=FALSE)$mids
    fitY <- hist(smpl, breaks, plot=FALSE)$density

    par <- c(0.3, 0.95, 0.1, 0.05)
    lower <- c(0.2, 0.9, 0.01, minNonEmpty)
    upper <- c(0.5, 1, 0.3, 1 - minNonEmpty)
    names(lower) <- names(upper) <- names(par) <-
      c("mu1", "mu2", "sigma", "lambda")
    optim_out <- optim(par, diffFun, lower=lower, upper=upper,
                       fitX=fitX, fitY=fitY,
                       method="L-BFGS-B")
    opar <- optim_out$par
    x <- seq(0,1,by=0.001)
    mixd <- cbind(dnorm(x, opar["mu1"], opar["sigma"]) * opar["lambda"],
                  dnorm(x, opar["mu2"], opar["sigma"]) * (1 - opar["lambda"]))
    thresh <-
      x[which(mixd[,1] < mixd[,2])[1]]
    if(verbose)
        message("mix model pars: ", paste(sapply(c(opar, thresh),
                                                 signif, 2), collapse=", "))
  }
  emptiness <- sum(platevec <= thresh) / sum(platevec > thresh)
  if(emptiness < minNonEmpty)
    return(NULL)
  im_mat <- matrix(0, nrow=dim(im)[1], ncol=dim(im)[2])
  im_mat[imageData(im) < thresh] <- 1
  Image(im_mat, dim=dim(im))
}

plateDiffFun <- function(par, mat, radius, threshold) {
  .Call("phenotyping_count_outside", mat,
        as.double(radius),
        as.double(par[[1]]),
        as.double(par[[2]]),
        as.double(threshold), PACKAGE=PKG)
}

#' Find a plate in an image
#'
#' Uses Nelder-Mead optimization to find the position of a circle
#' that captures as many pixels below a given threshold as possible.
#' @param im an image matrix
#' @param radius the radius of the circle
#' @param threshold the threshold pixel value
#' @return a list with optimal deltax and deltay parameters
#' (deviations from origo)
#' @author Henning Redestig
findPlate <- function(im, radius, threshold) {
  opar <-
      optim(list(0, 0), plateDiffFun, mat=im, radius=radius,
            threshold=threshold, method="Nelder-Mead")
  list(deltax=-2*opar$par[[1]], deltay=-2*opar$par[[2]])
}

repositionImage <- function(im, deltax, deltay) {
  if(abs(deltax) > 0 | abs(deltay) > 0) {
    if(length(dim(im)) == 2) {
      if(deltay > 0) {
        im <- im[,(1:ncol(im)) < (ncol(im) - deltay)]
      }
      if(deltay < 0) {
        im <- im[,(1:ncol(im)) > abs(deltay)]
      }
      if(deltax > 0) {
        im <- im[(1:nrow(im)) < (nrow(im) - deltax),]
      }
      if(deltax < 0) {
        im <- im[(1:nrow(im)) > abs(deltax),]
      }
    }
    if(length(dim(im)) == 3) {
      if(deltay > 0) {
        im <- im[,(1:ncol(im)) < (ncol(im) - deltay),]
      }
      if(deltay < 0) {
        im <- im[,(1:ncol(im)) > abs(deltay),]
      }
      if(deltax > 0) {
        im <- im[(1:nrow(im)) < (nrow(im) - deltax),,]
      }
      if(deltax < 0) {
        im <- im[(1:nrow(im)) > abs(deltax),,]
      }
    }
  }
  im
}

makeGrey <- function(im, channels=3) {
  x <- imageData(im)[,,1] * 0
  for(i in channels)
    x <- x + (imageData(im)[,,i])
  x <- x / length(channels)
  Image(x, dim(im)[1:2], "Grayscale")
}

makeSquare <- function(im, cx=NULL) {
  if(!is.null(cx)) if(cx == 0) return(im)
  if(nrow(im) > ncol(im)) {
    if(is.null(cx))
      cx <- (nrow(im) / 2)
    height <- (ncol(im) / 2)
    crop_row <- (cx - height):(cx + height)
    im <- im[crop_row,,]
  }
  if(nrow(im) < ncol(im)) {
    if(is.null(cx))
      cx <- (ncol(im) / 2)
    width <- (nrow(im) / 2)
    crop_col <- (cx - width):(cx + width)
    im <- im[,crop_col,]
  }
  im
}

emptyResult <- function(df, qcpath="", doqc=TRUE) {
  df$centerx <- df$centery <- df$rotation <- df$deltax_mm <-
      df$deltay_mm <- df$ambig_box <- df$total_area_pixels <-
          df$AREA <- df$nfeats <- df$qc_picture <- NA
  if(doqc) {
    file.copy(system.file("examples/empty-qc.jpg", package=PKG), qcpath)
    splitPath <- strsplit(qcpath, .Platform$file.sep)[[1]]
    df$qc_picture <-
      paste(c(".", splitPath[(length(splitPath) - 3):length(splitPath)]),
            collapse=.Platform$file.sep)
  }
  df
}

#' Process an image
#'
#' Does all steps from checking plate location, rotating the picture
#' to ensure grid is horizontal/vertical, gridding, calculating plant
#' areas and generating quality control image.
#' @param file the path of the image to analyze
#' @param griddf data frame that describes the grid. Has three
#' columns, the \code{RANGE} which lists all columns by an integer,
#' \code{ROW} that lists all rows by a factor and \code{box_num} that
#' indicates the number of the box (running integer from 1 to number
#' of boxes). One row per box.
#' @param pixelsmm number of pixels that correspond to 1 mm.
#' @param boxWidth the width (and height) of a single box in the grid in mm.
#' @param nBoxGrid the number of boxes in the grid.
#' @param plateRadius the radius of the plate in millimeters.
#' @param thresh the threshold to appply for deciding what is plant
#' (close zero) or background (close to 1)
#' @param deltax a shift of the places location in millimeters to
#' apply in to the right if negative or left if positive
#' @param deltay a shift of the plates location in millimeters to
#' apply in to the upwards (if negative number) or downwards (if
#' positive number)
#' @param savedir directory to save pictures for quality control purpose
#' @param hires the resolution to use when calculating plant
#' areas. The higher the more precise but also slower and more memory
#' intensive.
#' @param lowres the resolution in pixels along the width the
#' image at which to detect rotation and plate eccentricity. The
#' higher, the more exact but also severly increases computation
#' time. Does not affect area calculations as this is always done at
#' maximum resolution.
#' @param homeratio the ratio of area of the box in which a feature
#' occupies the largest area (home box) and the occupance of the box
#' the feature occupies second most area (neighbor box). With 1/3,
#' the home box must account for 3 times as much area as the neighbor
#' box or the feature will be marked ambiguous and its home will be
#' subjected to grid-segmentation.
#' @param channels the channels to use when converting to grey scale
#' (R, G, B)
#' @param rimThreshold a step of sorting features into different
#' boxes is to identify features that are too far from any box center
#' that may be part of the plate edges or detections of the physical
#' plate itself. Features that are farther away from any box center
#' than half the box size times this constant are considered bad and
#' removed.
#' @param minArea minimum area of a plant. A box with a total
#' feature area less than this number will have area set to missing
#' value.
#' @param rotation the rotation of the plate in degrees
#' @param checkrotation check that plate is placed with grid
#' perpendicular to edges and correct the image if it the rotation
#' angle is more than a given small angle (defined by
#' \code{\link{optimAngle}})
#' @param checklocation check that the plate is in the center of the
#' image and adjust the image otherwise. Ignored if either deltax or
#' deltay are set to a non-zero value.
#' @param doqc generate quality control picture or not
#' @param verbose print some messages about progress
#' @param overwriteQc overwrite any already exisiting QC
#' images.
#' @param squareWidth The desired width of the image after initial
#' cropping. Prior to any further analysis, the image is made square
#' by cropping the left and right margins assuming that the max
#' horizontal excentricity of the plate is smaller than top / bottom
#' margin. Set this argument to zero for very excentric plates in
#' which case cropping is skipped.
#' @param minNonEmpty the minimum fraction of pixels above the threshold for a
#' plate to be considered empty
#' @param ... passed on to \code{\link{findPlate}} and
#' \code{\link{optimAngle}}.
#' @return a data frame with information about plant areas, locations
#' @export
#' @examples
#' \dontrun{
#' file <- system.file('examples/plate', 'plate.jpg', package='BCS.Phenotyping')
#' reps <- cbind(as.vector(sapply(1:6, rep, 6)), 1:6)
#' griddf <- subset(data.frame(RANGE=reps[,2], ROW=LETTERS[reps[,1]]),
#' !(ROW %in% c('A','F') & RANGE %in% c(1,6)))
#' griddf$box_num <- 1:32
#' griddf$too_small <- FALSE
#' df <- analyzeImage(file, griddf, 3.7454, 20, 6, 75, savedir=tempdir(),
#' verbose=TRUE)
#' library(EBImage)
#' display(readImage(df$qc_picture[1]), method='raster')
#' }
#' @author Henning Redestig
analyzeImage <- function(file, griddf, pixelsmm, boxWidth, nBoxGrid, plateRadius,
                         thresh=NULL,
                         deltax=0, deltay=0,
                         savedir=".",
                         hires=1500, lowres=500, homeratio=1/2,
                         channels=3, rimThreshold=0.9, minArea=1.5, rotation=0,
                         checkrotation=TRUE, checklocation=TRUE, doqc=TRUE,
                         verbose=FALSE, overwriteQc=FALSE, squareWidth=NULL,
                         minNonEmpty=1e-3, ...) {
  qcpath <- file.path(savedir, paste("qc_", basename(file), sep=""))
  if(file.exists(qcpath) && !overwriteQc)
    qcpath <-
        tempfile(gsub(".jpg", "", basename(qcpath)),
                 tmpdir=dirname(qcpath), fileext=".jpg")
  if(verbose) message("reading / cropping / scaling picture")
  largeRgb <- readImage(file)
  if(dim(largeRgb)[1] > hires) {
    pixelsmm <- pixelsmm * (hires / nrow(largeRgb))
    largeRgb <- resize(largeRgb, hires)
  }
  largeRgb <- makeSquare(largeRgb, cx=squareWidth)
  largeGrey <- makeGrey(largeRgb, channels=channels)
  smallGrey <- resize(largeGrey, lowres, filter="bilinear")
  scaling <- lowres / nrow(largeRgb)
  boxWidthP <- pixelsmm * boxWidth
  pixelsmmSmall <- pixelsmm * scaling
  if(verbose) message("make binary")
  largeBin <- makeBinary(largeGrey, pixelsmm, thresh, verbose=verbose,
                           minNonEmpty=minNonEmpty)
  if(is.null(largeBin))
    return(emptyResult(griddf, qcpath, doqc))
  smallBin <- resize(largeBin, lowres, filter="bilinear")

  location <- list(deltax=deltax, deltay=deltay)
  if(checklocation & any(location$deltax != 0, location$deltay != 0))
    checklocation <- FALSE
  if(checklocation)
    location <- findPlate(smallGrey, plateRadius * pixelsmmSmall,
                           quantile(smallGrey, 0.15))
  if(verbose) message("plate location: ", location[[1]], ", ", location[[2]])
  largeGrey <-
      repositionImage(largeGrey, deltax * pixelsmm, deltay * pixelsmm)
  smallBin <- repositionImage(smallBin,
                                location$deltax, location$deltay)
  largeBin <- repositionImage(largeBin, location$deltax /
                                scaling, location$deltay / scaling)
  largeRgb <- repositionImage(largeRgb, location$deltax /
                                scaling, location$deltay / scaling)

  if(checkrotation & any(rotation != 0))
    checkrotation <- FALSE
  if(checkrotation)
    rotation <- optimAngle(smallBin, pixelsmmSmall, boxWidth, nBoxGrid, ...)
  if(verbose) message("plate rotation: ", rotation)
  largeRgb <- rotate(largeRgb, rotation)
  largeBin <- roundImage(rotate(largeBin, rotation))

  df <- addWellCenters(griddf, largeBin, pixelsmm, boxWidth, nBoxGrid)
  df$rotation <- rotation
  df$deltax_mm <- location$deltax / pixelsmm
  df$deltay_mm <- location$deltay / pixelsmm
  df$ambig_box <- FALSE
  ntries <- 0
  if(verbose) message("sorting features to boxes")
  repeat {
    ntries <- ntries + 1
    largeLabeled <- bwlabel(largeBin)
    feats <- computeFeatures.shape(largeLabeled)
    if(ntries > nrow(df)) {
      stop("(bug) failed to resolve all ambiguities")
    }
    feats <- as.data.frame(feats)
    boxocc <- wellOccupance(largeLabeled, df, feats, boxWidthP)
    feats$home_box <- apply(boxocc, 1, which.max)
    overlapping <- apply(boxocc, 1, function(oc) {
      if(sum(oc) == 0) return(NA) #completely outside the grid
      o <- order(oc, decreasing=TRUE)
      oc[o[2]] / oc[o[1]] > homeratio
    })
    mat <- imageData(largeLabeled)
    points <- as.matrix(df[,c("centerx", "centery")])
    closestCenters <-
      .Call("phenotyping_closest_point", mat, points, nrow(feats), PACKAGE=PKG)
    far_away <- apply(closestCenters, 1, min) > (boxWidthP / 2) * 0.9
    overlapping[far_away] <- NA
    if(any(na.omit(overlapping))) {
      problem_boxes <-
          unique(apply(boxocc[which(overlapping),,drop=FALSE], 1,
                       function(x) order(x, decreasing=TRUE)[1:2]))
      df$ambig_box[problem_boxes] <- TRUE
      for(i in problem_boxes) {
        bi <- wellIndex(df, i, boxWidthP)
        box <- largeBin[bi$rows, bi$cols]
        largeBin[bi$rows, bi$cols] <- frameWell(box, 0)
      }
    }
    if(any(is.na(overlapping))) {
      index_mat <- imageData(largeLabeled) %in% which(is.na(overlapping))
      imageData(largeBin)[index_mat] <- 0
    }
    if(!any(is.na(overlapping), overlapping))
        break
  }
  df <- ddply(df, c("ROW", "RANGE"), function(dd) {
    totpix <- 0
    inThisBox <- feats$home_box == dd$box_num
    totpix <- sum(feats$s.area[inThisBox])
    dd$total_area_pixels <- totpix
    dd$AREA <- dd$total_area_pixels / pixelsmm^2
    dd$nfeats <- sum(inThisBox)
    if(sum(inThisBox) == 0)
      dd$AREA <- NA
    if(!is.na(dd$AREA)) {
      if(dd$AREA < minArea)
        dd$too_small <- TRUE
    }
    dd
  })

  if(doqc) {
    if(verbose) message("doing qc")
    mat <- imageData(largeLabeled)
    nfeats <- max(largeLabeled)
    for(i in 1:nrow(feats)) {
      box_idx <- match(feats$home_box[i], df$box_num)
      if(is.na(df$AREA[box_idx]) | df$too_small[box_idx])
        mat[mat == i] <- nfeats + 1
      else
        mat[mat == i] <- nfeats + 2 + feats$home_box[i]
    }
    mat <- mat - nfeats
    mat[mat < 0] <- 0
    mat <- colorWellFrames(mat, df, boxWidthP)
    mat <- markPlate(mat, plateRadius, pixelsmm)
    cols <- c("black", "black", "red", sample(rainbow(nrow(df))))
    qcpic <- Image(matrix(cols[1 + imageData(mat)], nrow=nrow(mat),
                          ncol=ncol(mat)),
                   colormode="color")
    largeRgb[mat != 0] <- 0
    qcpic <- qcpic + largeRgb
    qcpic <- resize(qcpic, lowres)
    scaling <- lowres / nrow(largeRgb)
    pixelsmmSmall <- pixelsmm * scaling
    xy <- as.matrix(ddply(df, c("ROW", "RANGE"), function(dd) {
      c((dd$centerx * scaling) - (boxWidth * pixelsmmSmall) / 3,
        (dd$centery * scaling) - (boxWidth * pixelsmmSmall) / 3)
    })[,c("V1", "V2")])
    qcpic <- drawText(qcpic, x=xy[,1], y=xy[,2],
                      labels=paste(df$ROW, df$RANGE, sep=":"),
                      col="black")
    param_lab <- sprintf("dx=%.2f\ndy=%.2f\nr=%.2f",
                         location$deltax, location$deltay, rotation)
    param_xy <- c(min(xy[,1]) - .5 * (boxWidth * pixelsmmSmall),
                  min(xy[,2]) - .5 * (boxWidth * pixelsmmSmall))
    qcpic <- drawText(qcpic, x=param_xy[1], y=param_xy[2],
                      labels=param_lab)
    if(!file.exists(savedir))
      dir.create(savedir, recursive=TRUE)
    writeImage(qcpic, qcpath)
    splitPath <- strsplit(qcpath, .Platform$file.sep)[[1]]
    df$qc_picture <-
      paste(c(".", splitPath[(length(splitPath) - 3):length(splitPath)]),
            collapse=.Platform$file.sep)
    
  }
  rownames(df) <- paste(df$ROW, df$RANGE, sep='')
  df
}

#' Calibrate the scale of an image
#'
#' Used to figure out how many pixels there are on 1 mm by allowing
#' the user to click on an image. Updates the meta data with the
#' calculated pixels per millimeter
#' @param path path to the experiment to calibrate the scale for
#' @return the number of pixels per mm
#' @export
calibrateScale <- function(path) {
  meta <- readMeta(path)
  mf <- readManifest(path)
  im <- readImage(file.path(path, mf$image[1]))
  done <- FALSE
  while(!done) {
    plotimage(im)
    answer <- readline("choose the distance to indicate [default=80mm]: ")
    selectedDistance <- ifelse(answer == "", 80, as.integer(answer))
    message("left-click on two points separated ", selectedDistance, "mm")
    pp <- locator(n=2, type="l")
    npixels <- sqrt(diff(pp$x)^2 + diff(pp$y)^2)
    answer <- readline("try again? [default=n]: ")
    done <- ifelse(answer %in% c("", "No", "NO", "n", "nO"), TRUE, FALSE)
  }
  meta$pixelsmm <- npixels / selectedDistance
  writeMeta(meta, path)
  message(npixels / selectedDistance, " pixels per mm, meta has been updated")
  invisible(npixels / selectedDistance)
}

#' Process a set of plate images
#'
#' Apply \code{analyzeImage} to all images in an experiment.
#'
#' Parallelization is done using the the plyr package plugin to the
#' \code{foreach} package. \code{doParallel} is used on Linux but the
#' multi-core functionality is (deliberately) defunct on
#' Windows.. Missing required arguments for \code{analyzeImage} are
#' taken from the loaded configuration.
#' @param path the path to the experiment
#' @param mf the manifest data frame for the experiment that must
#' have a column named 'image' that specifies the path to each image
#' as well as 'timepoint' that specifies the timepoint at which the
#' image was taken.
#' @param verbose display progressbar or not. Disabled for parallel
#' computation as it does not work reliably.
#' @param save should the created phenodata object and a csv file
#' with the final dataset be placed in the output directory of the
#' experiment or not
#' @param ... passed on to \code{\link{analyzeImage}} and
#' \code{\link{doProcessPlateImages}}
#' @return invisibly, the \code{phenodata} object with the results.
#' @export
#' @author Henning Redestig
processPlateImages <- function(path, mf=readManifest(path), meta=readMeta(path),
                               verbose=FALSE, save=TRUE, ...) { 
  phenodata <-
    doProcessPlateImages(path=path, mf=mf, meta=meta, verbose=verbose, ...)

  phenodata <- postProcessPlatePhenodata(phenodata, ...)
  if(save)
    writePhenodata(path, phenodata)
  invisible(phenodata)
}

#' The function that actually process the plate images
#'
#' Not intended for user-level usage. 
#' @param path 
#' @param mf the manifest - a data frame describing all images to be processed
#' @param meta a list with parameters to use for \code{\link{analyzeImage}}, see
#' \code{\link{metaTemplate}}
#' @param verbose if true, display progress bar
#' @param .progress passed to \code{\link{ddply}}
#' @param .parallel passed to \code{\link{ddply}}
#' @param ... passed on to \code{analyzeImage}. Named arguments have precedence
#' over arguments in the \code{meta} list.
#' @return the resulting data frame
#' @export
#' @author Henning Redestig
doProcessPlateImages <- function(path, mf=readManifest(path),
                                 meta=readMeta(path), verbose=FALSE, 
                                 .progress="text", .parallel=FALSE, ...) {
  meta <- updateMeta(..., meta=meta)
  ddply(mf, "image", function(dd) {
    file <- file.path(path, unique(as.character(dd$image)))
    tpt <- unique(as.integer((as.character(dd$timepoint))))
    if(length(file) != 1 | length(tpt) != 1)
      stop("multiple files or time points for the same image")
    meta$file <- file
    meta$savedir <- file.path(path, "Output", "qc", sprintf("D%02d", tpt))
    pp <- do.call("analyzeImage", meta)
    merge(pp, dd)
  }, .progress=.progress, .parallel=.parallel)
}

#' Post-process phenodata
#'
#' Apply post-process steps to the phenodata obtained from
#' \code{analyzeImage}. Currently this consists of marking all
#' measurements for plants that at at least one day are not larger
#' than the specified minimum area as 'removed'.
#' @param phenodata a phenodata object to post process
#' @param minArea minimum area of a plant. A box with a total
#' feature area less than this number will have area set to missing
#' value.
#' @param ... not used
#' @return phenodata object
#' @author Henning Redestig
postProcessPlatePhenodata <- function(phenodata, minArea=1.5, ...) {
  ddply(phenodata, c("plate", "ROW", "RANGE"), function(dd) {
    tmp_area <- dd$AREA
    tmp_area[is.na(tmp_area)] <- Inf
    dd$too_small <- dd$too_small | tmp_area < minArea
    if(any(dd$too_small))
      dd$removed <- TRUE
    dd
  })
}

#' Re-process a (sub)set of previously processed plate images
#'
#' Apply \code{\link{analyzeImage}} to a selected set of images in
#' an experiment and combine the results with previously computed
#' results.
#' @param path the path to the experiment to reprocess
#' @param mf a manifest data frame for the selected images.
#' @param meta a list with parameters to configure the image analysis (passed as
#' arguments to \code{\link{analyzeImage}}
#' @param verbose display progressbar or not. Disabled for parallel
#' computation as it does not work reliably.
#' @param save should the created phenodata object and csv-file with
#' the complemented dataset be placed in the output directory of the
#' experiment or not. Previous results for selected images are lost.
#' @param ... passed on to \code{\link{analyzeImage}} (precedence over the
#' parameters defined in \code{meta}) and \code{\link{doProcessPlateImages}}
#' @param cores number of cores to use.
#' @seealso \code{\link{processPlateImages}} that has the same arguments and
#' functionality but is used for all images in the experiment.
#' @return invisibly, the complete (all images, not only those that
#' were re-processed) \code{phenodata} object with the results.
#' @examples
#' \dontrun{
#' untar(system.file('examples/plate/EXP242.tar.gz', package='BCS.Phenotyping'), exdir='.')
#' load_config('plate-6x6.top_hh_bottom_HH')
#' ## plate 1 is not in the centre so it fails at cell 3:3 and 3:4
#' pda <- process_plate_experiment('EXP242_HR', checkrotation=FALSE, verbose=TRUE)
#' qcpic <- subset(dataset(pda), timepoint == 18 & plate == 'plate001.jpg')$qc[1]
#' display(readImage(qcpic), method='raster')
#' mf <- read_manifest('EXP242_HR')
#' mf <- subset(mf, timepoint == 18 & plate  == 'plate001.jpg')
#' ## reprocess the plate applying a shift of 5 mm in the horizontal direction
#' pda <- reprocess_plate_images('EXP242_HR', mf, deltax=5, checkrotation=FALSE)
#' qcpic <- subset(dataset(pda), timepoint == 18 & plate == 'plate001.jpg')$qc[1]
#' display(readImage(qcpic), method='raster')
#' }
#' @export
#' @author Henning Redestig
reprocessPlateImages <- function(path, mf, meta=readMeta(path), verbose=FALSE,
                                 save=TRUE, ...) {
  meta <- updateMeta(..., meta)
  df <- readPhenodata(path)
  newphenodata <-
    doProcessPlateImages(path, mf, meta, ...)
  df <- rbind(df[!(df$image %in% newphenodata$image), , drop=FALSE],
              newphenodata)
  if(save)
    writePhenodata(path, df)
  invisible(phenodata)
}

rotationSigma <- function(dat, center, boxWidthP, nBoxGrid, sigma=10, eps=1e-2) {
  mu0 <- center - boxWidthP / 2 - ((nBoxGrid / 2) - 1) * boxWidthP
  mu <- sapply(1:nBoxGrid, function(i) mu0 + (i - 1) * boxWidthP)
  lambda <- rep(1 / nBoxGrid, nBoxGrid)
  n <- length(dat)
  change <- 1
  while(change > eps) {
    psigma <- sigma
    sigmavec <- rep(sigma, nBoxGrid)
    zc <- .Call("phenotyping_normalpostp", dat, mu, sigmavec, lambda,
                PACKAGE=PKG)
    sigma <-
        sqrt(sum(sapply(1:nBoxGrid,
                        function(i) sum(zc$scaled_resp[,i] *
                                          (dat - mu[i])^2))) / n)
    change <- (psigma - sigma)^2
  }
  sigma
}

rotationSigmaFun <- function(par, im, pxlsmm, boxWidth, nBoxGrid) {
  boxWidthP <- pxlsmm * boxWidth
  rim <- rotate(im, par[[1]])
  rim <- roundImage(rim)
  centerX <- floor(nrow(rim) / 2)
  centerY <- floor(ncol(rim) / 2)
  dataX <- row(rim)[rim > 0]
  dataY <- col(rim)[rim > 0]
  sigmax <- rotationSigma(dataX, centerX, boxWidthP, nBoxGrid)
  sigmay <- rotationSigma(dataY, centerY, boxWidthP, nBoxGrid)
  sigmax + sigmay
}

#' Round the values in a binary single-frame image to 0 or 1.
#' @param im an Image object
#' @return the rounded Image
#' @author Henning Redestig
roundImage <- function(im) {
  if(length(dim(imageData(im))) != 2)
      stop("can only round single frame / channel images")
  mat <- imageData(im)
  imageData(im) <- .Call("phenotyping_round_binary", mat, PACKAGE=PKG)
  im
}

#' Simple threshold image
#'
#' Make an image binary by setting all values above a threshold to 1 and all
#' values below to zero
#' @param im an image
#' @param thresh a threshold
#' @return the binary image
#' @export 
#' @author Henning Redestig
#' library(EBImage)
#' lena <- readImage(system.file("images", "lena-color.png", package="EBImage"))
#' threshImage(lena, 0.9)
threshImage <- function(im, thresh) {
  im[im > thresh] <- 1
  im[im < thresh] <- 0
  im
}

#' Draw text on images
#'
#' Copied from EBImage version 3.13.1 since it was removed in version 4
#' @param img An \code{Image} object or an array.
#' @param xy Matrix (or a list of matrices if \code{img} contains
#' multiple frames) of coordinates of labels.
#' @param labels A character vector (or a list of vectors if
#' \code{img} contains multiple frames) containing the labels to be
#' output.
#' @param font A font object, returned by \code{\link{drawFont}}. If
#' missing, a default OS-dependent font will be chosen.
#' @param col A character vector of font colors.
#' @return An \code{Image} object or an array, containing the
#' transformed version of \code{img}.
#' @export
#' @examples
#' library(EBImage)
#' plate <- readImage(system.file("examples", "plate_merged.jpg", package="rosettR"))
#' font <- drawFont(weight=600, size=28)
#' plate <- drawText(plate, xy=c(250, 450), labels="rosettR!", font=font, col="red")
#' display(plate, method="raster")
#' @author Oleg Sklyar and copied by Henning Redestig
drawText <- function(img, x, y, labels, col="black", cex=1) {
  tmpFile <- tempfile()
  on.exit(unlink(tmpFile))
  jpeg(filename=tmpFile, width=ncol(img), height=nrow(img))
  display(img, method="raster")
  text(x, y, labels=labels, col=col, cex=cex)
  dev.off()
  readImage(tmpFile, type="jpeg")
}

#' Frame a well
#'
#' Set the pixels around perimeter of a well to a given value
#' @param well a matrix (or similar object) that represents the well
#' @param value the value to set the pixels to
#' @return the framed well
#' @export 
#' @author Henning Redestig
frameWell <- function(well, value) {
  well[1:2,] <- value
  well[,1:2] <- value
  well[(nrow(well) - 1):nrow(well),] <- value
  well[,(nrow(well) - 1):ncol(well)] <- value
  well
}

#' The indices in a well
#'
#' Get the matrix indices for all pixels in a given well
#' @param griddf the data frame that defines the plate grid 
#' @param i the index of the well to return (row in the grid data frame )
#' @param wellWidth the width (in pixels) of the well
#' @return a list of row and column indices
#' @export 
#' @author Henning Redestig
wellIndex <- function(griddf, i, wellWidth) {
  dd <- griddf[i, ]
  list(rows=floor(dd$centerx - (wellWidth / 2)) :
         floor(dd$centerx + (wellWidth / 2)),
       cols=floor(dd$centery - (wellWidth / 2)) :
         floor(dd$centery + (wellWidth / 2)))
}

#' Color the frame of all wells
#'
#' Set the permiters of all wells to 1 or 2 depending on whether they are
#' indicated as being ambiguous 
#' @param mat a matrix representing the image
#' @param griddf a data frame that defines the grid 
#' @param wellWidth 
#' @return the image with frames indicated
#' @export
#' @author Henning Redestig
colorWellFrames <- function(mat, griddf, wellWidth) {
  for(i in which(!griddf$ambig_box)) {
    bi <- wellIndex(griddf, i, wellWidth)
    mat[bi$rows, bi$cols] <-
        frameWell(mat[bi$rows, bi$cols], 1)
  }
  for(i in which(griddf$ambig_box)) {
    bi <- wellIndex(griddf, i, wellWidth)
    mat[bi$rows, bi$cols] <-
        frameWell(mat[bi$rows, bi$cols], 2)
  }
  mat
}
