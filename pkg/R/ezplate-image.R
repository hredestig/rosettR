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

colorWellFrames <- function(mat, df, dp) {
  for(i in which(!df$ambig_box)) {
    bi <- box_index(df, i, dp)
    mat[bi$rows, bi$cols] <-
        frame_box(mat[bi$rows, bi$cols], 1)
  }
  for(i in which(df$ambig_box)) {
    bi <- box_index(df, i, dp)
    mat[bi$rows, bi$cols] <-
        frame_box(mat[bi$rows, bi$cols], 2)
  }
  mat
}

wellOccupance <- function(large_labeled, df, feats, dp) {
  nfeat <- max(large_labeled)
  res <- matrix(0, nrow=nfeat, ncol=nrow(df))
  rownames(res) <- as.character(1:nfeat)
  for(i in 1:nrow(df)) {
    dd <- df[i,]
    rows <- floor(dd$centerx - (dp / 2)):floor(dd$centerx + (dp / 2))
    cols <- floor(dd$centery - (dp / 2)):floor(dd$centery + (dp / 2))
    box <- as.vector(imageData(large_labeled[rows, cols]))
    if(any(box != 0)) {
      tab <- table(box[box > 0])
      res[names(tab),i] <- res[names(tab),i] + tab
    }
  }
  res
}

addWellCenters <- function(griddf, im, pxlsmm, d, r, plot=FALSE) {
  if(r %% 2 > 0) stop("grid must be symmetric")
  dp <- pxlsmm * d
  center_x <- nrow(im) / 2
  center_y <- ncol(im) / 2
  vertical <- sapply(-(r/2):(r/2), function(i) center_x + i * dp)
  horizontal <- sapply(-(r/2):(r/2), function(i) center_y + i * dp)
  if(any(is.na(c(vertical, horizontal))))
      stop("lost coordinates, perhaps wrong pixels / mm")
  if(any(is.na(c(vertical, horizontal))))
      stop("got negative coordinates, perhaps wrong pixels / mm")
  if(plot) {
    plotimage(im)
    abline(v=vertical, col="red")
    abline(h=horizontal, col="blue")
    abline(v=center_x, col="green")
    abline(h=center_y, col="green")
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
#' @param d the number of rows and columns in the expected symmetric grid
#' @param r the width of a box
#' @param angles min and max values for rotation correction
#' @return the suggested rotation correction in degrees
#' @author Henning Redestig
optimAngle <- function(im, pxlsmm, d, r, angles=c(-3,3)) {
  par <- list(0)
  opt <- optim(par, rotationSigmaFun,
               im=im, pxlsmm=pxlsmm, d=d, r=r,
               method='Brent', lower=min(angles),
               upper=max(angles))
  opt$par[[1]]
}

rotationSigmaFun <- function(par, im, pxlsmm, d, r) {
  dp <- pxlsmm * d
  rim <- rotate(im, par[[1]])
  rim <- roundImage(rim)
  center_x <- floor(nrow(rim) / 2)
  center_y <- floor(ncol(rim) / 2)
  data_x <- row(rim)[rim > 0]
  data_y <- col(rim)[rim > 0]
  sigmax <- rotationSigma(data_x, center_x, dp, r)
  sigmay <- rotationSigma(data_y, center_y, dp, r)
  sigmax + sigmay
}

findCombinedCurve <- function(mu1=NULL, mu2=NULL, sigma=NULL, lambda=NULL, x=NULL) {
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
                        fraction=0.25, min_non_empty=1e-3, max_prop=10,
                        step=0.01, verbose=FALSE) {
  plate_area <-
    (row(im) - (nrow(im) / 2))^2 + (col(im) - (ncol(im) / 2))^2 < (pixelsmm * radius)^2
  platevec <- as.vector(imageData(im)[plate_area])
  if(is.null(thresh)) {
    smpl <- sample(platevec, length(platevec) * fraction)
    breaks <- seq(from = min(smpl) - step / 2, to = max(smpl) + 3 * step / 2, by=step)
    max_bin_size <- length(platevec) / length(breaks) / max_prop
    smpl <- unlist(tapply(smpl, cut(smpl, breaks),
                          function(x) sample(x, min(length(x), max_bin_size))))
    fitX <- hist(smpl, breaks, plot=FALSE)$mids
    fitY <- hist(smpl, breaks, plot=FALSE)$density

    par <- c(0.3, 0.95, 0.1, 0.05)
    lower <- c(0.2, 0.9, 0.01, min_non_empty)
    upper <- c(0.5, 1, 0.3, 1 - min_non_empty)
    names(lower) <- names(upper) <- names(par) <-
      c("mu1", "mu2", "sigma", "lambda")
    optim_out <- optim(par, diffFun, lower=lower, upper=upper, fitX=fitX, fitY=fitY,
                       method="L-BFGS-B")
    opar <- optim_out$par
    x <- seq(0,1,by=0.001)
    mixd <- cbind(dnorm(x, opar["mu1"], opar["sigma"]) * opar["lambda"],
                  dnorm(x, opar["mu2"], opar["sigma"]) * (1 - opar["lambda"]))
    thresh <-
      x[which(mixd[,1] < mixd[,2])[1]]
    if(verbose)
      message("mix model pars: ", paste(sapply(c(opar, thresh), signif, 2), collapse=", "))
  }
  emptiness <- sum(platevec <= thresh) / sum(platevec > thresh)
  if(emptiness < min_non_empty)
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
    df$qc_picture <- qcpath
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
#' @param d the width (and height) of a single box in the grid in mm.
#' @param r the number of boxes in the grid.
#' @param plate_radius the radius of the plate in millimeters.
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
#' @param rim_threshold a step of sorting features into different
#' boxes is to identify features that are too far from any box center
#' that may be part of the plate edges or detections of the physical
#' plate itself. Features that are farther away from any box center
#' than half the box size times this constant are considered bad and
#' removed.
#' @param min_area minimum area of a plant. A box with a total
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
#' @param overwrite_qc overwrite any already exisiting QC
#' images.
#' @param square_width The desired width of the image after initial
#' cropping. Prior to any further analysis, the image is made square
#' by cropping the left and right margins assuming that the max
#' horizontal excentricity of the plate is smaller than top / bottom
#' margin. Set this argument to zero for very excentric plates in
#' which case cropping is skipped.
#' @param min_non_empty the minimum fraction of pixels above the threshold for a
#' plate to be considered empty
#' @param ... passed on to \code{\link{findPlate}} and
#' \code{\link{optimAngle}}.
#' @return a data frame with information about plant areas, locations
#' @export
#' @examples
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
#' @author Henning Redestig
analyzeImage <- function(file, griddf, pixelsmm, d, r, plate_radius, thresh=NULL,
                         deltax=0, deltay=0,
                         savedir=".",
                         hires=1500, lowres=500, homeratio=1/2,
                         channels=3, rim_threshold=0.9, min_area=1.5, rotation=0,
                         checkrotation=TRUE, checklocation=TRUE, doqc=TRUE,
                         verbose=FALSE, overwrite_qc=FALSE, square_width=NULL,
                         min_non_empty=1e-3, ...) {

  qcpath <- file.path(savedir, paste("qc_", basename(file), sep=""))
  if(file.exists(qcpath) && !overwrite_qc)
    qcpath <-
      tempfile(gsub(".jpg", "", basename(qcpath)), tmpdir=dirname(qcpath), fileext=".jpg")

  if(verbose) message("reading / cropping / scaling picture")
  large_rgb <- readImage(file)
  if(dim(large_rgb)[1] > hires) {
    pixelsmm <- pixelsmm * (hires / nrow(large_rgb))
    large_rgb <- resize(large_rgb, hires)
  }
  large_rgb <- makeSquare(large_rgb, cx=square_width)
  large_grey <- makeGrey(large_rgb, channels=channels)
  small_grey <- resize(large_grey, lowres, filter="bilinear")
  scaling <- lowres / nrow(large_rgb)
  dp <- pixelsmm * d
  pixelsmm_small <- pixelsmm * scaling
  if(verbose) message("make binary")
  large_bin <- makeBinary(large_grey, pixelsmm, thresh, verbose=verbose,
                           min_non_empty=min_non_empty)
  if(is.null(large_bin))
    return(emptyResult(griddf, qcpath, doqc))
  small_bin <- resize(large_bin, lowres, filter="bilinear")

  location <- list(deltax=deltax, deltay=deltay)
  if(checklocation & any(location$deltax != 0, location$deltay != 0))
    checklocation <- FALSE
  if(checklocation)
    location <- findPlate(small_grey, plate_radius * pixelsmm_small,
                           quantile(small_grey, 0.15))
  if(verbose) message("plate location: ", location[[1]], ", ", location[[2]])
  large_grey <-
      repositionImage(large_grey, deltax * pixelsmm, deltay * pixelsmm)
  small_bin <- repositionImage(small_bin,
                                location$deltax, location$deltay)
  large_bin <- repositionImage(large_bin, location$deltax /
                                scaling, location$deltay / scaling)
  large_rgb <- repositionImage(large_rgb, location$deltax /
                                scaling, location$deltay / scaling)

  if(checkrotation & any(rotation != 0))
    checkrotation <- FALSE
  if(checkrotation)
    rotation <- optimAngle(small_bin, pixelsmm_small, d, r, ...)
  if(verbose) message("plate rotation: ", rotation)
  large_rgb <- rotate(large_rgb, rotation)
  large_bin <- roundImage(rotate(large_bin, rotation))

  df <- addWellCenters(griddf, large_bin, pixelsmm, d, r)
  df$rotation <- rotation
  df$deltax_mm <- location$deltax / pixelsmm
  df$deltay_mm <- location$deltay / pixelsmm
  df$ambig_box <- FALSE
  ntries <- 0
  if(verbose) message("sorting features to boxes")
  repeat {
    ntries <- ntries + 1
    large_labeled <- bwlabel(large_bin)
    feats <- computeFeatures.shape(large_labeled)
    if(ntries > nrow(df)) {
      stop("(bug) failed to resolve all ambiguities")
    }
    feats <- as.data.frame(feats)
    boxocc <- wellOccupance(large_labeled, df, feats, dp)
    feats$home_box <- apply(boxocc, 1, which.max)
    overlapping <- apply(boxocc, 1, function(oc) {
      if(sum(oc) == 0) return(NA) #completely outside the grid
      o <- order(oc, decreasing=TRUE)
      oc[o[2]] / oc[o[1]] > homeratio
    })
    mat <- imageData(large_labeled)
    points <- as.matrix(df[,c("centerx", "centery")])
    closest_centers <-
      .Call("phenotyping_closest_point", mat, points, nrow(feats), PACKAGE=PKG)
    far_away <- apply(closest_centers, 1, min) > (dp / 2) * 0.9
    overlapping[far_away] <- NA
    if(any(na.omit(overlapping))) {
      problem_boxes <-
          unique(apply(boxocc[which(overlapping),,drop=FALSE], 1,
                       function(x) order(x, decreasing=TRUE)[1:2]))
      df$ambig_box[problem_boxes] <- TRUE
      for(i in problem_boxes) {
        bi <- wellIndex(df, i, dp)
        box <- large_bin[bi$rows, bi$cols]
        large_bin[bi$rows, bi$cols] <- frameWell(box, 0)
      }
    }
    if(any(is.na(overlapping))) {
      index_mat <- imageData(large_labeled) %in% which(is.na(overlapping))
      imageData(large_bin)[index_mat] <- 0
    }
    if(!any(is.na(overlapping), overlapping))
        break
  }
  df <- ddply(df, c("ROW", "RANGE"), function(dd) {
    totpix <- 0
    in_this_box <- feats$home_box == dd$box_num
    is_bad <- any(feats$ambig_feat[in_this_box])
    totpix <- sum(feats$s.area[in_this_box])
    dd$total_area_pixels <- totpix
    dd$AREA <- dd$total_area_pixels / pixelsmm^2
    dd$nfeats <- sum(in_this_box)
    if(sum(in_this_box) == 0)
      dd$AREA <- NA
    if(!is.na(dd$AREA)) {
      if(dd$AREA < min_area)
        dd$too_small <- TRUE
    }
    dd
  })

  if(doqc) {
    if(verbose) message("doing qc")
    mat <- imageData(large_labeled)
    nfeats <- max(large_labeled)
    for(i in 1:nrow(feats)) {
      box_idx <- match(feats$home_box[i], df$box_num)
      if(is.na(df$AREA[box_idx]) | df$too_small[box_idx])
        mat[mat == i] <- nfeats + 1
      else
        mat[mat == i] <- nfeats + 2 + feats$home_box[i]
    }
    mat <- mat - nfeats
    mat[mat < 0] <- 0
    mat <- colorWellFrames(mat, df, dp)
    mat <- markPlate(mat, plate_radius, pixelsmm)
    cols <- c("black", "black", "red", sample(rainbow(nrow(df))))
    qcpic <- Image(cols[1 + imageData(mat)],
                            dim=dim(mat))
    large_rgb[mat != 0] <- 0
    qcpic <- qcpic + large_rgb
    qcpic <- resize(qcpic, lowres)
    scaling <- lowres / nrow(large_rgb)
    pixelsmm_small <- pixelsmm * scaling
    font_11 <- drawFont(weight=500, size=11)
    xy <- as.matrix(ddply(df, c("ROW", "RANGE"), function(dd) {
      c((dd$centerx * scaling) - (d * pixelsmm_small) / 3,
        (dd$centery * scaling) - (d * pixelsmm_small) / 3)
    })[,c("V1", "V2")])
    qcpic <- drawText(qcpic, xy=xy,
                      labels=paste(df$ROW, df$RANGE, sep=":"),
                      font=font_11, col="black")
    param_lab <- sprintf("dx=%.2f\ndy=%.2f\nr=%.2f",
                         location$deltax, location$deltay, rotation)
    param_xy <- c(min(xy[,1]) - .5 * (d * pixelsmm_small),
                  min(xy[,2]) - .5 * (d * pixelsmm_small))
    qcpic <- drawText(qcpic, xy=param_xy,
                      labels=param_lab,
                      font=font_11, col="black")
    if(!file.exists(savedir))
      dir.create(savedir, recursive=TRUE)
    writeImage(qcpic, qcpath)
    df$qc_picture <- qcpath
  }
  rownames(df) <- paste(df$ROW, df$RANGE, sep='')
  df
}

neighborWell <- function(griddf, i, direction=c("north", "south", "east", "west")) {
  direction <- match.arg(direction)
  rowcol <- switch(direction,
                   north = {
                     c(subset(griddf, griddf$box_num == i)$ROW - 1,
                       subset(griddf, griddf$box_num == i)$RANGE)
                   },
                   south = {
                     c(subset(griddf, griddf$box_num == i)$ROW + 1,
                       subset(griddf, griddf$box_num == i)$RANGE)
                   },
                   east = {
                     c(subset(griddf, griddf$box_num == i)$ROW,
                       subset(griddf, griddf$box_num == i)$RANGE - 1)
                   },
                   south = {
                     c(subset(griddf, griddf$box_num == i)$ROW,
                       subset(griddf, griddf$box_num == i)$RANGE + 1)
                   })
  if(any(griddf$ROW == rowcol[1] & griddf$RANGE == rowcol[2]))
    return(griddf$box_num[griddf$ROW == rowcol[1] & griddf$RANGE == rowcol[2]])
  else
    return(NA)
}

#' Calibrate the scale of an image
#'
#' Used to figure out how many pixels there are on 1 mm by allowing
#' the user to click on an image
#' @param file an image file to use
#' @return the number of pixels per mm
#' @export
#' @author Henning Redestig
calibrateScale <- function(file) {
  im <- readImage(file)
  done <- FALSE
  while(!done) {
    plotimage(im)
    answer <- readline("choose the distance to indicate [default=80mm]: ")
    d <- ifelse(answer == "", 80, as.integer(answer))
    message("left-click on two points separated ", d, "mm")
    pp <- locator(n=2, type="l")
    npixels <- sqrt(diff(pp$x)^2 + diff(pp$y)^2)
    answer <- readline("try again? [default=n]: ")
    done <- ifelse(answer %in% c("", "No", "NO", "n", "nO"), TRUE, FALSE)
  }
  npixels / d
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
#' @param cores number of cores to use. Can be left \code{NA} for
#' using a separately registered parallel backend (see \code{ddply})
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
  meta <- updateMeta(..., meta)
  ddply(mf, "image", function(dd) {
    file <- unique(as.character(dd$image))
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
#' @param min_area minimum area of a plant. A box with a total
#' feature area less than this number will have area set to missing
#' value.
#' @param ... not used
#' @return phenodata object
#' @author Henning Redestig
postProcessPlatePhenodata <- function(phenodata, min_area=1.5, ...) {
  ddply(phenodata, c("plate", "ROW", "RANGE"), function(dd) {
    tmp_area <- dd$AREA
    tmp_area[is.na(tmp_area)] <- Inf
    dd$too_small <- dd$too_small | tmp_area < min_area
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

rotationSigma <- function(dat, center, dp, r, sigma=10, eps=1e-2) {
  mu0 <- center - dp / 2 - ((r / 2) - 1) * dp
  mu <- sapply(1:r, function(i) mu0 + (i - 1) * dp)
  lambda <- rep(1/r, r)
  n <- length(dat)
  change <- 1
  while(change > eps) {
    psigma <- sigma
    sigmavec <- rep(sigma, r)
    zc <- .Call("phenotyping_normalpostp", dat, mu, sigmavec, lambda,
                PACKAGE=PKG)
    sigma <-
        sqrt(sum(sapply(1:r,
                        function(i) sum(zc$scaled_resp[,i] * (dat - mu[i])^2))) / n)
    change <- (psigma - sigma)^2
  }
  sigma
}

rotationSigmaFun <- function(par, im, pxlsmm, d, r) {
  dp <- pxlsmm * d
  rim <- rotate(im, par[[1]])
  rim <- roundImage(rim)
  center_x <- floor(nrow(rim) / 2)
  center_y <- floor(ncol(rim) / 2)
  data_x <- row(rim)[rim > 0]
  data_y <- col(rim)[rim > 0]
  sigmax <- rotationSigma(data_x, center_x, dp, r)
  sigmay <- rotationSigma(data_y, center_y, dp, r)
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

#' Plot an image using R-graphics device for easy annotation.
#'
#' Copied from the now defunct biOps package.
#' @param im an Image object or a character path to an image file
#' @param ... passed on to \code{image}
#' @return Nothing, used for the side effect
#' @export
#' @seealso \code{display} which is much faster and
#' better for just showing a picture but does not work well in vignettes etc
#' @examples
#' im <- system.file("examples/plate_merged.jpg", package="ezplate")
#' plotimage(im)
#' @author Henning Redestig
plotimage <- function(im, ...) {
  if(is.character(im))
    im <- readImage(im)
  colvec <- switch(ifelse(colorMode(im) == 0, "grey", "rgb"),
                   grey=grey(imageData(im)),
                   rgb=rgb(imageData(im)[,,1],
                   imageData(im)[,,2],
                   imageData(im)[,,3]))
  colors <- unique(colvec)
  colmat <- t(array(match(colvec, colors), dim = dim(im)[1:2]))
  image(x=0:(dim(colmat)[2]), y=0:(dim(colmat)[1]),
        z=t(colmat[nrow(colmat):1, ]),
        col=colors, xlab="", ylab="", axes=FALSE,
        asp=1, ...)
}

#' Label wells on image
#'
#' Create an image with highlighted well perimeters and labels in the middle of
#' each well indicating the well indices.
#' @param im an image
#' @param griddf a data frame detailing the grid
#' @return an image with highlighted well perimeters and indices
#' @export 
#' @author Henning Redestig
labelGridImage <- function(im, griddf, wellWidth) {
  imageData(im) <- colorWellFrames(imageData(im), griddf, wellWidth)
  font <- drawFont(size=15)
  xy <- as.matrix(cbind(griddf$centerx - wellWidth * 0.3, griddf$centery))
  drawText(im, xy=xy,
           labels=paste(griddf$row, griddf$col, sep=":"),
           font=font, col="white")
}

#' Create a font for drawing text on images
#'
#' Adopted from EBImage version 3.13.1 since it was removed in version 4
#' @param family A character value indicating the font family to
#' use. Valid examples on Linux/UNIX systems include \code{helvetica},
#' \code{times}, \code{courier} and \code{symbol}. Valid examples on
#' Windows machines include TrueType like \code{Arial} and
#' \code{Verdana}.
#' @param style A character value specifying the font style to use.
#' Supported styles are: \code{normal} (default), \code{italic}, and
#' \code{oblique}.
#' @param size Font size in points.
#' @param weight A numeric value indicating the font weight (bold
#' font). Supported values range between 100 and 900.
#' @param antialias A logical value indicating whether the font should
#' be anti-aliased.
#' @return An \code{Image} object or an array, containing the
#' transformed version of \code{img}.
#' @export
#' @seealso \code{\link{drawText}}
#' @author Oleg Sklyar and copied by Henning Redestig
drawFont <- function(family=switch(.Platform$OS.type, windows="Arial",
                       "helvetica"),
                     style="n", size=14, weight=200, antialias=TRUE) {
  ## drawFont
  res <- list(family=family, style=style, size=size, weight=weight, antialias=antialias)
  class(res) <- "DrawFont"
  res
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
#' lena <- readImage(system.file("images", "lena-color.png", package="EBImage"))
#' font <- drawFont(weight=600, size=28)
#' lena <- drawText(lena, xy=c(250, 450), labels="Lena", font=font, col="white")
#' display(lena, method="raster")
#' @author Oleg Sklyar and copied by Henning Redestig
drawText <- function(img, xy, labels, font, col) {
  ## drawText
  if(is.numeric(xy)) {
    xy <- list(as.numeric(xy))
    labels <- list(labels)
  }
  if(missing(font)) font <- drawFont()
  if(missing(col)) col <- "white"
  
  if (length(xy) != length(labels) || length(xy) != getNumberOfFrames(img, "render"))
    stop("lists of coordinates 'xy' labels 'labels' must",
         " be of the same length as the number of render frames")
  xy <- lapply(xy, as.numeric)
  for(i in seq_along(labels)) 
    if(!is.character(labels[[i]]))
      stop("all elements of 'labels' must be of class 'character'")
  if(!is(font, "DrawFont") )
    stop("to set the font use the 'drawFont' function which ",
         "returns an S3 class 'DrawFont', modify the slots as needed")

  font$style <- as.integer(switch(tolower(substr(font$style,1,1)), i=1, o=2, 0))
  font$size <- as.numeric(font$size)
  font$weight <- as.numeric(font$weight)
  font$antialias <- as.logical(font$antialias)
  return(.Call("lib_drawText", castImage(img), xy, labels, font, col,
               PACKAGE=PKG))
}

castImage <- function(x) {
  if (storage.mode(imageData(x)) != "double")
    storage.mode(imageData(x)) <- "double"
  x
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
