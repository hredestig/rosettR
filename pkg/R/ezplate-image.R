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
  imageData(im) <- colorWellFrames(mat, df, wellWidth)
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
  ## draw_font
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
#' @param font A font object, returned by \code{\link{draw_font}}. If
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
  ## draw_text
  if (is.numeric(xy)) {
    xy <- list(as.numeric(xy))
    labels <- list(labels)
  }
  if (missing(font)) font <- draw_font()
  if (missing(col)) col <- "white"
  
  if (length(xy) != length(labels) || length(xy) != getNumberOfFrames(img,'render'))
    stop("lists of coordinates 'xy' labels 'labels' must",
         " be of the same length as the number of render frames")
  xy <- lapply(xy, as.numeric)
  for ( i in seq_along(labels)) 
    if (!is.character(labels[[i]]))
      stop("all elements of 'labels' must be of class 'character'")
  if ( !is(font, "DrawFont") )
    stop("to set the font use the 'draw_font' function which ",
         "returns an S3 class 'DrawFont', modify the slots as needed")

  font$style <- as.integer(switch(tolower(substr(font$style,1,1)), i=1, o=2, 0))
  font$size <- as.numeric(font$size)
  font$weight <- as.numeric(font$weight)
  font$antialias <- as.logical(font$antialias)
  return(.Call("lib_drawText", castImage(img), xy, labels, font, col,
               PACKAGE="ezplate"))
}

castImage <- function(x) {
  ## cast_image
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
  ## frame_box
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
  ## box_index
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
  ## color_box_frames
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
