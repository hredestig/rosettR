library(EBImage)
test_that("experiment can be created and sowing report compiled", {
  expect_that({
    path <- makeTestExperiment()
    makeReport(path, "layout", quiet=TRUE)
    file.exists(file.path(path, "Output/layout/sowing-summary.csv"))
  }, is_true())
})

test_that("a plate rotation can be compensated", {
  expect_that({
    file <- system.file("examples", "plate.jpg", package="rosettR")
    pixelsmm <- 3.7454
    im <- resize(rosettR:::makeGrey(readImage(file), channels=3), 500)
    binary <- rosettR:::makeBinary(im, pixelsmm / 2)
    binaryRotated <- rotate(binary, 4)
    origAngle <- rosettR:::optimAngle(binary, pixelsmm / 2, 20, 6)
    estAngle <- rosettR:::optimAngle(binaryRotated, pixelsmm / 2, 20, 6,
                                     c(-5,5))
    binaryCorrected <- rotate(binaryRotated, estAngle)
    ## display(binaryCorrected, method="raster")
    ## abline(v=260, col="red")
    ## curiosly some difference in results linux / windows
    abs(4 + origAngle + estAngle) < 2
  }, is_true())
})

test_that("major plate dislocation can be compensated", {
  expect_that({
    file <- system.file("examples", "plate.jpg", package="rosettR")
    largeRgb <- EBImage::readImage(file)
    largeRgb <- rosettR:::makeSquare(largeRgb)
    largeGrey <- rosettR:::makeGrey(largeRgb, 3)
    smallGrey <- EBImage::resize(largeGrey, 500)
    smallGreyCut <- smallGrey[-(1:30),]
    scaling <- 500 / nrow(largeRgb)
    pixelsmmSmall <- 3.7454 * scaling
    radius <- 76.5 * pixelsmmSmall
    loc <- rosettR:::findPlate(smallGreyCut, radius, 0.95)
    ## marked <- markCircle(smallGreyCut,
    ##                      0.5 * nrow(smallGreyCut) - loc$deltax / 2,
    ##                      0.5 * ncol(smallGreyCut) - loc$deltay / 2,
    ##                      radius,
    ##                      edgeValue=0)
    ## display(marked, method="raster")
    abs(loc$deltax - 24 + loc$deltay  + 2) < 4
  }, is_true())
})


test_that("major plate dislocation can be compensated", {
  expect_that({
    file <- system.file("examples", "plate.jpg", package="rosettR")
    largeRgb <- EBImage::readImage(file)
    largeRgb <- rosettR:::makeSquare(largeRgb)
    largeGrey <- rosettR:::makeGrey(largeRgb, 3)
    smallGrey <- EBImage::resize(largeGrey, 500)
    smallGreyCut <- smallGrey[-(1:30),]
    scaling <- 500 / nrow(largeRgb)
    pixelsmmSmall <- 3.7454 * scaling
    radius <- 76.5 * pixelsmmSmall
    loc <- rosettR:::findPlate(smallGreyCut, radius, 0.95)
    ## marked <- markCircle(smallGreyCut,
    ##                      0.5 * nrow(smallGreyCut) - loc$deltax / 2,
    ##                      0.5 * ncol(smallGreyCut) - loc$deltay / 2,
    ##                      radius,
    ##                      edgeValue=0)
    ## display(marked, method="raster")
    abs(loc$deltax - 24 + loc$deltay  + 2) < 4
  }, is_true())
})

test_that("plant areas can be estimated", {
  expect_that({
    file <- system.file("examples", "plate.jpg", package="rosettR")
    meta <- metaTemplate(letters[1:4], LETTERS[1:2], reference="a")
    df <- analyzeImage(file, meta$griddf, 3.7454, 20, 6, 75, verbose=TRUE)
    answer <- c(0.43, 25.73, 23.03, 28.66, 17.89, 2.92, 18.46, 13.19, 18.25, 
                21.31, 10.69, 10.91, 16.75, 16.04, 23.88, 14.61, 20.67, 18.75, 
                19.32, 16.61, 25.95, 26.95, 22.17, 21.74, 19.18, 18.18, 24.52, 
                17.61, 21.81, 20.32, 18.32, 25.09)
    cor(df$AREA, answer) > 0.95
  }, is_true())
})
