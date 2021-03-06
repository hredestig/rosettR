% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rosettr-reporting.R
\name{createPlateTestDf}
\alias{createPlateTestDf}
\title{Data frame for analysis of plant areas from a plate experiment}
\usage{
createPlateTestDf(df)
}
\arguments{
\item{df}{a data frame from \code{\link{processPlateImages}}}
}
\value{
a data frame with added RGR (all consecutive timepoints as
well as between the first and te last timepoint) and \code{dcast}d
to have timepoints as separate columns.
}
\description{
Create a data frame suitable for performing hypothesis testing
containing the relative measurements, RGR and extra factors.
}
\author{
Henning Redestig
}

