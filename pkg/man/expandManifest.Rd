% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rosettr-utils.R
\name{expandManifest}
\alias{expandManifest}
\title{Expand experiment design}
\usage{
expandManifest(meta, plateName = "plate\%03d.jpg", plateOffset = 0, ...)
}
\arguments{
\item{meta}{the meta data object that defines the experiment
design; the number of repeats, timepoints, germplasms and
treatments}

\item{plateName}{canonical file name of the plates to use with
\code{\link{sprintf}}.}

\item{plateOffset}{the number that the plate numbering should be
offset to. The number of the first plate is 1 + plateoffset.}

\item{...}{not used}
}
\value{
a data frame representing the suggested design of the
experiment
}
\description{
Given the parameters defined in the metadata, expand the balanced
design matrix with each germplasm and each treatment represented
in each repeat of the experiment.
}
\details{
If you want to be able to regenerate the exact same randomize
allocation of genotypes to plates, then make sure to set the
randomization seed prior to calling this function.
}
\examples{
data(exampleMetadata)
exampleMetadata
set.seed(123) # for reproducibility
expandManifest(exampleMetadata)
}
\author{
Henning Redestig
}

