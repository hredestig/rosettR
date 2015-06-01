#!/bin/env Rscript

exampleMetadata <-
  metaTemplate(name="6x6.abcd",
               germplasms=letters[1:12],
               nrepeats=5,
               timepoints=c(11, 14, 16, 18),
               treatments=c("control", "stress"))

save(exampleMetadata, file="../pkg/data/exampleMetadata.rda")


