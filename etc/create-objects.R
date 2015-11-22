#!/bin/env Rscript

exampleMetadata <-
  metaTemplate(name="6x6.abcd",
               genotypes=c("foo", "bar", "baz", "qux"),
               nblocks=3,
               timepoints=c(11, 14, 16, 18),
               treatments=c("control", "stress"))

save(exampleMetadata, file="pkg/data/exampleMetadata.rda")


