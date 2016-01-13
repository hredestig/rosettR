# rosettR

rosettR is a protocol for a high-throughput growth phenotyping assay
for *Arabidopsis* implemented as an R-package that details all steps
and provides data analysis with minimal required user interaction.

[Main project website.](http://hredestig.github.io/rosettR)

## Installation

This package is not yet on CRAN so installation is best done using the
[devtools package](https://github.com/hadley/devtools). First, install
EBImage which is a dependency from bioconductor that needs to be
installed separately. In R do:

```R
source("http://www.bioconductor.org/biocLite.R")
biocLite("EBImage")
```

Then, after installing devtools, do

```R
library(devtools)
install_github("hredestig/rosettR", subdir="pkg")
```

Then get started by opening the introduction vignette and following
the examples provided there.

```R
library(rosettR)
vignette("introduction", "rosettR")
```
