# TSAtools

## Overview

`TSAtools` is an R package developed by Omid Gholamalamdari from the Andrew Belmont lab at the University of Illinois at Urbana-Champaign. This package provides a suite of tools for working with TSA-seq BigWig files and other binned genomic data, facilitating efficient data processing and analysis in genomic studies.

## Installation

To install the latest version of `TSAtools` from GitHub, use the following commands in R:

```{r, eval = FALSE}
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("belmontLab/TSAtools")
```

## Usage

After installation, load TSAtools into your R session:

```{r}
library(TSAtools)
```

## Basic Functions
Here are some of the basic functions provided by TSAtools:

readBigWig(): Read data from BigWig files.
binGenomicData(): Bin genomic data for analysis.
plotTSA(): Plot TSA-seq data for visualization.

## License

This package is licensed under the MIT License - see the LICENSE.md file for details.
