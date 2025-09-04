# immLynx
Linking python packages and hugging face models for immune repertoire analysis in R

<img align="right" src="https://github.com/BorchLab/immLynx/blob/main/www/immlynx_hex.png" width="305" height="352">

## Introduction

immLynx provides a suite of functions to interface with popular Python-based tools for immune repertoire analysis. It handles the Python environment setup automatically using the `basilisk` package, so you can call these powerful tools directly from R without worrying about dependencies.

## Features

immLynx currently supports the following Python packages:

*   **tcrdist3**: Calculate pairwise distances between T-cell receptors using `calculate.tcrDist`.
*   **DeepTCR**: Perform unsupervised feature extraction from TCR sequences with `calculate.deepTCR`.
*   **OLGA**: Compute the generation probability of CDR3 sequences or generate new sequences with `calculate.olga`.
*   **soNNia**: Infer selection pressures on TCRs using `calculate.sonia`.
*   **clusTCR**: Cluster large sets of CDR3 sequences with `calculate.clustcr`.

For more details on each function, please refer to the R documentation (e.g., `?calculate.tcrDist`).

## System requirements 

immLynxhas been tested on R versions >= 4.0. Please consult the DESCRIPTION file for more details on required R packages - it is specifically designed to work with single-cell objects that have had BCR/TCRs added using [scRepertoire](https://github.com/BorchLab/scRepertoire). immLynx has been tested on OS X and Linux platforms.

## Installation

To run immLynx, open R and install immApex from github: 

```r
devtools::install_github("BorchLab/immLynx")
```

## Quick Start 

Here is a quick example of how to use `immLynx` to cluster CDR3 sequences with `clusTCR`.

```r
library(immLynx)

# Sample CDR3 sequences
seqs <- c("CASSLAGGREQYF", "CASSLSFGREQYF", "CASSIWSGREQYF", "CASSLGGRYNEQFF")

# Cluster the sequences using the MCL algorithm
clusters <- calculate.clustcr(sequences = seqs, method = "mcl")

# View the results
print(clusters)
```

The first time you run a function that uses a Python package, `immLynx` will automatically create a dedicated Conda environment with all the necessary dependencies. This might take a few minutes, but it only happens once.

## Bug Reports/New Features

#### If you run into any issues or bugs please submit a [GitHub issue](https://github.com/BorchLab/immLynx/issues) with details of the issue.

- If possible please include a [reproducible example](https://reprex.tidyverse.org/). 

#### Any requests for new features or enhancements can also be submitted as [GitHub issues](https://github.com/BorchLab/immLynx/issues).

#### [Pull Requests](https://github.com/BorchLab/immLynx/pulls) are welcome for bug fixes, new features, or enhancements.
