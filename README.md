# immLynx
Linking advanced TCR python pipelines and Hugging Face models in R

<!-- badges: start -->
  [![R-CMD-check](https://github.com/BorchLab/immLynx/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/BorchLab/immLynx/actions/workflows/R-CMD-check.yaml)
  [![Codecov test coverage](https://codecov.io/gh/BorchLab/immLynx/graph/badge.svg)](https://app.codecov.io/gh/BorchLab/immLynx)
<!-- badges: end -->


<img align="right" src="https://github.com/BorchLab/immLynx/blob/main/www/immlynx_hex.png" width="305" height="352">

immLynx provides a unified R interface for running multiple state-of-the-art TCR analysis
pipelines on single-cell TCR sequencing data. The package seamlessly integrates
with Seurat and scRepertoire workflows, wrapping popular Python-based tools to enable:

*   **tcrdist3**: Calculate pairwise distances between T-cell receptors using `runTCRdist`
*   **OLGA**: Compute the generation probability of CDR3 sequences or generate new sequences with `runOLGA`
*   **soNNia**: Infer selection pressures on TCRs using `runSoNNia`
*   **clusTCR**: Cluster large sets of CDR3 sequences with `runClustTCR`
*   **metaclonotypist**: Identify TCR metaclones with `runMetaclonotypist`
*   **ESM-2**: Generate protein language model embeddings with `runEmbeddings`

For more details on each function, please refer to the R documentation (e.g., `?runTCRdist`).

## System Requirements

immLynx has been tested on R versions >= 4.3. Please consult the DESCRIPTION file
for more details on required R packages - it is specifically designed to work with
single-cell objects that have had BCR/TCRs added using
[scRepertoire](https://github.com/BorchLab/scRepertoire). immLynx has been tested
on OS X and Linux platforms.

## Installation

**Install immLynx:**
```r
# Install from GitHub
remotes::install_github("BorchLab/immLynx")
```

The first time you use immLynx, it will automatically install the required Python
packages in an isolated environment via basilisk. This may take several minutes.

## Quick Start

```r
library(immLynx)
library(Seurat)

# Load example data
data("immLynx_example")

# Summarize TCR repertoire
summary <- summarizeTCRrepertoire(immLynx_example)
print(summary)

# Cluster TCRs
seurat_obj <- runClustTCR(immLynx_example, chains = "TRB", method = "mcl")

# Calculate generation probability
seurat_obj <- runOLGA(seurat_obj, chains = "TRB")

# Generate protein embeddings
seurat_obj <- runEmbeddings(seurat_obj, chains = "TRB")

# Visualize embeddings
seurat_obj <- RunUMAP(seurat_obj, reduction = "tcr_esm", dims = 1:30)
DimPlot(seurat_obj, reduction = "umap")
```

## Main Functions

### Utility Functions

Extract and validate TCR data:

```r
# Extract TCR data
tcr_data <- extractTCRdata(seurat_obj, chains = "TRB")

# Validate data format
validation <- validateTCRdata(tcr_data)

# Convert to tcrdist3 format
tcrdist_format <- convertToTcrdist(tcr_data)

# Generate repertoire summary
summary <- summarizeTCRrepertoire(seurat_obj)
```

### TCR Clustering

Cluster TCRs based on sequence similarity using clusTCR:

```r
# MCL clustering (default)
seurat_obj <- runClustTCR(seurat_obj,
                          chains = "TRB",
                          method = "mcl",
                          inflation = 2.0)

# DBSCAN clustering
seurat_obj <- runClustTCR(seurat_obj,
                          chains = "TRB",
                          method = "dbscan",
                          eps = 0.5)
```

### Metaclone Discovery

Identify metaclones using metaclonotypist:

```r
# Run metaclonotypist with TCRdist
seurat_obj <- runMetaclonotypist(seurat_obj,
                                  chains = "beta",
                                  method = "tcrdist",
                                  max_edits = 2,
                                  max_dist = 20)

# Use SCEPTR distance metric
seurat_obj <- runMetaclonotypist(seurat_obj,
                                  method = "sceptr",
                                  max_dist = 1.0)
```

### Generation Probability

Calculate how likely each TCR sequence is to be generated naturally:

```r
# Calculate Pgen for TRB sequences
seurat_obj <- runOLGA(seurat_obj,
                      chains = "TRB",
                      model = "humanTRB")

# Generate random TCR sequences
random_tcrs <- generateOLGA(n = 1000, model = "humanTRB")
```

### Protein Embeddings

Generate dense vector representations using ESM-2:

```r
# Default: ESM-2 35M model
seurat_obj <- runEmbeddings(seurat_obj,
                            chains = "TRB",
                            pool = "mean")

# Use larger model for better embeddings
seurat_obj <- runEmbeddings(seurat_obj,
                            model_name = "facebook/esm2_t33_650M_UR50D")

# Visualize in UMAP space
seurat_obj <- RunUMAP(seurat_obj, reduction = "tcr_esm", dims = 1:30)
DimPlot(seurat_obj, reduction = "umap")
```

### TCR Distance Calculation

Compute pairwise distances between TCRs:

```r
# Calculate TRB distances
dist_results <- runTCRdist(seurat_obj,
                           chains = "beta",
                           organism = "human")

# Access distance matrices
beta_dist <- dist_results$distances$pw_beta
```

### Selection Inference

**soNNia Selection:**
```r
# 1. Generate background
background <- generateOLGA(n = 10000, model = "humanTRB")
write.csv(background, "background.csv", row.names = FALSE)

# 2. Run soNNia
seurat_obj <- runSoNNia(seurat_obj,
                        background_file = "background.csv")
```

## Data Format

immLynx expects Seurat objects with scRepertoire TCR data in the metadata.
The data should include columns like:

- `CTgene`: Gene information (V/J genes)
- `CTaa`: CDR3 amino acid sequences
- `CTnt`: CDR3 nucleotide sequences
- `CTstrict`: Combined TCR information

The package uses `immApex::getIR()` internally to extract:
- `barcode`: Cell identifier
- `cdr3_aa`: CDR3 amino acid sequence
- `v`, `d`, `j`, `c`: Gene segments
- `chain`: TRA or TRB

## Citation

If you use immLynx in your research, please cite the underlying tools:

- **clusTCR**: [Valkiers et al. (2021)](https://pubmed.ncbi.nlm.nih.gov/34132766/)
- **tcrdist3**: [Mayer-Blackwell et al. (2021)](https://pubmed.ncbi.nlm.nih.gov/36087210/)
- **OLGA**: [Sethna et al. (2019)](https://pubmed.ncbi.nlm.nih.gov/30657870/)
- **soNNia**: [Isacchini et al. (2021)](https://pubmed.ncbi.nlm.nih.gov/33795515/)
- **metaclonotypist**: [qimmuno](https://github.com/qimmuno/metaclonotypist)
- **ESM-2**: [Lin et al. (2023)](https://pubmed.ncbi.nlm.nih.gov/36927031/)

## Related Packages

If you are interested in GLIPH2 specificity group analysis, please see [immGLIPH](https://github.com/BorchLab/immGLIPH), a dedicated R package for running GLIPH:

```r
remotes::install_github("BorchLab/immGLIPH")
```

## Bug Reports/New Features

#### If you run into any issues or bugs please submit a [GitHub issue](https://github.com/BorchLab/immLynx/issues) with details of the issue.

- If possible please include a [reproducible example](https://reprex.tidyverse.org/).

#### Any requests for new features or enhancements can also be submitted as [GitHub issues](https://github.com/BorchLab/immLynx/issues).

#### [Pull Requests](https://github.com/BorchLab/immLynx/pulls) are welcome for bug fixes, new features, or enhancements.
