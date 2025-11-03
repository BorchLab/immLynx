# immLynx
Linking advanced TCR python pipelines and Hugging Face models in R

<img align="right" src="https://github.com/BorchLab/immLynx/blob/main/www/immlynx_hex.png" width="305" height="352">

## Introduction

provides a unified R interface for running multiple state-of-the-art TCR analysis 
pipelines on single-cell TCR sequencing data. The package seamlessly integrates 
with Seurat and scRepertoire workflows, wrapping popular Python-based tools to enable:

*   **tcrdist3**: Calculate pairwise distances between T-cell receptors using `runTCRdist`.
*   **DeepTCR**: Perform unsupervised feature extraction from TCR sequences with `runDeepTCR`.
*   **OLGA**: Compute the generation probability of CDR3 sequences or generate new sequences with `runOLGA`.
*   **soNNia**: Infer selection pressures on TCRs using `runSoNNia`.
*   **clusTCR**: Cluster large sets of CDR3 sequences with `runClusTCR`.

For more details on each function, please refer to the R documentation (e.g., `?runTCRdist`).

## System requirements 

immLynx has been tested on R versions >= 4.0. Please consult the DESCRIPTION file 
for more details on required R packages - it is specifically designed to work with 
single-cell objects that have had BCR/TCRs added using 
[scRepertoire](https://github.com/BorchLab/scRepertoire). immLynx has been tested 
on OS X and Linux platforms.

## Installation

**Install immLynx:**
```r
remotes::install_github("yourusername/immLynx")
```

The first time you use immLynx, it will automatically install the required Python 
packages in an isolated environment. This may take several minutes.

## Main Functions

### TCR Clustering

Cluster TCRs based on sequence similarity using clusTCR:

```r
# MCL clustering (default)
seurat_obj <- runClusTCR(seurat_obj, 
                         chains = "TRB", 
                         method = "mcl",
                         inflation = 2.0)

# DBSCAN clustering
seurat_obj <- runClusTCR(seurat_obj,
                         chains = "TRB",
                         method = "dbscan",
                         eps = 0.5,
                         min_samples = 5)

# Cluster paired chains
seurat_obj <- runClusTCR(seurat_obj,
                         chains = "both",
                         combine_chains = TRUE)
```

### Generation Probability

Calculate how likely each TCR sequence is to be generated naturally:

```r
# Calculate Pgen for TRB sequences
seurat_obj <- runOLGA(seurat_obj, 
                      chains = "TRB",
                      model = "humanTRB")

# Include V/J gene information
seurat_obj <- runOLGA(seurat_obj,
                      chains = "TRB",
                      use_vj_genes = TRUE)

# Generate random TCR sequences
random_tcrs <- generateOLGA(n = 1000, model = "humanTRB")
```

### Protein Embeddings

Generate dense vector representations using protein language models:

```r
# Default: ESM-2 35M model
seurat_obj <- runEmbeddings(seurat_obj, 
                            chains = "TRB",
                            chunk_size = 32)

# Use larger model for better embeddings
seurat_obj <- runEmbeddings(seurat_obj,
                            chains = "TRB",
                            model_name = "facebook/esm2_t33_650M_UR50D")

# Visualize in UMAP space
seurat_obj <- RunUMAP(seurat_obj, reduction = "tcr_esm", dims = 1:30)
DimPlot(seurat_obj, reduction = "umap", group.by = "seurat_clusters")
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
cdr3_dist <- dist_results$distances$pw_cdr3_b_aa

# Perform hierarchical clustering on distances
hc <- hclust(as.dist(cdr3_dist))
plot(hc)
```

### Advanced: Deep Learning and Selection

**DeepTCR** (requires pre-formatted data files):
```r
# Run DeepTCR VAE
features <- runDeepTCR(seurat_obj,
                       output_dir = "deeptcr_output",
                       latent_dim = 100,
                       epochs = 100)
```

**soNNia Selection** (requires background sequences from OLGA):
```r
# 1. Generate background
background <- generateOLGA(n = 10000, model = "humanTRB")
write.csv(background, "background_tcrs.csv")

# 2. Run soNNia
results <- runSonia(seurat_obj,
                    background_file = "background_tcrs.csv",
                    organism = "human")
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

If you use immLynx in your research, please cite the tools you are using in the package:

And the underlying tools:
- **clusTCR**: [Valkiers et al. (2021)](https://pubmed.ncbi.nlm.nih.gov/34132766/)
- **tcrdist3**: [Mayer-Blackwell et al. (2021)](https://pubmed.ncbi.nlm.nih.gov/36087210/)
- **OLGA**: [Sethna et al. (2019)](https://pubmed.ncbi.nlm.nih.gov/30657870/)
- **soNNia**: [Isacchini etd al. (2021)](https://pubmed.ncbi.nlm.nih.gov/33795515/)
- **DeepTCR**: [Sidhom et al. (2021)](https://pubmed.ncbi.nlm.nih.gov/33707415/)


## Bug Reports/New Features

#### If you run into any issues or bugs please submit a [GitHub issue](https://github.com/BorchLab/immLynx/issues) with details of the issue.

- If possible please include a [reproducible example](https://reprex.tidyverse.org/). 

#### Any requests for new features or enhancements can also be submitted as [GitHub issues](https://github.com/BorchLab/immLynx/issues).

#### [Pull Requests](https://github.com/BorchLab/immLynx/pulls) are welcome for bug fixes, new features, or enhancements.
