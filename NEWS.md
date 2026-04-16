# immLynx 0.99.3

* Addressed Bioconductor reviewer feedback
* Updated R dependency to >= 4.6.0
* Converted `TCR_summary` from S3 class to formal S4 class with `show` method
* Added ORCID to `Authors@R`
* Replaced `\dontrun{}` with `\donttest{}` in all examples except `runSoNNia`
  (kept as `\dontrun{}` due to upstream soNNia/numpy incompatibility)
* Removed redundant `tryCatch`/`stop()` pattern in `huggingModel()`;
  replaced with `on.exit()` cleanup
* Extracted shared `.run_in_basilisk()` helper to reduce code repetition
  across `calculate.*` functions
* Added `@importFrom SummarizedExperiment colData` and removed explicit `::`
  qualifiers in `runOLGA`, `runClustTCR`, and `runEmbeddings`
* Added `n == 0L` edge-condition guard for `seq.int()` in `proteinEmbeddings()`
* Expanded vignette introductions with biological context and detailed prose
  descriptions for each code section
* Removed GitHub installation instructions from vignettes
* Added additional unit tests for S4 class structure, edge cases, and
  input validation

# immLynx 0.99.2

* Removed the umbrella roxygen block from `R/utils.R`
* Added `skip_on_bioc_build()` check directly inside `skip_if_no_python()` in `tests/testthat/helper-immLynx.R`

# immLynx 0.99.1

* Switched example data from Seurat to SingleCellExperiment object
* Replaced Seurat with scran/scater in vignettes
* Removed Seurat dependency; functions now use SingleCellExperiment exclusively
* Renamed `return_seurat` parameter to `return_input` in `runMetaclonotypist()`

# immLynx 0.99.0

* Initial Bioconductor submission
* Added `runClustTCR()` for TCR clustering via clusTCR
* Added `runTCRdist()` for pairwise TCR distance calculations
* Added `runOLGA()` and `generateOLGA()` for generation probability
* Added `runEmbeddings()` for protein language model embeddings
* Added `runMetaclonotypist()` for metaclone discovery
* Added `runHLAassociation()` for HLA-metaclone associations
* Added `runSoNNia()` for selection inference
* Added utility functions: `extractTCRdata()`, `validateTCRdata()`,
  `convertToTcrdist()`, `summarizeTCRrepertoire()`
* Added `huggingModel()`, `tokenizeSequences()`, and
  `proteinEmbeddings()` for custom embedding workflows
