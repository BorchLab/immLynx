# immLynx 1.3.1

* Added `runScXpand()`, which predicts T-cell clonal expansion from gene
  expression alone using scXpand's pretrained pan-cancer models. Unlike the
  other wrapped tools it needs no receptor sequences, so it works on
  datasets with no paired TCR sequencing. Inference only; training and
  hyperparameter optimization stay in Python.
* Added `listScXpandModels()` to enumerate the available pretrained models
  without building the Python environment.
* Added a dedicated `scXpandEnv` basilisk environment (Python 3.11, CPU
  PyTorch, scxpand 0.4.6). It cannot share the existing environments
  because scxpand needs Python 3.11 and torch 2.5. The first call builds
  several gigabytes and downloads the selected model from figshare.
* `runScXpand()` resolves gene identifiers to the Ensembl IDs scXpand's
  models are indexed by, optionally mapping symbols through `org.Hs.eg.db`.
  Ambiguous symbols (`HLA-DRA` alone maps to eight Ensembl IDs) and
  collapsed duplicates are counted and reported rather than resolved
  silently, because scXpand zero-fills genes it cannot find.
* When scRepertoire clone calls are present, `runScXpand()` derives
  `clone_id_size`, `median_clone_size` and `expansion` per sample using
  scXpand's 1.5x-median rule, so a gene-expression-only prediction can be
  scored against the observed repertoire. Clone sizes are tabulated fresh
  rather than read from `clonalFrequency`, which `combineExpression()`
  computes under whatever grouping was in effect.
* Pretrained models are cached under `tools::R_user_dir("immLynx", "cache")`.
  scXpand's own default would write a `.scxpand_cache` directory into the
  current working directory, so `runScXpand()` downloads the model as an
  explicit step with an explicit cache location.
* Worked around an upstream download bug: scXpand's registry points at
  `figshare.com/ndownloader/articles/...`, which answers HTTP 202 with an
  empty body, and the failed download is cached so retries keep failing.
  The same archive on `ndownloader.figshare.com` serves correctly, so
  `runScXpand()` rewrites the host. The rewrite becomes a no-op once
  upstream fixes its URLs.

# immLynx 1.1.2

* Added `exportToScanpy()` to write a `SingleCellExperiment` or `Seurat`
  object (with optional scRepertoire immune-receptor metadata) to
  scanpy/scirpy-compatible H5AD or H5MU files, plus an AIRR
  rearrangement TSV sidecar. H5AD writing uses `zellkonverter`; H5MU
  assembly uses scirpy/muon inside a dedicated `scanpyExportEnv`
  basilisk environment. BCR loci (`IGH`, `IGK`, `IGL`) are mapped to
  `immApex::getIR()`'s `Heavy`/`Light` chains and disambiguated by
  V-gene prefix.

# immLynx 1.0.0

* Official Bioconductor release in 3.23

# immLynx 0.99.4

* Added `\value` section to `TCR_summary-class` documentation to resolve
  R CMD check warning about empty/missing `\value` sections

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
