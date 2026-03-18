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
