#' immLynx: Advanced Analysis Pipelines for single-cell immune repertoire data
#'
#' @description
#' immLynx provides a unified interface for running multiple analysis pipelines
#' on single-cell sequencing data in the scRepertoire format. The package wraps
#' several popular Python-based tools and makes them accessible within R workflows.
#'
#' @section Main Functions:
#'
#' **Clustering & Grouping:**
#' \itemize{
#'   \item \code{\link{runClustTCR}}: Cluster TCRs using clusTCR (MCL or DBSCAN)
#'   \item \code{\link{runMetaclonotypist}}: Identify metaclones with metaclonotypist
#' }
#'
#' **Distance Metrics:**
#' \itemize{
#'   \item \code{\link{runTCRdist}}: Calculate TCR distances with tcrdist3
#' }
#'
#' **Generation Probability:**
#' \itemize{
#'   \item \code{\link{runOLGA}}: Calculate Pgen using OLGA
#'   \item \code{\link{generateOLGA}}: Generate random TCR sequences
#' }
#'
#' **Embeddings:**
#' \itemize{
#'   \item \code{\link{runEmbeddings}}: Generate protein embeddings using ESM-2
#' }
#'
#' **Selection Inference:**
#' \itemize{
#'   \item \code{\link{runSoNNia}}: Infer selection with soNNia
#' }
#'
#' **Utility Functions:**
#' \itemize{
#'   \item \code{\link{extractTCRdata}}: Extract TCR data from Seurat objects
#'   \item \code{\link{validateTCRdata}}: Validate TCR data format
#'   \item \code{\link{convertToTcrdist}}: Convert to tcrdist3 format
#'   \item \code{\link{summarizeTCRrepertoire}}: Generate repertoire statistics
#' }
#'
#' @section Data Format:
#' All functions expect Seurat objects with scRepertoire TCR data stored in the metadata.
#' The functions use \code{immApex::getIR()} to extract TCR sequences in the format:
#' \itemize{
#'   \item \code{barcode}: Cell barcode
#'   \item \code{cdr3_aa}: CDR3 amino acid sequence
#'   \item \code{v}: V gene
#'   \item \code{d}: D gene (if applicable)
#'   \item \code{j}: J gene
#'   \item \code{c}: C gene
#'   \item \code{chain}: Chain type (TRA/TRB)
#' }
#'
#' @section Python Dependencies:
#' Python packages are managed automatically by basilisk. The following are included:
#' \itemize{
#'   \item clustcr - TCR clustering
#'   \item tcrdist3 - TCR distance calculations
#'   \item olga - Generation probability
#'   \item sonnia - Selection inference
#'   \item metaclonotypist - Metaclone discovery
#'   \item transformers - Hugging Face models
#'   \item torch - PyTorch for GPU support
#' }
#'
#' @seealso
#' Useful links:
#' \itemize{
#'   \item scRepertoire: \url{https://github.com/BorchLab/scRepertoire}
#'   \item immApex: \url{https://github.com/BorchLab/immApex}
#' }
#'
#' @docType package
#' @name immLynx-package
#' @aliases immLynx
#' @keywords internal
"_PACKAGE"
