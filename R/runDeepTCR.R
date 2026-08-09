#' Featurize TCR Sequences with the DeepTCR Variational Autoencoder
#'
#' @description Extracts CDR3 amino acid sequences from a
#'   SingleCellExperiment object and learns a low-dimensional representation
#'   using DeepTCR's unsupervised variational autoencoder (DeepTCR_U).
#'
#' @param input A SingleCellExperiment object containing scRepertoire TCR data.
#' @param chains Which chain(s) to featurize: "TRB", "TRA", or "both".
#'   Default is "TRB".
#' @param latent_dim Width of the VAE latent space. Default is 256.
#' @param use_genes Logical. If TRUE, includes V and J gene usage alongside the
#'   CDR3 sequence. Default is FALSE, sequence only.
#' @param stop_criterion Training convergence threshold passed to Train_VAE.
#'   Default is 0.01.
#' @param combine_chains Logical. If TRUE and chains="both", concatenates alpha
#'   and beta sequences with "_". Default is FALSE.
#' @param reduction_name Name for the dimensional reduction.
#'   Default is "tcr_deeptcr".
#' @param reduction_key Key prefix for the reduction columns.
#'   Default is "DeepTCR_".
#' @param return_object Logical. If TRUE, adds features as a dimensional
#'   reduction. If FALSE, returns the raw features. Default is TRUE.
#' @param seed Optional integer passed to Train_VAE as `graph_seed` for
#'   reproducible training. Default is NULL.
#' @param verbose Logical. If FALSE, suppresses DeepTCR's training output.
#'   Default is TRUE.
#' @param ... Additional arguments passed to calculate.deepTCR().
#'
#' @return If return_object=TRUE, the input object with features added as the
#'   `reduction_name` reduction. Cells with no sequence for the requested chain
#'   get NA. If return_object=FALSE, a list with `features` (one row per unique
#'   sequence), `sequences`, `barcodes`, and `explained_variance_ratio`.
#'
#' @details The VAE is trained on unique sequences rather than on cells.
#'   DeepTCR's own documentation notes that `Load_Data` does not merge
#'   identical amino acid sequences, so passing a redundant repertoire would
#'   both waste training time and weight the loss toward expanded clones.
#'   Features are expanded back onto cells afterward, so two cells sharing a
#'   CDR3 share a feature vector.
#'
#'   DeepTCR runs in its own basilisk environment (`deepTCREnv`) because it
#'   pins TensorFlow 2.12 and an older numpy, pandas, and scipy stack that
#'   cannot coexist with the main immLynx environment. The first call builds
#'   that environment, which is a large download.
#'
#'   `latent_dim` sets the requested width. DeepTCR prunes uninformative latent
#'   features during training, so the returned matrix may be narrower.
#'
#' @export
#' @importFrom methods is
#' @importFrom SummarizedExperiment colData colData<-
#'
#' @examples
#' data(immLynx_example)
#' \donttest{
#'   # Learn a 64-dimensional representation of the TRB repertoire
#'   sce <- runDeepTCR(immLynx_example,
#'                     chains = "TRB",
#'                     latent_dim = 64)
#'
#'   # Include V and J gene usage
#'   sce <- runDeepTCR(immLynx_example,
#'                     chains = "TRB",
#'                     use_genes = TRUE)
#'
#'   # Get the raw feature matrix
#'   res <- runDeepTCR(immLynx_example,
#'                     chains = "TRB",
#'                     return_object = FALSE)
#'   dim(res$features)
#' }
runDeepTCR <- function(input,
                       chains = c("TRB", "TRA", "both"),
                       latent_dim = 256,
                       use_genes = FALSE,
                       stop_criterion = 0.01,
                       combine_chains = FALSE,
                       reduction_name = "tcr_deeptcr",
                       reduction_key = "DeepTCR_",
                       return_object = TRUE,
                       seed = NULL,
                       verbose = TRUE,
                       ...) {

  chains <- match.arg(chains)

  if (!methods::is(input, "SingleCellExperiment")) {
    stop("Input must be a SingleCellExperiment object")
  }

  if (!is.numeric(latent_dim) || length(latent_dim) != 1L ||
      is.na(latent_dim) || latent_dim < 1) {
    stop("latent_dim must be a single positive integer")
  }
  latent_dim <- as.integer(latent_dim)

  message("Extracting TCR sequences from object...")
  seq_map <- .extractChainSeqs(input, chains = chains,
                               combine_chains = combine_chains,
                               with_genes = use_genes)

  keep <- !duplicated(seq_map$sequences)
  unique_seqs <- seq_map$sequences[keep]

  if (length(unique_seqs) == 0L) {
    stop("No CDR3 sequences found for the requested chain(s).")
  }

  message("Training DeepTCR VAE on ", length(unique_seqs),
          " unique sequences (latent_dim = ", latent_dim, ")...")

  res <- calculate.deepTCR(
    sequences = unique_seqs,
    v_genes   = if (use_genes) seq_map$v[keep] else NULL,
    j_genes   = if (use_genes) seq_map$j[keep] else NULL,
    latent_dim = latent_dim,
    stop_criterion = stop_criterion,
    seed = seed,
    verbose = verbose,
    ...
  )

  features <- res$features
  if (is.null(dim(features))) {
    features <- matrix(features, nrow = length(unique_seqs))
  }
  message("Learned ", ncol(features), " latent feature(s).")

  if (!return_object) {
    return(list(
      features = features,
      sequences = unique_seqs,
      barcodes = seq_map$barcodes,
      explained_variance_ratio = res$explained_variance_ratio
    ))
  }

  # Expand unique-sequence features onto the cells that carry them.
  cell_features <- .deeptcrFeaturesToCells(features, unique_seqs,
                                           seq_map$sequences)

  full <- matrix(NA_real_,
                 nrow = ncol(input),
                 ncol = ncol(cell_features),
                 dimnames = list(colnames(input),
                                 paste0(reduction_key,
                                        seq_len(ncol(cell_features)))))
  full[seq_map$barcodes, ] <- cell_features

  SingleCellExperiment::reducedDim(input, reduction_name) <- full

  message("Features added as '", reduction_name, "' reduction")
  input
}


#' Expand per-sequence features onto cells
#'
#' @param features Matrix with one row per unique sequence.
#' @param unique_seqs Character vector matching the rows of `features`.
#' @param cell_seqs Character vector of the sequence carried by each cell.
#' @return A matrix with one row per entry of `cell_seqs`.
#' @keywords internal
.deeptcrFeaturesToCells <- function(features, unique_seqs, cell_seqs) {

  if (nrow(features) != length(unique_seqs)) {
    stop("features must have as many rows as unique sequences")
  }

  idx <- match(cell_seqs, unique_seqs)
  if (anyNA(idx)) {
    stop("every cell sequence must appear in the unique sequence set")
  }

  features[idx, , drop = FALSE]
}
