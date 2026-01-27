#' Calculate TCR Distances on scRepertoire Data
#'
#' @description This function extracts TCR sequences from a Seurat or
#'   SingleCellExperiment object with scRepertoire data and calculates pairwise
#'   TCR distances using tcrdist3.
#'
#' @param input A Seurat or SingleCellExperiment object containing scRepertoire TCR data.
#' @param chains Character vector specifying chains: "alpha", "beta", or c("alpha", "beta").
#'   Default is "beta".
#' @param organism Organism: "human" or "mouse". Default is "human".
#' @param compute_distances Logical. Whether to compute full distance matrix. Default TRUE.
#' @param add_to_object Logical. If TRUE, attempts to add distance matrix to object
#'   (stored in misc slot for Seurat or metadata for SCE). Default is FALSE due to large
#'   matrix size.
#'
#' @return A list containing:
#'   \item{distances}{Distance matrices (pw_alpha, pw_beta, pw_cdr3_a_aa, pw_cdr3_b_aa)}
#'   \item{barcodes}{Cell barcodes corresponding to matrix rows/columns}
#'   \item{tcr_data}{The formatted TCR data.frame used for analysis}
#'   If add_to_object=TRUE, returns the input object with distances stored.
#'
#' @export
#' @importFrom immApex getIR
#' @importFrom methods is
#'
#' @examples
#' \dontrun{
#'   # Calculate beta chain distances
#'   dist_results <- runTCRdist(seurat_obj, chains = "beta", organism = "human")
#'
#'   # Access distance matrix
#'   beta_distances <- dist_results$distances$pw_beta
#'
#'   # Calculate both chains
#'   dist_results <- runTCRdist(seurat_obj, chains = c("alpha", "beta"))
#'
#'   # Works with SingleCellExperiment too
#'   dist_results <- runTCRdist(sce, chains = "beta")
#' }
runTCRdist <- function(input,
                       chains = "beta",
                       organism = "human",
                       compute_distances = TRUE,
                       add_to_object = FALSE) {

  # Determine input type
  .is_sce <- methods::is(input, "SingleCellExperiment")
  .is_seurat <- methods::is(input, "Seurat")

  if (!.is_sce && !.is_seurat) {
    stop("Input must be a Seurat or SingleCellExperiment object")
  }

  message("Extracting TCR sequences from object...")

  # Map chain names
  chain_map <- list(
    "alpha" = "TRA",
    "beta" = "TRB"
  )

  # Extract data for each requested chain
  tcr_list <- list()

  for (chain in chains) {
    immapex_chain <- chain_map[[chain]]
    if (is.null(immapex_chain)) {
      stop("Invalid chain: ", chain, ". Must be 'alpha' or 'beta'.")
    }

    chain_data <- immApex::getIR(input, chains = immapex_chain)
    chain_data <- chain_data[!is.na(chain_data$cdr3_aa), ]

    if (nrow(chain_data) == 0) {
      warning(paste0("No valid ", chain, " chain sequences found."))
      next
    }

    tcr_list[[chain]] <- chain_data
  }

  if (length(tcr_list) == 0) {
    stop("No valid TCR sequences found for the specified chains.")
  }

  # Format data for tcrdist3
  # tcrdist3 expects specific column names
  message("Formatting data for tcrdist3...")

  formatted_df <- data.frame(count = integer(), stringsAsFactors = FALSE)

  # Add alpha chain data if present
  if ("alpha" %in% names(tcr_list)) {
    alpha_data <- tcr_list[["alpha"]]
    formatted_df$v_a_gene <- alpha_data$v
    formatted_df$j_a_gene <- alpha_data$j
    formatted_df$cdr3_a_aa <- alpha_data$cdr3_aa
    formatted_df$count <- 1
  }

  # Add beta chain data if present
  if ("beta" %in% names(tcr_list)) {
    beta_data <- tcr_list[["beta"]]
    if (nrow(formatted_df) == 0) {
      formatted_df <- data.frame(count = rep(1, nrow(beta_data)),
                                stringsAsFactors = FALSE)
    }
    formatted_df$v_b_gene <- beta_data$v
    formatted_df$j_b_gene <- beta_data$j
    formatted_df$cdr3_b_aa <- beta_data$cdr3_aa
  }

  # Add cell barcodes (not used by tcrdist3 but useful for mapping back)
  if ("beta" %in% names(tcr_list)) {
    barcodes <- tcr_list[["beta"]]$barcode
  } else {
    barcodes <- tcr_list[["alpha"]]$barcode
  }

  message("Calculating TCR distances for ", nrow(formatted_df), " sequences...")

  # Calculate distances
  dist_results <- calculate.tcrDist(
    df = formatted_df,
    organism = organism,
    chains = chains,
    compute_distances = compute_distances
  )

  # Prepare output
  output <- list(
    distances = dist_results,
    barcodes = barcodes,
    tcr_data = formatted_df
  )

  if (add_to_object) {
    if (.is_seurat) {
      input@misc$tcrdist <- output
      message("TCR distances added to object at obj@misc$tcrdist")
    } else {
      # For SCE, store in metadata
      S4Vectors::metadata(input)$tcrdist <- output
      message("TCR distances added to object at metadata(obj)$tcrdist")
    }
    return(input)
  } else {
    return(output)
  }
}
