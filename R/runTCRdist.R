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
#' @param max_sequences Maximum number of sequences to analyze. If NULL (default), uses all.
#'   Recommended to set this for large datasets (>5000 sequences) to avoid memory issues.
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
#'   # Subsample for large datasets
#'   dist_results <- runTCRdist(seurat_obj, chains = "beta", max_sequences = 5000)
#' }
runTCRdist <- function(input,
                       chains = "beta",
                       organism = "human",
                       compute_distances = TRUE,
                       max_sequences = NULL,
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
  barcodes <- character()

  # Helper function to format gene names for tcrdist3

  # tcrdist3 expects format like "TRBV5-1*01" with allele info
  format_gene_for_tcrdist <- function(genes) {
    # Handle NA
    genes[is.na(genes)] <- ""

    # Add *01 allele if not present (tcrdist3 requires allele info)
    needs_allele <- nchar(genes) > 0 & !grepl("\\*", genes)
    genes[needs_allele] <- paste0(genes[needs_allele], "*01")

    genes
  }

  # Add beta chain data if present
  if ("beta" %in% names(tcr_list)) {
    beta_data <- tcr_list[["beta"]]

    # Format V/J gene names for tcrdist3
    clean_v <- format_gene_for_tcrdist(beta_data$v)
    clean_j <- format_gene_for_tcrdist(beta_data$j)

    formatted_df <- data.frame(
      count = rep(1L, nrow(beta_data)),
      v_b_gene = clean_v,
      j_b_gene = clean_j,
      cdr3_b_aa = beta_data$cdr3_aa,
      stringsAsFactors = FALSE
    )
    barcodes <- beta_data$barcode
  }

  # Add alpha chain data if present
  if ("alpha" %in% names(tcr_list)) {
    alpha_data <- tcr_list[["alpha"]]

    clean_v <- format_gene_for_tcrdist(alpha_data$v)
    clean_j <- format_gene_for_tcrdist(alpha_data$j)

    if (nrow(formatted_df) == 0) {
      formatted_df <- data.frame(
        count = rep(1L, nrow(alpha_data)),
        v_a_gene = clean_v,
        j_a_gene = clean_j,
        cdr3_a_aa = alpha_data$cdr3_aa,
        stringsAsFactors = FALSE
      )
      barcodes <- alpha_data$barcode
    } else {
      # Add alpha columns to existing beta data
      formatted_df$v_a_gene <- clean_v
      formatted_df$j_a_gene <- clean_j
      formatted_df$cdr3_a_aa <- alpha_data$cdr3_aa
    }
  }


  # Remove rows with empty CDR3 sequences or missing V/J genes
  # tcrdist3 requires valid gene names - empty strings cause errors
  if ("beta" %in% chains) {
    valid_cdr3 <- !is.na(formatted_df$cdr3_b_aa) & nchar(formatted_df$cdr3_b_aa) > 0
    valid_v <- !is.na(formatted_df$v_b_gene) & nchar(formatted_df$v_b_gene) > 0
    valid_j <- !is.na(formatted_df$j_b_gene) & nchar(formatted_df$j_b_gene) > 0
    valid_rows <- valid_cdr3 & valid_v & valid_j
  } else {
    valid_cdr3 <- !is.na(formatted_df$cdr3_a_aa) & nchar(formatted_df$cdr3_a_aa) > 0
    valid_v <- !is.na(formatted_df$v_a_gene) & nchar(formatted_df$v_a_gene) > 0
    valid_j <- !is.na(formatted_df$j_a_gene) & nchar(formatted_df$j_a_gene) > 0
    valid_rows <- valid_cdr3 & valid_v & valid_j
  }

  n_removed <- sum(!valid_rows)
  if (n_removed > 0) {
    message("Removing ", n_removed, " sequences with missing CDR3/V/J gene info")
  }

  formatted_df <- formatted_df[valid_rows, , drop = FALSE]
  barcodes <- barcodes[valid_rows]

  if (nrow(formatted_df) == 0) {
    stop("No valid TCR sequences remaining after filtering. ",
         "Ensure sequences have CDR3, V gene, and J gene annotations.")
  }

  # Subsample if requested or if dataset is very large
  if (!is.null(max_sequences) && nrow(formatted_df) > max_sequences) {
    message("Subsampling from ", nrow(formatted_df), " to ", max_sequences, " sequences...")
    set.seed(42)
    idx <- sample(nrow(formatted_df), max_sequences)
    formatted_df <- formatted_df[idx, , drop = FALSE]
    barcodes <- barcodes[idx]
  } else if (nrow(formatted_df) > 10000 && is.null(max_sequences)) {
    warning("Large dataset (", nrow(formatted_df), " sequences). ",
            "Consider using max_sequences parameter to subsample. ",
            "Distance matrix will be ", nrow(formatted_df), "x", nrow(formatted_df), " elements.")
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
