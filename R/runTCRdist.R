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
#' data(immLynx_example)
#' \donttest{
#'   # Calculate TCR distances for beta chain
#'   dist_results <- runTCRdist(immLynx_example,
#'                              chains = "beta")
#'
#'   # Access the distance matrix
#'   dist_matrix <- dist_results$distances
#'   barcodes <- dist_results$barcodes
#'
#'   # Calculate for both chains
#'   dist_both <- runTCRdist(immLynx_example,
#'                           chains = c("alpha", "beta"))
#'
#'   # Add distances directly to the Seurat object
#'   seurat_obj <- runTCRdist(immLynx_example,
#'                            chains = "beta",
#'                            add_to_object = TRUE)
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
      warning("No valid ", chain, " chain sequences found.")
      next
    }

    tcr_list[[chain]] <- chain_data
  }

  if (length(tcr_list) == 0) {
    stop("No valid TCR sequences found for the specified chains.")
  }

  # Format data for tcrdist3
  # tcrdist3 expects specific column names and IMGT allele format
  message("Formatting data for tcrdist3...")

  # tcrdist3 requires allele suffixes (e.g., TRBV10-3*01)
  # scRepertoire data typically lacks allele info
  .add_allele <- function(genes) {
    ifelse(!is.na(genes) & !grepl("\\*", genes),
           paste0(genes, "*01"),
           genes)
  }

  if (length(tcr_list) == 2) {
    # Both chains: merge by barcode to ensure matched rows
    alpha_data <- tcr_list[["alpha"]]
    beta_data <- tcr_list[["beta"]]
    common_barcodes <- intersect(alpha_data$barcode, beta_data$barcode)

    if (length(common_barcodes) == 0) {
      stop("No cells have both alpha and beta chain sequences.")
    }

    alpha_data <- alpha_data[match(common_barcodes, alpha_data$barcode), ]
    beta_data <- beta_data[match(common_barcodes, beta_data$barcode), ]

    formatted_df <- data.frame(
      count = rep(1L, length(common_barcodes)),
      v_a_gene = .add_allele(alpha_data$v),
      j_a_gene = .add_allele(alpha_data$j),
      cdr3_a_aa = alpha_data$cdr3_aa,
      v_b_gene = .add_allele(beta_data$v),
      j_b_gene = .add_allele(beta_data$j),
      cdr3_b_aa = beta_data$cdr3_aa,
      stringsAsFactors = FALSE
    )
    barcodes <- common_barcodes
  } else {
    # Single chain
    chain_name <- names(tcr_list)[1]
    chain_data <- tcr_list[[chain_name]]
    barcodes <- chain_data$barcode

    if (chain_name == "alpha") {
      formatted_df <- data.frame(
        count = rep(1L, nrow(chain_data)),
        v_a_gene = .add_allele(chain_data$v),
        j_a_gene = .add_allele(chain_data$j),
        cdr3_a_aa = chain_data$cdr3_aa,
        stringsAsFactors = FALSE
      )
    } else {
      formatted_df <- data.frame(
        count = rep(1L, nrow(chain_data)),
        v_b_gene = .add_allele(chain_data$v),
        j_b_gene = .add_allele(chain_data$j),
        cdr3_b_aa = chain_data$cdr3_aa,
        stringsAsFactors = FALSE
      )
    }
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
      Seurat::Misc(input, slot = "tcrdist") <- output
      message("TCR distances added to object misc slot")
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
