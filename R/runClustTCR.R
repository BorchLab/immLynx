#' Run clusTCR Clustering on scRepertoire Data
#'
#' @description This function extracts TCR sequences from a Seurat or
#'   SingleCellExperiment object with scRepertoire data and performs clustering
#'   using the clusTCR algorithm.
#'
#' @param input A Seurat or SingleCellExperiment object containing scRepertoire
#'   TCR data in the metadata.
#' @param chains Character string specifying which chains to use: "TRA", "TRB", or "both".
#'   Default is "TRB".
#' @param method Clustering method passed to clusTCR. Default is "mcl"
#'   (Markov Clustering), which is accurate for typical repertoire datasets.
#' @param combine_chains Logical. If TRUE and chains="both", concatenates alpha and beta
#'   sequences with "_". Default is FALSE (clusters chains separately).
#' @param return_object Logical. If TRUE, adds cluster assignments back to the input object
#'   metadata. If FALSE, returns only the clustering results. Default is TRUE.
#' @param column_prefix Prefix for the new metadata column(s). Default is "clustcr".
#' @param ... Additional arguments passed to calculate.clustcr (e.g., inflation).
#'
#' @return If return_object=TRUE, returns the input object with cluster assignments added
#'   to metadata. If return_object=FALSE, returns a data.frame with barcodes and cluster
#'   assignments.
#'
#' @export
#' @importFrom immApex getIR
#' @importFrom methods is
#'
#' @examples
#' data(immLynx_example)
#' \donttest{
#'   # Cluster TRB chain using MCL algorithm
#'   seurat_obj <- runClustTCR(immLynx_example,
#'                             chains = "TRB")
#'
#'   # Adjust MCL inflation parameter
#'   seurat_obj <- runClustTCR(immLynx_example,
#'                             chains = "TRB",
#'                             inflation = 3.0)
#'
#'   # Cluster both chains separately
#'   seurat_obj <- runClustTCR(immLynx_example,
#'                             chains = "both")
#'
#'   # Combine alpha and beta chains before clustering
#'   seurat_obj <- runClustTCR(immLynx_example,
#'                             chains = "both",
#'                             combine_chains = TRUE)
#'
#'   # Get results as data.frame
#'   clusters_df <- runClustTCR(immLynx_example,
#'                              chains = "TRB",
#'                              return_object = FALSE)
#' }
runClustTCR <- function(input,
                        chains = c("TRB", "TRA", "both"),
                        method = "mcl",
                        combine_chains = FALSE,
                        return_object = TRUE,
                        column_prefix = "clustcr",
                       ...) {

  chains <- match.arg(chains)

  # Determine input type
.is_sce <- methods::is(input, "SingleCellExperiment")
  .is_seurat <- methods::is(input, "Seurat")

  if (!.is_sce && !.is_seurat) {
    stop("Input must be a Seurat or SingleCellExperiment object")
  }

  # Get cell names based on object type
  .get_cells <- function(obj) {
    if (methods::is(obj, "SingleCellExperiment")) {
      colnames(obj)
    } else {
      colnames(obj)
    }
  }

  # Add metadata based on object type
  .add_metadata <- function(obj, col_name, values, cell_names) {
    if (methods::is(obj, "SingleCellExperiment")) {
      col_vec <- rep(NA, ncol(obj))
      names(col_vec) <- colnames(obj)
      col_vec[cell_names] <- values
      SummarizedExperiment::colData(obj)[[col_name]] <- col_vec
    } else {
      col_vec <- rep(NA, ncol(obj))
      names(col_vec) <- colnames(obj)
      col_vec[cell_names] <- values
      obj[[col_name]] <- col_vec
    }
    obj
  }

  # Extract TCR data using immApex
  message("Extracting TCR sequences from object...")

  if (chains == "both") {
    tra_data <- immApex::getIR(input, chains = "TRA")
    trb_data <- immApex::getIR(input, chains = "TRB")

    if (combine_chains) {
      # Create combined sequences for cells with both chains
      message("Combining TRA and TRB sequences...")
      combined_data <- merge(tra_data, trb_data,
                            by = "barcode",
                            suffixes = c("_TRA", "_TRB"))

      # Filter cells with both chains
      combined_data <- combined_data[!is.na(combined_data$cdr3_aa_TRA) &
                                     !is.na(combined_data$cdr3_aa_TRB), ]

      if (nrow(combined_data) == 0) {
        stop("No cells found with both TRA and TRB chains.")
      }

      # Concatenate sequences
      sequences <- paste0(combined_data$cdr3_aa_TRA, "_", combined_data$cdr3_aa_TRB)
      barcodes <- combined_data$barcode

      # Cluster
      clusters <- calculate.clustcr(sequences = sequences, method = method, ...)

      # Create result data frame
      result_df <- data.frame(
        barcode = barcodes,
        combined_sequence = sequences,
        cluster = clusters,
        stringsAsFactors = FALSE
      )

      if (return_object) {
        input <- .add_metadata(input, paste0(column_prefix, "_combined"),
                              result_df$cluster, result_df$barcode)
        col_name <- paste0(column_prefix, "_combined")
        message("Cluster assignments added to metadata as '",
                col_name, "'")
        return(input)
      } else {
        return(result_df)
      }

    } else {
      # Cluster chains separately
      message("Clustering TRA and TRB separately...")
      results <- list()

      for (chain_name in c("TRA", "TRB")) {
        chain_data <- if (chain_name == "TRA") tra_data else trb_data
        chain_data <- chain_data[!is.na(chain_data$cdr3_aa), ]

        if (nrow(chain_data) == 0) {
          warning("No valid ", chain_name,
                  " sequences found. Skipping...")
          next
        }

        # Cluster
        clusters <- calculate.clustcr(
          sequences = chain_data$cdr3_aa,
          method = method, ...
        )

        result_df <- data.frame(
          barcode = chain_data$barcode,
          cdr3_aa = chain_data$cdr3_aa,
          chain = chain_name,
          cluster = clusters,
          stringsAsFactors = FALSE
        )

        results[[chain_name]] <- result_df

        if (return_object) {
          input <- .add_metadata(
            input, paste0(column_prefix, "_", chain_name),
            result_df$cluster, result_df$barcode
          )
        }
      }

      if (return_object) {
        message(
          "Cluster assignments added to metadata as '",
          column_prefix, "_TRA' and '",
          column_prefix, "_TRB'"
        )
        return(input)
      } else {
        return(do.call(rbind, results))
      }
    }

  } else {
    # Single chain analysis
    chain_data <- immApex::getIR(input, chains = chains)
    chain_data <- chain_data[!is.na(chain_data$cdr3_aa), ]

    if (nrow(chain_data) == 0) {
      stop("No valid ", chains, " sequences found.")
    }

    message("Found ", nrow(chain_data), " valid ",
            chains, " sequences")

    # Cluster
    clusters <- calculate.clustcr(
      sequences = chain_data$cdr3_aa,
      method = method, ...
    )

    result_df <- data.frame(
      barcode = chain_data$barcode,
      cdr3_aa = chain_data$cdr3_aa,
      chain = chains,
      cluster = clusters,
      stringsAsFactors = FALSE
    )

    if (return_object) {
      col_name <- paste0(column_prefix, "_", chains)
      input <- .add_metadata(input, col_name,
                            result_df$cluster, result_df$barcode)
      message("Cluster assignments added to metadata as '",
              col_name, "'")
      return(input)
    } else {
      return(result_df)
    }
  }
}
