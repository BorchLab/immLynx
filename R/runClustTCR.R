#' Run clusTCR Clustering on scRepertoire Data
#'
#' @description This function extracts TCR sequences from a Seurat object with
#'   scRepertoire data and performs clustering using the clusTCR algorithm.
#'
#' @param seurat_obj A Seurat object containing scRepertoire TCR data in the metadata.
#' @param chains Character string specifying which chains to use: "TRA", "TRB", or "both".
#'   Default is "TRB".
#' @param method Clustering method to use: 'mcl' (default) or 'dbscan'.
#' @param combine_chains Logical. If TRUE and chains="both", concatenates alpha and beta
#'   sequences with "_". Default is FALSE (clusters chains separately).
#' @param return_seurat Logical. If TRUE, adds cluster assignments back to the Seurat object
#'   metadata. If FALSE, returns only the clustering results. Default is TRUE.
#' @param column_prefix Prefix for the new metadata column(s). Default is "clustcr".
#' @param ... Additional arguments passed to calculate.clustcr (e.g., inflation, eps, min_samples).
#'
#' @return If return_seurat=TRUE, returns the Seurat object with cluster assignments added
#'   to metadata. If return_seurat=FALSE, returns a data.frame with barcodes and cluster
#'   assignments.
#'
#' @export
#' @importFrom immApex getIR
#'
#' @examples
#'   # Cluster TRB sequences using MCL
#'   seurat_obj <- runClustCR(seurat_obj, chains = "TRB", method = "mcl", inflation = 2.5)
#'
#'   # Cluster both chains separately using DBSCAN
#'   seurat_obj <- runClustCR(seurat_obj, chains = "both", method = "dbscan", eps = 0.6)
#'
#'   # Cluster concatenated alpha-beta pairs
#'   seurat_obj <- runClustCR(seurat_obj, chains = "both", combine_chains = TRUE)
runClustTCR <- function(seurat_obj,
                        chains = c("TRB", "TRA", "both"),
                        method = "mcl",
                        combine_chains = FALSE,
                        return_seurat = TRUE,
                        column_prefix = "clustcr",
                       ...) {
  
  chains <- match.arg(chains)
  
  # Extract TCR data using immApex
  message("Extracting TCR sequences from Seurat object...")
  
  if (chains == "both") {
    tra_data <- immApex::getIR(seurat_obj, chains = "TRA")
    trb_data <- immApex::getIR(seurat_obj, chains = "TRB")
    
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
      
      if (return_seurat) {
        # Add to metadata
        cluster_col <- rep(NA, ncol(seurat_obj))
        names(cluster_col) <- colnames(seurat_obj)
        cluster_col[result_df$barcode] <- result_df$cluster
        seurat_obj[[paste0(column_prefix, "_combined")]] <- cluster_col
        message("Cluster assignments added to metadata as '", paste0(column_prefix, "_combined"), "'")
        return(seurat_obj)
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
          warning(paste0("No valid ", chain_name, " sequences found. Skipping..."))
          next
        }
        
        # Cluster
        clusters <- calculate.clustcr(sequences = chain_data$cdr3_aa, method = method, ...)
        
        result_df <- data.frame(
          barcode = chain_data$barcode,
          cdr3_aa = chain_data$cdr3_aa,
          chain = chain_name,
          cluster = clusters,
          stringsAsFactors = FALSE
        )
        
        results[[chain_name]] <- result_df
        
        if (return_seurat) {
          cluster_col <- rep(NA, ncol(seurat_obj))
          names(cluster_col) <- colnames(seurat_obj)
          cluster_col[result_df$barcode] <- result_df$cluster
          seurat_obj[[paste0(column_prefix, "_", chain_name)]] <- cluster_col
        }
      }
      
      if (return_seurat) {
        message("Cluster assignments added to metadata as '", 
                paste0(column_prefix, "_TRA"), "' and '", 
                paste0(column_prefix, "_TRB"), "'")
        return(seurat_obj)
      } else {
        return(do.call(rbind, results))
      }
    }
    
  } else {
    # Single chain analysis
    chain_data <- immApex::getIR(seurat_obj, chains = chains)
    chain_data <- chain_data[!is.na(chain_data$cdr3_aa), ]
    
    if (nrow(chain_data) == 0) {
      stop(paste0("No valid ", chains, " sequences found."))
    }
    
    message(paste0("Found ", nrow(chain_data), " valid ", chains, " sequences"))
    
    # Cluster
    clusters <- calculate.clustcr(sequences = chain_data$cdr3_aa, method = method, ...)
    
    result_df <- data.frame(
      barcode = chain_data$barcode,
      cdr3_aa = chain_data$cdr3_aa,
      chain = chains,
      cluster = clusters,
      stringsAsFactors = FALSE
    )
    
    if (return_seurat) {
      cluster_col <- rep(NA, ncol(seurat_obj))
      names(cluster_col) <- colnames(seurat_obj)
      cluster_col[result_df$barcode] <- result_df$cluster
      seurat_obj[[paste0(column_prefix, "_", chains)]] <- cluster_col
      message("Cluster assignments added to metadata as '", 
              paste0(column_prefix, "_", chains), "'")
      return(seurat_obj)
    } else {
      return(result_df)
    }
  }
}
