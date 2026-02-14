#' Run GLIPH2 Analysis for TCR Specificity Groups
#'
#' @description Identifies TCR specificity groups using the GLIPH (Grouping of
#'   Lymphocyte Interactions by Paratope Hotspots) algorithm. This function
#'   provides an interface to the turboGliph R package which implements GLIPH2.
#'
#' @param input A Seurat object with scRepertoire data, or a data.frame with TCR data.
#' @param chains Which chain to analyze. Currently only "TRB"/"beta" is supported
#'   by GLIPH. Default is "TRB".
#' @param local_similarities Logical. If TRUE, performs local similarity analysis
#'   (motif-based clustering). Default is TRUE.
#' @param global_similarities Logical. If TRUE, performs global similarity analysis
#'   (Hamming distance-based). Default is TRUE.
#' @param local_method Method for local enrichment testing: "fisher" (recommended)
#'   or "rrs" (repeated random sampling). Default is "fisher".
#' @param global_method Method for global enrichment: "fisher" or "cutoff".
#'   Default is "fisher".
#' @param motif_length Length of CDR3 motifs to analyze. Default is 3.
#' @param gccutoff Global cluster distance cutoff (for cutoff method). Default is 1.
#' @param vgene_match Logical. If TRUE, requires V gene matching in clusters.
#'   Default is FALSE.
#' @param n_cores Number of cores for parallel processing. Default is 1.
#' @param return_seurat Logical. If TRUE and input is Seurat object, adds cluster
#'   assignments to metadata. Default is TRUE.
#' @param column_prefix Prefix for metadata columns. Default is "gliph".
#'
#' @return If return_seurat=TRUE and input is Seurat, returns the object with
#'   GLIPH cluster assignments. Otherwise returns a list containing:
#'   \describe{
#'     \item{clusters}{Data.frame with cluster membership and properties}
#'     \item{motifs}{Data.frame with enriched motifs}
#'     \item{network}{Edge list for cluster network}
#'   }
#'
#' @details
#' GLIPH identifies TCR specificity groups based on:
#' \itemize{
#'   \item Local similarity: Shared CDR3 motifs that are enriched vs. reference
#'   \item Global similarity: CDR3 sequences within a Hamming distance threshold
#' }
#'
#' The turboGliph implementation is significantly faster than the original Perl
#' version (50-500x depending on dataset size).
#'
#' @section Requirements:
#' This function requires the turboGliph package to be installed:
#' \code{remotes::install_github("HetzDra/turboGliph")}
#'
#' @references
#' Glanville et al. (2017). Identifying specificity groups in the T cell receptor
#' repertoire. Nature, 547(7661), 94-98.
#'
#' @seealso \code{\link{runMetaclonotypist}} for an alternative metaclone discovery method
#'
#' @export
#' @importFrom immApex getIR
#' @importFrom methods is
#'
#' @examples
#' \dontrun{
#'   # Run GLIPH2 analysis
#'   seurat_obj <- runGLIPH(seurat_obj, chains = "TRB")
#'
#'   # Run with V gene matching
#'   seurat_obj <- runGLIPH(seurat_obj, vgene_match = TRUE)
#'
#'   # Get full results
#'   gliph_results <- runGLIPH(seurat_obj, return_seurat = FALSE)
#'   head(gliph_results$clusters)
#' }
runGLIPH <- function(input,
                     chains = c("TRB", "beta"),
                     local_similarities = TRUE,
                     global_similarities = TRUE,
                     local_method = c("fisher", "rrs"),
                     global_method = c("fisher", "cutoff"),
                     motif_length = 3,
                     gccutoff = 1,
                     vgene_match = FALSE,
                     n_cores = 1,
                     return_seurat = TRUE,
                     column_prefix = "gliph") {

  chains <- match.arg(chains)
  local_method <- match.arg(local_method)
  global_method <- match.arg(global_method)

  # Check for turboGliph
  if (!requireNamespace("turboGliph", quietly = TRUE)) {
    stop("Package 'turboGliph' is required for GLIPH analysis.\n",
         "Install it with: remotes::install_github('HetzDra/turboGliph')")
  }

  # Extract TCR data
  is_seurat <- methods::is(input, "Seurat") || methods::is(input, "SingleCellExperiment")

  if (is_seurat) {
    tcr_data <- immApex::getIR(input, chains = "TRB")
    tcr_data <- tcr_data[!is.na(tcr_data$cdr3_aa), ]
  } else {
    tcr_data <- input
    if (!"cdr3_aa" %in% names(tcr_data) && !"CDR3b" %in% names(tcr_data)) {
      stop("Input must contain 'cdr3_aa' or 'CDR3b' column")
    }
  }

  if (nrow(tcr_data) == 0) {
    stop("No valid TCR sequences found")
  }

  message("Running GLIPH2 analysis on ", nrow(tcr_data), " sequences...")

  # Prepare data for turboGliph
  # turboGliph expects column named 'CDR3b' and optionally 'TRBV'
  gliph_input <- data.frame(
    CDR3b = tcr_data$cdr3_aa,
    stringsAsFactors = FALSE
  )

  # Add V gene if available and required
  if (vgene_match && "v" %in% names(tcr_data)) {
    gliph_input$TRBV <- tcr_data$v
  }

  # Add barcode for mapping back
  if ("barcode" %in% names(tcr_data)) {
    gliph_input$barcode <- tcr_data$barcode
  }

  # Run GLIPH2 via turboGliph
  gliph_result <- turboGliph::gliph_combined(
    cdr3_sequences = gliph_input,
    local_similarities = local_similarities,
    local_method = local_method,
    global_similarities = global_similarities,
    global_method = global_method,
    motif_length = motif_length,
    gccutoff = gccutoff,
    vgene_match = vgene_match,
    n_cores = n_cores,
    result_folder = NULL  # Don't save files
  )

  # Extract results
  # turboGliph returns a list with multiple components
  clusters <- gliph_result$cluster_properties
  network <- gliph_result$network

  # Create summary of results
  n_clusters <- if (!is.null(clusters)) nrow(clusters) else 0
  message("Identified ", n_clusters, " GLIPH specificity groups")

  if (return_seurat && is_seurat && "barcode" %in% names(tcr_data)) {
    # Map cluster assignments back to cells
    cluster_col <- rep(NA_character_, ncol(input))
    names(cluster_col) <- colnames(input)

    # Process cluster membership
    if (!is.null(gliph_result$cluster_membership)) {
      membership <- gliph_result$cluster_membership

      # Map via barcode
      for (i in seq_len(nrow(membership))) {
        cdr3 <- membership$CDR3b[i]
        cluster_id <- membership$cluster_id[i]

        # Find barcodes with this CDR3
        matching_idx <- which(gliph_input$CDR3b == cdr3)
        if (length(matching_idx) > 0) {
          matching_barcodes <- gliph_input$barcode[matching_idx]
          cluster_col[matching_barcodes] <- as.character(cluster_id)
        }
      }
    }

    input[[paste0(column_prefix, "_cluster")]] <- cluster_col

    message("GLIPH cluster assignments added to metadata as '",
            paste0(column_prefix, "_cluster"), "'")

    return(input)

  } else {
    # Return full results
    return(list(
      clusters = clusters,
      membership = gliph_result$cluster_membership,
      motifs = gliph_result$motif_enrichment,
      global = gliph_result$global_enrichment,
      network = network,
      parameters = list(
        local_similarities = local_similarities,
        global_similarities = global_similarities,
        local_method = local_method,
        global_method = global_method,
        motif_length = motif_length,
        vgene_match = vgene_match
      )
    ))
  }
}


#' Score GLIPH Clusters
#'
#' @description Calculates quality scores for GLIPH specificity groups based on
#'   cluster properties like size, motif enrichment, and V gene usage.
#'
#' @param gliph_results Output from runGLIPH() with return_seurat=FALSE
#' @param min_size Minimum cluster size to consider. Default is 2.
#' @param weight_size Weight for cluster size in scoring. Default is 0.3.
#' @param weight_motif Weight for motif enrichment in scoring. Default is 0.4.
#' @param weight_vgene Weight for V gene consistency in scoring. Default is 0.3.
#'
#' @return A data.frame with cluster scores and quality metrics
#'
#' @export
#'
#' @examples
#' \dontrun{
#'   gliph_results <- runGLIPH(seurat_obj, return_seurat = FALSE)
#'   cluster_scores <- scoreGLIPH(gliph_results)
#' }
scoreGLIPH <- function(gliph_results,
                       min_size = 2,
                       weight_size = 0.3,
                       weight_motif = 0.4,
                       weight_vgene = 0.3) {

  if (is.null(gliph_results$clusters)) {
    stop("No clusters found in gliph_results")
  }

  clusters <- gliph_results$clusters

  # Filter by size
  clusters <- clusters[clusters$size >= min_size, ]

  if (nrow(clusters) == 0) {
    warning("No clusters meet minimum size requirement")
    return(NULL)
  }

  # Calculate component scores
  # Size score: log-scaled
  max_size <- max(clusters$size, na.rm = TRUE)
  clusters$size_score <- log10(clusters$size + 1) / log10(max_size + 1)

  # Motif score: based on enrichment p-value if available
  if ("motif_pvalue" %in% names(clusters)) {
    clusters$motif_score <- -log10(clusters$motif_pvalue + 1e-100)
    max_motif <- max(clusters$motif_score, na.rm = TRUE)
    clusters$motif_score <- clusters$motif_score / max_motif
  } else {
    clusters$motif_score <- 1  # Default if no p-value
  }

  # V gene score: proportion of dominant V gene
  if ("vgene_fraction" %in% names(clusters)) {
    clusters$vgene_score <- clusters$vgene_fraction
  } else {
    clusters$vgene_score <- 1  # Default if no V gene info
  }

  # Combined score
  clusters$quality_score <- (
    weight_size * clusters$size_score +
    weight_motif * clusters$motif_score +
    weight_vgene * clusters$vgene_score
  )

  # Normalize to 0-1
  clusters$quality_score <- clusters$quality_score / max(clusters$quality_score, na.rm = TRUE)

  # Sort by score
  clusters <- clusters[order(clusters$quality_score, decreasing = TRUE), ]
  rownames(clusters) <- NULL

  return(clusters)
}


#' Extract GLIPH Motifs
#'
#' @description Extracts enriched motifs from GLIPH analysis results with
#'   associated statistics.
#'
#' @param gliph_results Output from runGLIPH() with return_seurat=FALSE
#' @param fdr_threshold FDR threshold for significant motifs. Default is 0.05.
#' @param min_observations Minimum number of observations for a motif. Default is 2.
#'
#' @return A data.frame with enriched motifs:
#'   \describe{
#'     \item{motif}{The CDR3 motif sequence}
#'     \item{count}{Number of occurrences}
#'     \item{expected}{Expected count under null}
#'     \item{fold_enrichment}{Observed/Expected ratio}
#'     \item{pvalue}{Enrichment p-value}
#'     \item{fdr}{FDR-adjusted p-value}
#'   }
#'
#' @export
#'
#' @examples
#' \dontrun{
#'   gliph_results <- runGLIPH(seurat_obj, return_seurat = FALSE)
#'   motifs <- extractGLIPHmotifs(gliph_results)
#' }
extractGLIPHmotifs <- function(gliph_results,
                               fdr_threshold = 0.05,
                               min_observations = 2) {

  if (is.null(gliph_results$motifs)) {
    warning("No motif data found in gliph_results")
    return(NULL)
  }

  motifs <- gliph_results$motifs

  # Filter by observations
  motifs <- motifs[motifs$count >= min_observations, ]

  # Add FDR if not present
  if (!"fdr" %in% names(motifs) && "pvalue" %in% names(motifs)) {
    motifs$fdr <- stats::p.adjust(motifs$pvalue, method = "BH")
  }

  # Filter by FDR
  if ("fdr" %in% names(motifs)) {
    motifs <- motifs[motifs$fdr < fdr_threshold, ]
  }

  # Sort by enrichment
  if ("fold_enrichment" %in% names(motifs)) {
    motifs <- motifs[order(motifs$fold_enrichment, decreasing = TRUE), ]
  } else if ("pvalue" %in% names(motifs)) {
    motifs <- motifs[order(motifs$pvalue), ]
  }

  rownames(motifs) <- NULL

  return(motifs)
}
