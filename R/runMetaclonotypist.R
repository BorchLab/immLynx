#' Run Metaclonotypist for TCR Metaclone Discovery
#'
#' @description Identifies TCR metaclones (groups of related T cell receptors)
#'   using the metaclonotypist pipeline. Metaclonotypist uses a two-stage approach:
#'   fast edit-distance-based screening followed by TCRdist or SCEPTR refinement.
#'
#' @param input A Seurat object with scRepertoire data, or a data.frame with TCR data.
#' @param chains Which chain to analyze: "alpha" or "beta". Default is "beta".
#' @param method Distance metric for refinement: "tcrdist" (default) or "sceptr".
#' @param max_edits Maximum CDR3 edit distance for initial screening. Default is 2.
#' @param max_dist Maximum distance threshold for clustering. Default is 20 for
#'   tcrdist, 1.5 for sceptr.
#' @param clustering Clustering algorithm: "cc" (connected components, default),
#'   "leiden", "louvain", or "mcl".
#' @param resolution Resolution parameter for leiden/louvain clustering. Default is 1.0.
#' @param return_seurat Logical. If TRUE and input is Seurat object, adds metaclone
#'   assignments to metadata. Default is TRUE.
#' @param column_name Name for the metadata column. Default is "metaclone".
#'
#' @return If return_seurat=TRUE and input is Seurat, returns the object with
#'   metaclone assignments added to metadata. Otherwise returns a data.frame with:
#'   \describe{
#'     \item{barcode}{Cell barcode}
#'     \item{cdr3_aa}{CDR3 amino acid sequence}
#'     \item{metaclone}{Metaclone cluster assignment}
#'     \item{metaclone_size}{Number of cells in the metaclone}
#'   }
#'
#' @details
#' Metaclonotypist identifies groups of related TCRs that may recognize similar

#' antigens. The algorithm:
#' 1. Uses the Symdel algorithm for fast edit-distance-based candidate identification
#' 2. Refines candidates using TCRdist or SCEPTR similarity metrics
#' 3. Applies graph-based clustering (Leiden by default) to identify metaclones
#'
#' @references
#' Metaclonotypist: \url{https://github.com/qimmuno/metaclonotypist}
#'
#' @export
#' @importFrom immApex getIR
#' @importFrom methods is
#'
#' @examples
#' data(immLynx_example)
#' \dontrun{
#'   # Run metaclonotypist on beta chain
#'   seurat_obj <- runMetaclonotypist(immLynx_example, chains = "beta")
#'
#'   # Get results as data.frame instead of adding to object
#'   metaclones <- runMetaclonotypist(immLynx_example,
#'                                    return_seurat = FALSE)
#'
#'   # Adjust edit distance threshold
#'   seurat_obj <- runMetaclonotypist(immLynx_example,
#'                                    max_edits = 3,
#'                                    max_dist = 50)
#' }
runMetaclonotypist <- function(input,
                               chains = c("beta", "alpha"),
                               method = c("tcrdist", "sceptr"),
                               max_edits = 2,
                               max_dist = NULL,
                               clustering = c("cc", "leiden", "louvain", "mcl"),
                               resolution = 1.0,
                               return_seurat = TRUE,
                               column_name = "metaclone") {

  chains <- match.arg(chains)
  method <- match.arg(method)
  clustering <- match.arg(clustering)

  # Set default max_dist based on method

if (is.null(max_dist)) {
    max_dist <- if (method == "tcrdist") 20 else 1.5
  }

  # Extract TCR data
  is_seurat <- methods::is(input, "Seurat") || methods::is(input, "SingleCellExperiment")

  if (is_seurat) {
    chain_code <- if (chains == "alpha") "TRA" else "TRB"
    tcr_data <- immApex::getIR(input, chains = chain_code)
    tcr_data <- tcr_data[!is.na(tcr_data$cdr3_aa), ]
  } else {
    tcr_data <- input
    if (!all(c("barcode", "cdr3_aa") %in% names(tcr_data))) {
      stop("Input data.frame must contain 'barcode' and 'cdr3_aa' columns")
    }
    tcr_data <- tcr_data[!is.na(tcr_data$cdr3_aa), ]
  }

  if (nrow(tcr_data) == 0) {
    stop("No valid TCR sequences found")
  }

  message("Running metaclonotypist on ", nrow(tcr_data), " ", chains, " chain sequences...")
  message("Method: ", method, ", Max edits: ", max_edits, ", Max distance: ", max_dist)

  # Run metaclonotypist
  results <- calculate.metaclonotypist(
    cdr3_sequences = tcr_data$cdr3_aa,
    v_genes = if ("v" %in% names(tcr_data)) tcr_data$v else NULL,
    j_genes = if ("j" %in% names(tcr_data)) tcr_data$j else NULL,
    chain = chains,
    method = method,
    max_edits = max_edits,
    max_dist = max_dist,
    clustering = clustering,
    resolution = resolution
  )

  # Map cluster results back to input sequences
  # metaclonotypist only returns clustered sequences (singletons dropped)
  # Build bioidentity for mapping (same as in calculate.metaclonotypist)
  if ("v" %in% names(tcr_data) && !is.null(tcr_data$v)) {
    input_bioidentity <- paste0(tcr_data$v, ":", tcr_data$cdr3_aa)
  } else {
    input_bioidentity <- tcr_data$cdr3_aa
  }

  cluster_map <- stats::setNames(
    as.character(results$cluster),
    as.character(results$node)
  )

  mapped_clusters <- unname(cluster_map[input_bioidentity])

  result_df <- data.frame(
    barcode = tcr_data$barcode,
    cdr3_aa = tcr_data$cdr3_aa,
    metaclone = mapped_clusters,
    stringsAsFactors = FALSE
  )

  # Add metaclone size (only for clustered cells)
  non_na <- result_df$metaclone[!is.na(result_df$metaclone)]
  metaclone_sizes <- table(non_na)
  result_df$metaclone_size <- as.integer(
    metaclone_sizes[as.character(result_df$metaclone)]
  )

  n_clustered <- sum(!is.na(mapped_clusters))
  n_metaclones <- length(unique(non_na))
  message("Identified ", n_metaclones, " metaclones covering ", n_clustered,
          " of ", nrow(tcr_data), " sequences")

  if (return_seurat && is_seurat) {
    # Add to Seurat metadata
    metaclone_col <- rep(NA_character_, ncol(input))
    names(metaclone_col) <- colnames(input)
    metaclone_col[result_df$barcode] <- as.character(result_df$metaclone)

    size_col <- rep(NA_integer_, ncol(input))
    names(size_col) <- colnames(input)
    size_col[result_df$barcode] <- result_df$metaclone_size

    input[[column_name]] <- metaclone_col
    input[[paste0(column_name, "_size")]] <- size_col

    message("Metaclone assignments added to metadata as '", column_name, "'")
    return(input)
  } else {
    return(result_df)
  }
}


#' Internal Metaclonotypist Calculation Function
#' @description Calls metaclonotypist Python package via basilisk
#' @param cdr3_sequences Character vector of CDR3 amino acid sequences
#' @param v_genes Character vector of V gene names (optional)
#' @param j_genes Character vector of J gene names (optional)
#' @param chain Chain type: "alpha" or "beta"
#' @param method Distance method: "tcrdist" or "sceptr"
#' @param max_edits Maximum edit distance for screening
#' @param max_dist Maximum distance for clustering
#' @param clustering Clustering algorithm
#' @param resolution Clustering resolution
#' @return List with cluster assignments
#' @keywords internal
calculate.metaclonotypist <- function(cdr3_sequences,
                                      v_genes = NULL,
                                      j_genes = NULL,
                                      chain = "beta",
                                      method = "tcrdist",
                                      max_edits = 2,
                                      max_dist = 20,
                                      clustering = "leiden",
                                      resolution = 1.0) {

  proc <- basilisk::basiliskStart(immLynxEnv)
  on.exit(basilisk::basiliskStop(proc))

  results <- basilisk::basiliskRun(proc, function(cdr3_sequences, v_genes, j_genes, chain,
                                                   method, max_edits, max_dist,
                                                   clustering, resolution) {
    pd <- reticulate::import("pandas")
    mc <- reticulate::import("metaclonotypist")

    # Metaclonotypist uses pyrepseq standard column names:
    # CDR3B/CDR3A for CDR3 sequences, TRBV/TRAV for V genes, TRBJ/TRAJ for J genes
    chain_letter <- toupper(substr(chain, 1, 1))
    cdr3_col <- paste0("CDR3", chain_letter)
    v_col <- paste0("TR", chain_letter, "V")
    j_col <- paste0("TR", chain_letter, "J")

    df_data <- list()
    df_data[[cdr3_col]] <- cdr3_sequences

    if (!is.null(v_genes)) {
      df_data[[v_col]] <- v_genes
    }
    if (!is.null(j_genes)) {
      df_data[[j_col]] <- j_genes
    }

    # Build bioidentity column (required for node labeling)
    # bioidentity is typically V gene + CDR3 combination
    if (!is.null(v_genes)) {
      df_data[["bioidentity"]] <- paste0(v_genes, ":", cdr3_sequences)
    } else {
      df_data[["bioidentity"]] <- cdr3_sequences
    }

    df <- pd$DataFrame(df_data)

    # Run metaclonotypist based on method
    if (method == "tcrdist") {
      result <- mc$metaclonotypist(
        df = df,
        chain = chain,
        max_edits = as.integer(max_edits),
        max_tcrdist = as.numeric(max_dist),
        clustering = clustering
      )
    } else {
      result <- mc$metaclonotypist_sceptr(
        df = df,
        chain = chain,
        max_edits = as.integer(max_edits),
        max_sceptrdist = as.numeric(max_dist),
        clustering = clustering
      )
    }

    # Convert result to R data.frame
    result_r <- reticulate::py_to_r(result)

    list(
      cluster = result_r$cluster,
      node = result_r$node
    )

  }, cdr3_sequences = cdr3_sequences, v_genes = v_genes, j_genes = j_genes,
     chain = chain, method = method, max_edits = max_edits, max_dist = max_dist,
     clustering = clustering, resolution = resolution)

  return(results)
}


#' Perform HLA Association Analysis on Metaclones
#'
#' @description Tests associations between metaclone membership and HLA alleles
#'   using Fisher's exact test with FDR correction.
#'
#' @param metaclone_data A data.frame with metaclone assignments, typically from
#'   runMetaclonotypist() with return_seurat=FALSE.
#' @param hla_data A data.frame with HLA typing information. Must have a 'barcode'
#'   or 'sample' column for matching and columns for each HLA allele.
#' @param by Column name to use for matching between datasets. Default is "barcode".
#' @param fdr_threshold FDR threshold for significance. Default is 0.05.
#'
#' @return A data.frame with HLA association results:
#'   \describe{
#'     \item{metaclone}{Metaclone identifier}
#'     \item{hla_allele}{HLA allele tested}
#'     \item{odds_ratio}{Association odds ratio}
#'     \item{pvalue}{Raw p-value from Fisher's exact test}
#'     \item{fdr}{FDR-adjusted p-value}
#'     \item{significant}{Whether FDR < threshold}
#'   }
#'
#' @export
#'
#' @examples
#' # Create example metaclone and HLA data
#' metaclone_data <- data.frame(
#'   barcode = paste0("cell_", 1:20),
#'   metaclone = rep(c("MC1", "MC2"), each = 10),
#'   stringsAsFactors = FALSE
#' )
#' hla_data <- data.frame(
#'   barcode = paste0("cell_", 1:20),
#'   HLA_A = c(rep("A*02:01", 8), rep("A*01:01", 4),
#'             rep("A*02:01", 3), rep("A*01:01", 5)),
#'   stringsAsFactors = FALSE
#' )
#' results <- runHLAassociation(metaclone_data, hla_data)
#'
runHLAassociation <- function(metaclone_data,
                              hla_data,
                              by = "barcode",
                              fdr_threshold = 0.05) {

  if (!by %in% names(metaclone_data)) {
    stop("Column '", by, "' not found in metaclone_data")
  }
  if (!by %in% names(hla_data)) {
    stop("Column '", by, "' not found in hla_data")
  }

  # Merge data
  merged <- merge(metaclone_data, hla_data, by = by)

  if (nrow(merged) == 0) {
    stop("No matching samples between metaclone_data and hla_data")
  }

  # Identify HLA columns (typically A, B, C, DRB1, etc.)
  hla_cols <- grep("^HLA|^A\\*|^B\\*|^C\\*|^DR|^DQ|^DP", names(hla_data), value = TRUE)
  hla_cols <- setdiff(hla_cols, by)

  if (length(hla_cols) == 0) {
    stop("No HLA columns found in hla_data")
  }

  message("Testing associations for ", length(unique(merged$metaclone)),
          " metaclones and ", length(hla_cols), " HLA alleles")

  # Run association tests
  results <- list()

  # Only test non-NA metaclones
  unique_mc <- unique(merged$metaclone)
  unique_mc <- unique_mc[!is.na(unique_mc)]

  for (mc in unique_mc) {
    mc_data <- merged[!is.na(merged$metaclone) & merged$metaclone == mc, ]

    for (hla in hla_cols) {
      # Create contingency table
      in_mc <- !is.na(merged$metaclone) & merged$metaclone == mc
      hla_vals <- merged[[hla]]
      if (is.logical(hla_vals)) {
        has_hla <- !is.na(hla_vals) & hla_vals
      } else {
        has_hla <- !is.na(hla_vals) & hla_vals != ""
      }

      cont_table <- table(in_mc, has_hla)

      # Skip if table is degenerate
      if (nrow(cont_table) < 2 || ncol(cont_table) < 2) next

      # Fisher's exact test
      test <- tryCatch(
        stats::fisher.test(cont_table),
        error = function(e) NULL
      )

      if (!is.null(test)) {
        results[[length(results) + 1]] <- data.frame(
          metaclone = mc,
          hla_allele = hla,
          odds_ratio = test$estimate,
          pvalue = test$p.value,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if (length(results) == 0) {
    warning("No valid association tests could be performed")
    return(NULL)
  }

  result_df <- do.call(rbind, results)

  # FDR correction
  result_df$fdr <- stats::p.adjust(result_df$pvalue, method = "BH")
  result_df$significant <- result_df$fdr < fdr_threshold

  # Sort by FDR
  result_df <- result_df[order(result_df$fdr), ]
  rownames(result_df) <- NULL

  return(result_df)
}
