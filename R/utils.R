#' @title Utility Functions for TCR Data Handling
#' @name utils
#' @description Helper functions for extracting, validating, and converting
#'   TCR data from various sources.
NULL


#' Extract TCR Data from Seurat Object
#'
#' @description Extracts T-cell receptor data from a Seurat object that has been
#'   processed with scRepertoire. This is a convenience wrapper around
#'   \code{immApex::getIR()} that provides additional formatting options.
#'
#' @param input A Seurat object containing scRepertoire TCR data in metadata,
#'   or a SingleCellExperiment object.
#' @param chains Which chains to extract: "TRA", "TRB", "TRG", "TRD",
#'   "IGH", "IGL", "IGK", or "both" (for TRA and TRB). Default is "TRB".
#' @param format Output format: "long" (one row per chain) or "wide" (one row
#'   per cell with columns for each chain). Default is "long".
#' @param remove_na Logical. If TRUE, removes rows with NA CDR3 sequences.
#'   Default is TRUE.
#'
#' @return A data.frame containing:
#'   \describe{
#'     \item{barcode}{Cell barcode}
#'     \item{cdr3_aa}{CDR3 amino acid sequence}
#'     \item{cdr3_nt}{CDR3 nucleotide sequence (if available)}
#'     \item{v}{V gene}
#'     \item{d}{D gene (if applicable)}
#'     \item{j}{J gene}
#'     \item{c}{C gene (if available)}
#'     \item{chain}{Chain type (TRA, TRB, etc.)}
#'   }
#'
#' @export
#' @importFrom immApex getIR
#'
#' @examples
#' # Extract beta chain data
#' data(immLynx_example)
#' tcr_data <- extractTCRdata(immLynx_example, chains = "TRB")
#' head(tcr_data)
#'
extractTCRdata <- function(input,
                           chains = c("TRB", "TRA", "TRG", "TRD",
                                     "IGH", "IGL", "IGK", "both"),
                           format = c("long", "wide"),
                           remove_na = TRUE) {

  chains <- match.arg(chains)
  format <- match.arg(format)

  if (chains == "both") {
    tra_data <- immApex::getIR(input, chains = "TRA")
    trb_data <- immApex::getIR(input, chains = "TRB")

    if (format == "wide") {
      # Merge into wide format
      result <- merge(tra_data, trb_data,
                     by = "barcode",
                     suffixes = c("_TRA", "_TRB"),
                     all = TRUE)
    } else {
      # Stack for long format
      tra_data$chain <- "TRA"
      trb_data$chain <- "TRB"
      result <- rbind(tra_data, trb_data)
    }
  } else {
    result <- immApex::getIR(input, chains = chains)
    result$chain <- chains
  }

  if (remove_na) {
    if (format == "wide" && chains == "both") {
      # For wide format, keep rows with at least one valid chain
      result <- result[!is.na(result$cdr3_aa_TRA) | !is.na(result$cdr3_aa_TRB), ]
    } else {
      result <- result[!is.na(result$cdr3_aa), ]
    }
  }

  return(result)
}


#' Validate TCR Data Format
#'
#' @description Validates that TCR data is in the correct format for immLynx
#'   analysis functions. Checks for required columns, valid sequence formats,
#'   and gene nomenclature.
#'
#' @param tcr_data A data.frame containing TCR data
#' @param check_genes Logical. If TRUE, validates gene names against IMGT
#'   nomenclature. Default is FALSE.
#' @param check_sequences Logical. If TRUE, validates that CDR3 sequences
#'   contain only valid amino acids. Default is TRUE.
#' @param strict Logical. If TRUE, stops with error on validation failure.
#'   If FALSE, returns validation report. Default is FALSE.
#'
#' @return If strict=FALSE, returns a list with:
#'   \describe{
#'     \item{valid}{Logical indicating overall validity}
#'     \item{errors}{Character vector of error messages}
#'     \item{warnings}{Character vector of warning messages}
#'     \item{summary}{Summary statistics of the data}
#'   }
#'   If strict=TRUE, returns TRUE invisibly on success or stops with error.
#'
#' @export
#'
#' @examples
#' # Create example TCR data
#' tcr_data <- data.frame(
#'   barcode = paste0("cell_", 1:5),
#'   cdr3_aa = c("CASSLGTGELFF", "CASSIRSSYEQYF", "CASSYSTGELFF",
#'               "CASNQGLNEKLFF", "CASSLDRNEQFF"),
#'   v = paste0("TRBV", c("7-2", "12-3", "5-1", "28", "7-9")),
#'   j = paste0("TRBJ", c("2-2", "1-1", "2-7", "1-5", "2-1")),
#'   chain = rep("TRB", 5),
#'   stringsAsFactors = FALSE
#' )
#' report <- validateTCRdata(tcr_data, strict = FALSE)
#' report$valid
#' report$summary
#'
validateTCRdata <- function(tcr_data,
                            check_genes = FALSE,
                            check_sequences = TRUE,
                            strict = FALSE) {

  errors <- character()
  warnings <- character()

  # Check required columns
  required_cols <- c("barcode", "cdr3_aa")
  missing_cols <- setdiff(required_cols, names(tcr_data))
  if (length(missing_cols) > 0) {
    errors <- c(errors, paste("Missing required columns:", paste(missing_cols, collapse = ", ")))
  }

  # Check optional but recommended columns
  recommended_cols <- c("v", "j", "chain")
  missing_rec <- setdiff(recommended_cols, names(tcr_data))
  if (length(missing_rec) > 0) {
    warnings <- c(warnings, paste("Missing recommended columns:", paste(missing_rec, collapse = ", ")))
  }

  # Check for valid amino acid sequences
  if (check_sequences && "cdr3_aa" %in% names(tcr_data)) {
    valid_aa <- "^[ACDEFGHIKLMNPQRSTVWY*]+$"
    non_na_seqs <- tcr_data$cdr3_aa[!is.na(tcr_data$cdr3_aa)]

    if (length(non_na_seqs) > 0) {
      invalid_seqs <- !grepl(valid_aa, toupper(non_na_seqs))
      n_invalid <- sum(invalid_seqs)

      if (n_invalid > 0) {
        warnings <- c(warnings, paste(n_invalid, "sequences contain non-standard amino acids"))
      }
    }
  }

  # Check V/J gene format if requested
  if (check_genes) {
    if ("v" %in% names(tcr_data)) {
      # Check for IMGT-style gene names
      v_pattern <- "^(TR[ABDG]V|IG[HKL]V)[0-9]+"
      non_na_v <- tcr_data$v[!is.na(tcr_data$v)]
      if (length(non_na_v) > 0) {
        invalid_v <- !grepl(v_pattern, non_na_v)
        if (sum(invalid_v) > 0) {
          warnings <- c(warnings, paste(sum(invalid_v), "V genes do not match IMGT nomenclature"))
        }
      }
    }

    if ("j" %in% names(tcr_data)) {
      j_pattern <- "^(TR[ABDG]J|IG[HKL]J)[0-9]+"
      non_na_j <- tcr_data$j[!is.na(tcr_data$j)]
      if (length(non_na_j) > 0) {
        invalid_j <- !grepl(j_pattern, non_na_j)
        if (sum(invalid_j) > 0) {
          warnings <- c(warnings, paste(sum(invalid_j), "J genes do not match IMGT nomenclature"))
        }
      }
    }
  }

  # Summary statistics
  summary_stats <- list(
    n_rows = nrow(tcr_data),
    n_unique_barcodes = if ("barcode" %in% names(tcr_data)) length(unique(tcr_data$barcode)) else NA,
    n_unique_cdr3 = if ("cdr3_aa" %in% names(tcr_data)) length(unique(tcr_data$cdr3_aa[!is.na(tcr_data$cdr3_aa)])) else NA,
    n_na_cdr3 = if ("cdr3_aa" %in% names(tcr_data)) sum(is.na(tcr_data$cdr3_aa)) else NA,
    chains = if ("chain" %in% names(tcr_data)) unique(tcr_data$chain) else NA
  )

  valid <- length(errors) == 0

  if (strict) {
    if (!valid) {
      stop("TCR data validation:\n",
           paste(errors, collapse = "\n"))
    }
    if (length(warnings) > 0) {
      warning("TCR data validation:\n",
              paste(warnings, collapse = "\n"))
    }
    return(invisible(TRUE))
  }

  list(
    valid = valid,
    errors = errors,
    warnings = warnings,
    summary = summary_stats
  )
}


#' Convert TCR Data to tcrdist3 Format
#'
#' @description Converts TCR data from immLynx/scRepertoire format to the format
#'   required by tcrdist3 for distance calculations.
#'
#' @param tcr_data A data.frame from extractTCRdata() or similar source
#' @param chains Which chains to include: "alpha", "beta", or "both".
#'   Default is "beta".
#' @param include_count Logical. If TRUE, adds a 'count' column (default 1).
#'   Default is TRUE.
#'
#' @return A data.frame in tcrdist3 format with columns:
#'   \describe{
#'     \item{count}{Clone count (default 1)}
#'     \item{v_a_gene}{Alpha V gene (if alpha chain included)}
#'     \item{j_a_gene}{Alpha J gene (if alpha chain included)}
#'     \item{cdr3_a_aa}{Alpha CDR3 amino acid sequence (if alpha chain included)}
#'     \item{v_b_gene}{Beta V gene (if beta chain included)}
#'     \item{j_b_gene}{Beta J gene (if beta chain included)}
#'     \item{cdr3_b_aa}{Beta CDR3 amino acid sequence (if beta chain included)}
#'   }
#'
#' @export
#'
#' @examples
#' # Convert long-format TCR data to tcrdist3 format
#' tcr_data <- data.frame(
#'   barcode = paste0("cell_", 1:5),
#'   cdr3_aa = c("CASSLGTGELFF", "CASSIRSSYEQYF", "CASSYSTGELFF",
#'               "CASNQGLNEKLFF", "CASSLDRNEQFF"),
#'   v = paste0("TRBV", c("7-2", "12-3", "5-1", "28", "7-9")),
#'   j = paste0("TRBJ", c("2-2", "1-1", "2-7", "1-5", "2-1")),
#'   chain = rep("TRB", 5),
#'   stringsAsFactors = FALSE
#' )
#' tcrdist_format <- convertToTcrdist(tcr_data, chains = "beta")
#' head(tcrdist_format)
#'
convertToTcrdist <- function(tcr_data,
                             chains = c("beta", "alpha", "both"),
                             include_count = TRUE) {

  chains <- match.arg(chains)

  # Determine if data is in wide or long format
  is_wide <- any(grepl("_TRA$|_TRB$", names(tcr_data)))

  if (is_wide) {
    # Wide format: columns like cdr3_aa_TRA, cdr3_aa_TRB
    result <- data.frame(row.names = seq_len(nrow(tcr_data)))

    if (chains %in% c("alpha", "both")) {
      result$v_a_gene <- tcr_data$v_TRA
      result$j_a_gene <- tcr_data$j_TRA
      result$cdr3_a_aa <- tcr_data$cdr3_aa_TRA
    }
    if (chains %in% c("beta", "both")) {
      result$v_b_gene <- tcr_data$v_TRB
      result$j_b_gene <- tcr_data$j_TRB
      result$cdr3_b_aa <- tcr_data$cdr3_aa_TRB
    }

  } else {
    # Long format: need to pivot
    if (!("chain" %in% names(tcr_data))) {
      # Assume single chain based on request
      if (chains == "alpha") {
        result <- data.frame(
          v_a_gene = tcr_data$v,
          j_a_gene = tcr_data$j,
          cdr3_a_aa = tcr_data$cdr3_aa,
          stringsAsFactors = FALSE
        )
      } else if (chains == "beta") {
        result <- data.frame(
          v_b_gene = tcr_data$v,
          j_b_gene = tcr_data$j,
          cdr3_b_aa = tcr_data$cdr3_aa,
          stringsAsFactors = FALSE
        )
      } else {
        stop("Cannot determine chain type. Please ensure data has 'chain' column or use format='wide'")
      }
    } else {
      # Has chain column - filter and pivot
      if (chains == "alpha") {
        alpha_data <- tcr_data[tcr_data$chain == "TRA", ]
        result <- data.frame(
          v_a_gene = alpha_data$v,
          j_a_gene = alpha_data$j,
          cdr3_a_aa = alpha_data$cdr3_aa,
          stringsAsFactors = FALSE
        )
      } else if (chains == "beta") {
        beta_data <- tcr_data[tcr_data$chain == "TRB", ]
        result <- data.frame(
          v_b_gene = beta_data$v,
          j_b_gene = beta_data$j,
          cdr3_b_aa = beta_data$cdr3_aa,
          stringsAsFactors = FALSE
        )
      } else {
        # Both chains - need to pivot
        alpha_data <- tcr_data[tcr_data$chain == "TRA", c("barcode", "v", "j", "cdr3_aa")]
        beta_data <- tcr_data[tcr_data$chain == "TRB", c("barcode", "v", "j", "cdr3_aa")]

        names(alpha_data) <- c("barcode", "v_a_gene", "j_a_gene", "cdr3_a_aa")
        names(beta_data) <- c("barcode", "v_b_gene", "j_b_gene", "cdr3_b_aa")

        result <- merge(alpha_data, beta_data, by = "barcode", all = TRUE)
        result$barcode <- NULL
      }
    }
  }

  if (include_count) {
    result$count <- 1L
    # Move count to first position
    result <- result[, c("count", setdiff(names(result), "count")), drop = FALSE]
  }

  rownames(result) <- NULL
  result
}


#' Summarize TCR Repertoire Statistics
#'
#' @description Generates summary statistics for TCR repertoire data including
#'   diversity metrics, clonality measures, and sequence characteristics.
#'
#' @param input A Seurat object with scRepertoire data, or a data.frame from
#'   extractTCRdata().
#' @param chains Which chains to summarize: "TRB", "TRA", or "both".
#'   Default is "TRB".
#' @param group.by Optional metadata column for grouping (Seurat objects only).
#' @param calculate_diversity Logical. If TRUE, calculates diversity indices.
#'   Default is TRUE.
#'
#' @return An object of class "TCR_summary" containing:
#'   \describe{
#'     \item{total_cells}{Total number of cells with TCR data}
#'     \item{unique_clonotypes}{Number of unique clonotypes}
#'     \item{clonality}{1 - normalized Shannon entropy}
#'     \item{diversity}{List of diversity indices (Shannon, Simpson, etc.)}
#'     \item{top_clones}{Top 10 most frequent clonotypes}
#'     \item{cdr3_length}{Summary of CDR3 sequence lengths}
#'     \item{gene_usage}{V and J gene usage frequencies}
#'   }
#'
#' @export
#' @importFrom methods is
#' @importFrom stats median sd
#'
#' @examples
#' # Summarize from a data.frame
#' tcr_data <- data.frame(
#'   barcode = paste0("cell_", 1:10),
#'   cdr3_aa = c("CASSLGTGELFF", "CASSIRSSYEQYF", "CASSLGTGELFF",
#'               "CASSYSTGELFF", "CASSIRSSYEQYF", "CASSLGTGELFF",
#'               "CASNQGLNEKLFF", "CASSYSTGELFF", "CASSLGTGELFF",
#'               "CASSIRSSYEQYF"),
#'   v = rep("TRBV7-2", 10),
#'   j = rep("TRBJ2-2", 10),
#'   chain = rep("TRB", 10),
#'   stringsAsFactors = FALSE
#' )
#' summary <- summarizeTCRrepertoire(tcr_data)
#' print(summary)
#'
summarizeTCRrepertoire <- function(input,
                                   chains = c("TRB", "TRA", "both"),
                                   group.by = NULL,
                                   calculate_diversity = TRUE) {

  chains <- match.arg(chains)

  # Extract data if Seurat object
  if (methods::is(input, "Seurat") || methods::is(input, "SingleCellExperiment")) {
    tcr_data <- extractTCRdata(input, chains = if (chains == "both") "both" else chains,
                               format = "long", remove_na = TRUE)
  } else {
    tcr_data <- input
  }

  # Calculate basic statistics
  total_cells <- nrow(tcr_data)
  unique_cdr3 <- unique(tcr_data$cdr3_aa)
  unique_clonotypes <- length(unique_cdr3)

  # Clonotype frequencies
  clone_freq <- table(tcr_data$cdr3_aa)
  clone_freq <- sort(clone_freq, decreasing = TRUE)

  # Top clones
  top_n <- min(10, length(clone_freq))
  top_idx <- seq_len(top_n)
  top_clones <- data.frame(
    cdr3_aa = names(clone_freq)[top_idx],
    count = as.integer(clone_freq[top_idx]),
    proportion = as.numeric(clone_freq[top_idx]) / total_cells,
    stringsAsFactors = FALSE
  )

  # CDR3 length distribution
  cdr3_lengths <- nchar(tcr_data$cdr3_aa)
  cdr3_length_summary <- list(
    min = min(cdr3_lengths, na.rm = TRUE),
    max = max(cdr3_lengths, na.rm = TRUE),
    mean = mean(cdr3_lengths, na.rm = TRUE),
    median = median(cdr3_lengths, na.rm = TRUE),
    sd = sd(cdr3_lengths, na.rm = TRUE)
  )

  # Gene usage
  gene_usage <- list()
  if ("v" %in% names(tcr_data)) {
    v_usage <- table(tcr_data$v[!is.na(tcr_data$v)])
    v_usage <- sort(v_usage, decreasing = TRUE)
    gene_usage$v_genes <- data.frame(
      gene = names(v_usage),
      count = as.integer(v_usage),
      proportion = as.numeric(v_usage) / sum(v_usage),
      stringsAsFactors = FALSE
    )
  }
  if ("j" %in% names(tcr_data)) {
    j_usage <- table(tcr_data$j[!is.na(tcr_data$j)])
    j_usage <- sort(j_usage, decreasing = TRUE)
    gene_usage$j_genes <- data.frame(
      gene = names(j_usage),
      count = as.integer(j_usage),
      proportion = as.numeric(j_usage) / sum(j_usage),
      stringsAsFactors = FALSE
    )
  }

  # Diversity metrics
  diversity <- NULL
  if (calculate_diversity) {
    p <- as.numeric(clone_freq) / sum(clone_freq)

    # Shannon entropy
    shannon <- -sum(p * log(p))

    # Normalized Shannon entropy (clonality complement)
    shannon_norm <- if (length(p) > 1) shannon / log(length(p)) else 0

    # Simpson index
    simpson <- sum(p^2)

    # Inverse Simpson
    inv_simpson <- 1 / simpson

    # Gini-Simpson index
    gini_simpson <- 1 - simpson

    # Clonality (1 - normalized entropy)
    clonality <- 1 - shannon_norm

    diversity <- list(
      shannon = shannon,
      shannon_normalized = shannon_norm,
      simpson = simpson,
      inverse_simpson = inv_simpson,
      gini_simpson = gini_simpson,
      clonality = clonality,
      richness = unique_clonotypes
    )
  }

  # Build result object
  result <- list(
    total_cells = total_cells,
    unique_clonotypes = unique_clonotypes,
    clonotype_ratio = unique_clonotypes / total_cells,
    diversity = diversity,
    top_clones = top_clones,
    cdr3_length = cdr3_length_summary,
    gene_usage = gene_usage,
    chains = chains
  )

  class(result) <- c("TCR_summary", "list")

  return(result)
}


#' Print Method for TCR_summary Objects
#'
#' @description Custom print method for TCR repertoire summary objects.
#'
#' @param x An object of class "TCR_summary"
#' @param ... Additional arguments (ignored)
#'
#' @return Invisibly returns the input object
#' @export
#' @method print TCR_summary
#' @importFrom utils head
#'
#' @examples
#' tcr_data <- data.frame(
#'   barcode = paste0("cell_", 1:10),
#'   cdr3_aa = c("CASSLGTGELFF", "CASSIRSSYEQYF", "CASSLGTGELFF",
#'               "CASSYSTGELFF", "CASSIRSSYEQYF", "CASSLGTGELFF",
#'               "CASNQGLNEKLFF", "CASSYSTGELFF", "CASSLGTGELFF",
#'               "CASSIRSSYEQYF"),
#'   v = rep("TRBV7-2", 10),
#'   j = rep("TRBJ2-2", 10),
#'   chain = rep("TRB", 10),
#'   stringsAsFactors = FALSE
#' )
#' summary_obj <- summarizeTCRrepertoire(tcr_data)
#' print(summary_obj)
#'
print.TCR_summary <- function(x, ...) {
  cat("=== TCR Repertoire Summary ===\n")
  cat("Chain(s):", x$chains, "\n\n")

  cat("--- Basic Statistics ---\n")
  cat("Total cells with TCR:", x$total_cells, "\n")
  cat("Unique clonotypes:", x$unique_clonotypes, "\n")
  cat("Clonotype ratio:", round(x$clonotype_ratio, 4), "\n\n")

  if (!is.null(x$diversity)) {
    cat("--- Diversity Metrics ---\n")
    cat("Shannon entropy:", round(x$diversity$shannon, 4), "\n")
    cat("Clonality:", round(x$diversity$clonality, 4), "\n")
    cat("Simpson index:", round(x$diversity$simpson, 4), "\n")
    cat("Inverse Simpson:", round(x$diversity$inverse_simpson, 2), "\n\n")
  }

  cat("--- Top Clonotypes ---\n")
  print(head(x$top_clones, 5), row.names = FALSE)
  cat("\n")

  cat("--- CDR3 Length Distribution ---\n")
  cat("Mean:", round(x$cdr3_length$mean, 1), "\n")
  cat("Median:", x$cdr3_length$median, "\n")
  cat("Range:", x$cdr3_length$min, "-", x$cdr3_length$max, "\n")

  invisible(x)
}
