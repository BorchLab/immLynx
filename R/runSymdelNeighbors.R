#' Find Near-Neighbor CDR3 Sequences by Symmetric Deletion Lookup
#'
#' @description Extracts CDR3 amino acid sequences from a
#'   SingleCellExperiment object and finds all pairs within a given edit
#'   distance using the symmetric deletion lookup ("symdel") algorithm from
#'   pyrepseq. This is the same algorithm described in the XT-neighbor
#'   preprint, running on CPU rather than GPU.
#'
#' @param input A SingleCellExperiment object containing scRepertoire TCR data.
#' @param chains Which chain(s) to search: "TRB", "TRA", or "both".
#'   Default is "TRB".
#' @param max_edits Maximum edit distance defining a neighbor. Default is 1.
#' @param max_returns Maximum number of neighbors to return per sequence.
#'   Default is NULL, meaning no limit.
#' @param n_cpu Number of CPU processes for the search. Default is 1.
#' @param combine_chains Logical. If TRUE and chains="both", concatenates alpha
#'   and beta sequences with "_" before searching. Default is FALSE.
#' @param return_object Logical. If TRUE, adds a per-cell neighbor count to the
#'   input object. If FALSE, returns the neighbor edge list. Default is TRUE.
#' @param column_prefix Prefix for the new metadata column. Default is "symdel",
#'   producing a "symdel_degree" column.
#' @param ... Additional arguments passed to calculate.symdel().
#'
#' @return If return_object=TRUE, the input object with a `<prefix>_degree`
#'   column added to colData giving the number of neighbors found for each
#'   cell's CDR3 sequence. Cells with no sequence for the requested chain get
#'   NA. If return_object=FALSE, a data.frame with columns `from_seq`,
#'   `to_seq`, and `distance`, one row per undirected neighbor pair.
#'
#' @details The search runs over unique sequences rather than over cells.
#'   Repertoires carry heavy clonal redundancy, so this is both substantially
#'   faster and what the symdel algorithm expects. Results are mapped back to
#'   cells afterward.
#'
#'   Degree is a neighborhood-density measure. High-degree sequences sit inside
#'   dense clusters of similar receptors, which is the signal used to identify
#'   convergent recombination and antigen-driven expansion. For explicit
#'   cluster assignments, see \code{\link{runClustTCR}}.
#'
#' @export
#' @importFrom immApex getIR
#' @importFrom methods is
#' @importFrom SummarizedExperiment colData colData<-
#'
#' @examples
#' data(immLynx_example)
#' \donttest{
#'   # Neighbors within one edit on the TRB chain
#'   sce <- runSymdelNeighbors(immLynx_example, chains = "TRB")
#'
#'   # Widen the search to two edits
#'   sce <- runSymdelNeighbors(immLynx_example,
#'                             chains = "TRB",
#'                             max_edits = 2)
#'
#'   # Get the neighbor edge list instead
#'   edges <- runSymdelNeighbors(immLynx_example,
#'                               chains = "TRB",
#'                               return_object = FALSE)
#'   head(edges)
#' }
runSymdelNeighbors <- function(input,
                               chains = c("TRB", "TRA", "both"),
                               max_edits = 1,
                               max_returns = NULL,
                               n_cpu = 1,
                               combine_chains = FALSE,
                               return_object = TRUE,
                               column_prefix = "symdel",
                               ...) {

  chains <- match.arg(chains)

  if (!methods::is(input, "SingleCellExperiment")) {
    stop("Input must be a SingleCellExperiment object")
  }

  if (!is.numeric(max_edits) || length(max_edits) != 1L ||
      is.na(max_edits) || max_edits < 1) {
    stop("max_edits must be a single positive integer")
  }
  max_edits <- as.integer(max_edits)

  message("Extracting TCR sequences from object...")
  seq_map <- .extractChainSeqs(input, chains = chains,
                               combine_chains = combine_chains)

  unique_seqs <- unique(seq_map$sequences[!is.na(seq_map$sequences)])
  if (length(unique_seqs) == 0L) {
    stop("No CDR3 sequences found for the requested chain(s).")
  }

  message("Searching ", length(unique_seqs), " unique sequences ",
          "(max_edits = ", max_edits, ")...")

  triplets <- calculate.symdel(unique_seqs,
                               max_edits = max_edits,
                               max_returns = max_returns,
                               n_cpu = n_cpu,
                               ...)

  edges <- .symdelTripletsToEdges(triplets, unique_seqs)
  message("Found ", nrow(edges), " neighbor pair(s).")

  if (!return_object) {
    return(edges)
  }

  deg <- .symdelDegree(edges, unique_seqs)

  # Map degrees from unique sequences back onto cells. Cells with no sequence
  # for the requested chain stay NA.
  cell_deg <- unname(deg[match(seq_map$sequences, unique_seqs)])

  col_vec <- rep(NA_integer_, ncol(input))
  names(col_vec) <- colnames(input)
  col_vec[seq_map$barcodes] <- cell_deg

  colData(input)[[paste0(column_prefix, "_degree")]] <- unname(col_vec)
  input
}


#' Convert symdel triplets into an undirected edge data.frame
#'
#' @description symdel returns (i, j, distance) using 0-based indices and emits
#'   both directions of every pair. This collapses them to one row per
#'   undirected pair and drops any self-pair.
#'
#' @param triplets A list of length-3 numeric vectors, or an n x 3 matrix.
#' @param unique_seqs Character vector the indices refer to.
#' @return A data.frame with `from_seq`, `to_seq`, and `distance`.
#' @keywords internal
.symdelTripletsToEdges <- function(triplets, unique_seqs) {

  empty <- data.frame(from_seq = character(0), to_seq = character(0),
                      distance = integer(0), stringsAsFactors = FALSE)

  if (is.null(triplets) || length(triplets) == 0L) {
    return(empty)
  }

  m <- if (is.matrix(triplets)) triplets else do.call(rbind, triplets)
  if (is.null(m) || nrow(m) == 0L) {
    return(empty)
  }

  # 0-based Python indices to 1-based R indices.
  i <- as.integer(m[, 1L]) + 1L
  j <- as.integer(m[, 2L]) + 1L
  d <- as.integer(m[, 3L])

  # Drop self-pairs, then canonicalize orientation so the two directions of a
  # pair collapse onto the same row.
  keep <- i != j
  i <- i[keep]; j <- j[keep]; d <- d[keep]

  lo <- pmin(i, j)
  hi <- pmax(i, j)
  dedup <- !duplicated(paste(lo, hi))

  out <- data.frame(
    from_seq = unique_seqs[lo[dedup]],
    to_seq   = unique_seqs[hi[dedup]],
    distance = d[dedup],
    stringsAsFactors = FALSE
  )
  rownames(out) <- NULL
  out
}


#' Neighbor count per unique sequence
#'
#' @param edges A data.frame from .symdelTripletsToEdges().
#' @param unique_seqs Character vector of every sequence searched.
#' @return A named integer vector of degrees, in the order of `unique_seqs`.
#' @keywords internal
.symdelDegree <- function(edges, unique_seqs) {

  deg <- integer(length(unique_seqs))
  names(deg) <- unique_seqs

  if (nrow(edges) > 0L) {
    # Each undirected edge contributes one to both endpoints.
    counts <- table(c(edges$from_seq, edges$to_seq))
    deg[names(counts)] <- as.integer(counts)
  }

  deg
}
