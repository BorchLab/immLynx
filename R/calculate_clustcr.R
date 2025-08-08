#' Cluster CDR3 Sequences using clusTCR
#'
#' This function is a wrapper for the `clusTCR` Python package. It performs
#' clustering of Complementarity-Determining Region 3 (CDR3) sequences.
#'
#' @param sequences A character vector of CDR3 amino acid sequences.
#' @param method A string specifying the clustering method to use.
#'   Supported values are 'mcl' (default) and 'dbscan'.
#' @param ... Additional arguments to be passed to the clustering method.
#'   For `method = 'mcl'`:
#'   \itemize{
#'     \item `inflation`: A float for the MCL inflation parameter (default: 2.0).
#'   }
#'   For `method = 'dbscan'`:
#'   \itemize{
#'     \item `eps`: A float for the DBSCAN epsilon parameter (default: 0.5).
#'     \item `min_samples`: An integer for the DBSCAN min_samples parameter (default: 5).
#'   }
#'
#' @return A data frame with clustering results, where each row corresponds to a
#'   sequence and columns include the cluster assignment.
#'
#' @export
#' @importFrom basilisk basiliskRun
#' @importFrom reticulate import
#'
#' @examples
#' \dontrun{
#'   # Sample CDR3 sequences
#'   seqs <- c("CASSLAGGREQYF", "CASSLSFGREQYF", "CASSIWSGREQYF", "CASSLGGRYNEQFF")
#'
#'   # Cluster using MCL
#'   clusters_mcl <- calculate.clustcr(sequences = seqs, method = "mcl", inflation = 2.5)
#'
#'   # Cluster using DBSCAN
#'   clusters_dbscan <- calculate.clustcr(sequences = seqs, method = "dbscan", eps = 0.6)
#' }
calculate.clustcr <- function(sequences, method = "mcl", ...) {

  message("Clustering sequences with clusTCR...")

  if (!is.character(sequences)) {
    stop("Input 'sequences' must be a character vector.")
  }

  args <- list(...)

  result <- tryCatch({
    basiliskRun(env = immLynxEnv, fun = function(sequences, method, args) {
      Clustering <- reticulate::import("clustcr")$Clustering

      cl <- Clustering(sequences = sequences)

      run_args <- c(list(method = method), args)

      do.call(cl$run, run_args)

      clusters <- cl$get_clusters()
      return(clusters)

    }, sequences = sequences, method = method, args = args)
  }, error = function(e) {
    stop(paste("An error occurred during clusTCR execution:", e$message))
  })

  message("clusTCR execution complete.")
  return(result)
}
