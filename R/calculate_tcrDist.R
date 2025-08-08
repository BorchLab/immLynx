#' Calculate TCR Distances using tcrdist3
#'
#' This function provides a wrapper around the `tcrdist3` Python package to compute
#' pairwise distances between T-cell receptors (TCRs).
#'
#' @param df A data frame containing TCR information. The data frame should
#'   include columns necessary for `tcrdist3`, such as 'v_b_gene', 'j_b_gene',
#'   'cdr3_b_aa', etc.
#' @param organism A string specifying the organism, either 'human' or 'mouse'.
#' @param chains A character vector indicating the TCR chains to be used for
#'   distance calculation, e.g., `c('alpha', 'beta')`.
#' @param compute_distances A logical indicating whether to compute the full
#'   TCRdist measure. Defaults to `TRUE`.
#'
#' @return A list containing the pairwise distance matrices for the specified chains.
#'   The list may include `pw_alpha`, `pw_beta`, `pw_cdr3_a_aa`, and `pw_cdr3_b_aa`.
#'
#' @export
#' @importFrom basilisk basiliskRun
#' @importFrom reticulate import r_to_py
#'
#' @examples
#' \dontrun{
#'   # Create a sample data frame
#'   tcr_df <- data.frame(
#'     count = c(1, 1),
#'     v_b_gene = c("TRBV19*01", "TRBV19*01"),
#'     j_b_gene = c("TRBJ2-7*01", "TRBJ2-7*01"),
#'     cdr3_b_aa = c("CASSLSAGGAYNEQFF", "CASSLSAGGAYNEQFF")
#'   )
#'
#'   # Calculate TCR distances
#'   distances <- calculate.tcrDist(
#'     df = tcr_df,
#'     organism = 'human',
#'     chains = 'beta'
#'   )
#' }
calculate.tcrDist <- function(df, organism, chains, compute_distances = TRUE) {

  message("Calculating TCR distances with tcrdist3...")

  if (!is.data.frame(df)) {
    stop("Input 'df' must be a data.frame.")
  }

  result <- tryCatch({
    basiliskRun(env = immLynxEnv, fun = function(df, organism, chains, compute_distances) {
      pd <- reticulate::import("pandas")
      TCRrep <- reticulate::import("tcrdist.repertoire")$TCRrep

      pd_df <- reticulate::r_to_py(df)

      tr <- TCRrep(cell_df = pd_df,
                   organism = organism,
                   chains = chains,
                   compute_distances = compute_distances)

      result_list <- list()
      if ('alpha' %in% chains) {
        result_list$pw_alpha <- tr$pw_alpha
        result_list$pw_cdr3_a_aa <- tr$pw_cdr3_a_aa
      }
      if ('beta' %in% chains) {
        result_list$pw_beta <- tr$pw_beta
        result_list$pw_cdr3_b_aa <- tr$pw_cdr3_b_aa
      }

      return(result_list)
    }, df = df, organism = organism, chains = chains, compute_distances = compute_distances)
  }, error = function(e) {
    stop(paste("An error occurred during tcrdist calculation:", e$message))
  })

  message("TCR distance calculation complete.")
  return(result)
}
