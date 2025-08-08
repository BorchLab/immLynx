#' Run DeepTCR Unsupervised Feature Extraction
#'
#' This function is a wrapper for a common workflow in the `DeepTCR` Python package.
#' It initializes a `DeepTCR` object, loads data, trains a Variational Autoencoder (VAE),
#' and extracts features from the data.
#'
#' @param output_dir A string specifying the directory to save `DeepTCR` models and results.
#' @param data_dir A string specifying the directory where the input data files are located.
#' @param files A character vector of file names to be loaded.
#' @param latent_dim An integer specifying the dimensionality of the latent space for the VAE.
#' @param epochs An integer specifying the number of training epochs.
#'
#' @return A data frame containing the extracted features.
#'
#' @export
#' @importFrom basilisk basiliskRun
#' @importFrom reticulate import
#'
#' @examples
#' \dontrun{
#'   # This is a placeholder example as it requires specific data files.
#'   # Create dummy directories and files for demonstration.
#'   dir.create("deeptcr_output")
#'   dir.create("deeptcr_data")
#'   write.csv(data.frame(beta = "CASSLSAGGAYNEQFF"), file = "deeptcr_data/sample1.csv")
#'
#'   features <- calculate.deepTCR(
#'     output_dir = "deeptcr_output",
#'     data_dir = "deeptcr_data",
#'     files = c("sample1.csv")
#'   )
#' }
calculate.deepTCR <- function(output_dir, data_dir, files, latent_dim = 100, epochs = 100) {

  message("Running DeepTCR workflow...")

  if (!dir.exists(output_dir)) {
    message("Output directory does not exist. Creating it.")
    dir.create(output_dir, recursive = TRUE)
  }

  result <- tryCatch({
    basiliskRun(env = immLynxEnv, fun = function(output_dir, data_dir, files, latent_dim, epochs) {
      DeepTCR_U <- reticulate::import("DeepTCR.DeepTCR_U")$DeepTCR_U

      # Initialize DeepTCR object
      DTCR <- DeepTCR_U(output_dir)

      # Load data
      DTCR$Load_Data(directory = data_dir, files = files)

      # Train VAE
      DTCR$Train_VAE(latent_dim = latent_dim, epochs = epochs)

      # Get features
      feats <- DTCR$Get_Feats()

      return(feats)
    }, output_dir = output_dir, data_dir = data_dir, files = files, latent_dim = latent_dim, epochs = epochs)
  }, error = function(e) {
    stop(paste("An error occurred during DeepTCR execution:", e$message))
  })

  message("DeepTCR workflow complete.")
  return(result)
}
