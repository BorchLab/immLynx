#' Infer TCR Selection Pressures using soNNia
#'
#' This function is a wrapper for the `soNNia` Python package, which supersedes
#' the original `sonia` package. It infers selection pressures on T-cell receptors (TCRs)
#' by comparing a dataset of selected TCRs against a background of unselected TCRs.
#'
#' Note: The `soNNia` package expects a background dataset of unselected TCRs,
#' which can be generated using the `calculate.olga` function with `action = "generate"`.
#' The output of `calculate.olga` should be saved to a file, and the path to this
#' file should be provided as the `pgen_filename`.
#'
#' @param data_folder A string specifying the directory where the input data file is located.
#' @param data_filename A string specifying the name of the data file with selected TCRs.
#' @param pgen_filename A string specifying the name of the file with the background TCRs.
#' @param organism A string specifying the organism, e.g., 'human' or 'mouse'.
#' @param dataset_type A string specifying the type of dataset, e.g., 'TCR' or 'BCR'.
#' @param v_beta A logical indicating whether to use V-beta gene information.
#' @param j_beta A logical indicating whether to use J-beta gene information.
#' @param CDR3_beta A logical indicating whether to use CDR3-beta sequence information.
#' @param n_epochs An integer specifying the number of training epochs.
#' @param save_folder A string specifying the directory to save the `soNNia` model and results.
#'
#' @return A named list of data frames containing the selection inference results.
#'
#' @export
#' @importFrom basilisk basiliskRun
#' @importFrom reticulate import
#'
#' @examples
#' \dontrun{
#'   # This is a placeholder example as it requires specific data files.
#'   # 1. Generate a background dataset with OLGA
#'   olga_data <- calculate.olga(action = "generate", model = "humanTRB", n = 1000)
#'   write.csv(olga_data, "pgen_data.csv")
#'
#'   # 2. Prepare your selected TCR dataset
#'   selected_data <- data.frame(cdr3_b_aa = "CASSL...", v_b_gene = "TRBV1...", ...)
#'   write.csv(selected_data, "selected_data.csv")
#'
#'   # 3. Run soNNia
#'   results <- calculate.sonia(
#'     data_folder = getwd(),
#'     data_filename = "selected_data.csv",
#'     pgen_filename = "pgen_data.csv",
#'     organism = "human",
#'     dataset_type = "TCR"
#'   )
#' }
calculate.sonia <- function(data_folder, data_filename, pgen_filename, organism,
                            dataset_type = "TCR", v_beta = TRUE, j_beta = TRUE, CDR3_beta = TRUE,
                            n_epochs = 100, save_folder = "sonia_results") {

  message("Running soNNia selection analysis...")

  if (!dir.exists(save_folder)) {
    message("Save folder does not exist. Creating it.")
    dir.create(save_folder, recursive = TRUE)
  }

  result <- tryCatch({
    basiliskRun(env = immLynxEnv, fun = function(data_folder, data_filename, pgen_filename,
                                                 organism, dataset_type, v_beta, j_beta,
                                                 CDR3_beta, n_epochs, save_folder) {

      sonnia_preprocess <- reticulate::import("sonnia.processing.pre_processing")
      sonnia_model <- reticulate::import("sonnia.model.sonnia")

      # Pre-process data
      processed <- sonnia_preprocess$pre_process_data(
        data_folder = data_folder,
        data_filename = data_filename,
        pgen_filename = pgen_filename,
        dataset_type = dataset_type,
        organism = organism
      )

      data <- processed[[1]]
      pgen <- processed[[2]]
      weights <- processed[[3]]

      # Define model
      model <- sonnia_model$SoNNia(
        data = data,
        pgen = pgen,
        weights = weights,
        name = "sonnia_model",
        v_beta = v_beta,
        j_beta = j_beta,
        CDR3_beta = CDR3_beta,
        save_folder = save_folder
      )

      # Train model
      model$train(
        n_epochs = n_epochs
      )

      # Get results
      results <- model$get_results()

      return(results)

    }, data_folder = data_folder, data_filename = data_filename, pgen_filename = pgen_filename,
       organism = organism, dataset_type = dataset_type, v_beta = v_beta, j_beta = j_beta,
       CDR3_beta = CDR3_beta, n_epochs = n_epochs, save_folder = save_folder)
  }, error = function(e) {
    stop(paste("An error occurred during soNNia execution:", e$message))
  })

  message("soNNia analysis complete.")
  return(result)
}
