#' Run OLGA for Pgen Calculation or Sequence Generation
#'
#' This function is a wrapper for the `OLGA` Python package. It can be used to
#' either compute the generation probability (Pgen) of CDR3 sequences or to
#' generate new sequences from a pre-trained model.
#'
#' @param action A string specifying the action to perform: 'pgen' or 'generate'.
#' @param model A string specifying the OLGA model to use.
#'   Supported values are 'humanTRB', 'humanTRA', 'humanIGH', 'mouseTRB'.
#' @param ... Additional arguments depending on the action.
#'   For 'pgen': `sequences` (required, character vector), `v_genes` (optional), `j_genes` (optional).
#'   For 'generate': `n` (required, integer, number of sequences to generate).
#'
#' @return For 'pgen', a numeric vector of generation probabilities.
#'   For 'generate', a data frame of generated sequences with their properties.
#'
#' @export
#' @importFrom basilisk basiliskRun
#' @importFrom reticulate import
#'
#' @examples
#' \dontrun{
#'   # Calculate Pgen
#'   pgens <- calculate.olga(
#'     action = "pgen",
#'     model = "humanTRB",
#'     sequences = c("CASSLGRDGGHEQYF", "CASSTGQANYGYTF")
#'   )
#'
#'   # Generate sequences
#'   new_seqs <- calculate.olga(
#'     action = "generate",
#'     model = "humanTRB",
#'     n = 10
#'   )
#' }
calculate.olga <- function(action, model, ...) {

  if (!action %in% c("pgen", "generate")) {
    stop("Invalid action. Must be 'pgen' or 'generate'.")
  }

  args <- list(...)
  message("Running OLGA with action: ", action, "...")

  result <- tryCatch({
    basiliskRun(env = immLynxEnv, fun = function(action, model, args) {
      olga_load <- reticulate::import("olga.load_model")
      olga_pgen <- reticulate::import("olga.generation_probability")
      olga_seqgen <- reticulate::import("olga.sequence_generation")
      olga_pkg <- reticulate::import("olga")

      model_map <- list(
        humanTRB = list(folder = "human_T_beta", type = "VDJ"),
        humanTRA = list(folder = "human_T_alpha", type = "VJ"),
        humanIGH = list(folder = "human_B_heavy", type = "VDJ"),
        mouseTRB = list(folder = "mouse_T_beta", type = "VDJ")
      )

      if (!model %in% names(model_map)) {
        stop("Unsupported model. Please use one of: ", paste(names(model_map), collapse = ", "))
      }

      model_info <- model_map[[model]]
      olga_path <- olga_pkg$`__path__`[[1]]
      model_folder <- file.path(olga_path, "default_models", model_info$folder)

      params_file <- file.path(model_folder, 'model_params.txt')
      marginals_file <- file.path(model_folder, 'model_marginals.txt')
      v_anchor_file <- file.path(model_folder, 'V_gene_CDR3_anchors.csv')
      j_anchor_file <- file.path(model_folder, 'J_gene_CDR3_anchors.csv')

      if (model_info$type == "VDJ") {
        genomic_data <- olga_load$GenomicDataVDJ()
        genomic_data$load_igor_genomic_data(params_file, v_anchor_file, j_anchor_file)
        generative_model <- olga_load$GenerativeModelVDJ()
        generative_model$load_and_process_igor_model(marginals_file)
        pgen_model <- olga_pgen$GenerationProbabilityVDJ(generative_model, genomic_data)
        seq_gen_model <- olga_seqgen$SequenceGenerationVDJ(generative_model, genomic_data)
      } else { # VJ
        genomic_data <- olga_load$GenomicDataVJ()
        genomic_data$load_igor_genomic_data(params_file, v_anchor_file, j_anchor_file)
        generative_model <- olga_load$GenerativeModelVJ()
        generative_model$load_and_process_igor_model(marginals_file)
        pgen_model <- olga_pgen$GenerationProbabilityVJ(generative_model, genomic_data)
        seq_gen_model <- olga_seqgen$SequenceGenerationVJ(generative_model, genomic_data)
      }

      if (action == "pgen") {
        sequences <- args$sequences
        v_genes <- args$v_genes
        j_genes <- args$j_genes

        if (is.null(sequences)) stop("sequences must be provided for action 'pgen'")

        num_seqs <- length(sequences)
        pgens <- numeric(num_seqs)
        for (i in 1:num_seqs) {
            v <- if (!is.null(v_genes)) v_genes[i] else NULL
            j <- if (!is.null(j_genes)) j_genes[i] else NULL
            pgens[i] <- pgen_model$compute_aa_CDR3_pgen(sequences[i], v, j)
        }
        return(pgens)

      } else if (action == "generate") {
        n <- args$n
        if (is.null(n)) stop("n must be provided for action 'generate'")

        results <- list()
        for (i in 1:n) {
            results[[i]] <- seq_gen_model$gen_rnd_prod_CDR3()
        }

        df <- do.call(rbind, lapply(results, function(x) {
            data.frame(nt_seq = x[[1]], aa_seq = x[[2]], v_index = x[[3]], j_index = x[[4]], stringsAsFactors = FALSE)
        }))
        return(df)
      }
    }, action = action, model = model, args = args)
  }, error = function(e) {
    stop(paste("An error occurred during OLGA execution:", e$message))
  })

  message("OLGA execution complete.")
  return(result)
}
