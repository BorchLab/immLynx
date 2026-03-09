#' Get Protein Embeddings from a Model
#'
#' @description Applies a pre-trained model to a batch of tokenized sequences
#'   to generate embeddings. This is the core embedding function used by
#'   \code{\link{runEmbeddings}}.
#'
#' @param model HF model (from AutoModel or similar), typically obtained via
#'   \code{\link{huggingModel}}.
#' @param tokenized.batch A *list* of tokenized tensors OR a list of such lists
#'   (i.e., already minibatched). If you pass a single big batch, set chunk_size.
#'   Typically obtained via \code{\link{tokenizeSequences}}.
#' @param pool One of "mean", "cls", or "none". "mean" is recommended for
#'   sequence-level embeddings.
#' @param chunk_size If tokenized.batch is a single big batch, split it into
#'   chunks of this many sequences. Ignored if you pre-batched upstream.
#' @param prefer_dtype One of "float16", "bfloat16", "float32". Lower precision
#'   uses less memory but may reduce accuracy.
#' @param prefer_device One of "auto", "cuda", "mps", "cpu". "auto" will select
#'   the best available device.
#'
#' @return An R matrix [n_sequences x hidden] if pool != "none".
#'         If pool == "none", returns a list of arrays per chunk.
#'
#' @export
#' @importFrom reticulate py_to_r array_reshape
#'
#' @seealso \code{\link{runEmbeddings}} for a higher-level wrapper that works
#'   directly with Seurat/SingleCellExperiment objects.
#'
#' @examples
#' sequences <- c("CASSLGTGELFF", "CASSIRSSYEQYF", "CASSYSTGELFF")
#' \dontrun{
#'   # Full workflow: load model, tokenize, embed
#'   hf_components <- huggingModel()
#'   tokenized <- tokenizeSequences(hf_components$tokenizer,
#'                                  sequences)
#'
#'   # Mean pooling (recommended for sequence-level tasks)
#'   embeddings <- proteinEmbeddings(hf_components$model,
#'                                   tokenized,
#'                                   pool = "mean",
#'                                   chunk_size = 32)
#'   dim(embeddings)  # [n_sequences x hidden_dim]
#'
#'   # CLS token embedding
#'   cls_emb <- proteinEmbeddings(hf_components$model,
#'                                tokenized,
#'                                pool = "cls",
#'                                chunk_size = 32)
#'
#'   # Per-token embeddings (no pooling)
#'   token_emb <- proteinEmbeddings(hf_components$model,
#'                                  tokenized,
#'                                  pool = "none",
#'                                  chunk_size = 32)
#'
#'   # Use GPU with half precision for speed
#'   embeddings_gpu <- proteinEmbeddings(
#'       hf_components$model, tokenized,
#'       pool = "mean", chunk_size = 64,
#'       prefer_device = "cuda",
#'       prefer_dtype = "float16")
#' }
proteinEmbeddings <- function(
    model,
    tokenized.batch,
    pool = c("none","mean", "cls"),
    chunk_size = NULL,
    prefer_dtype = c("float16","bfloat16","float32"),
    prefer_device = c("auto","cuda","mps","cpu")
) {
  stopifnot(!is.null(model), !is.null(tokenized.batch))
  pool <- match.arg(pool)
  prefer_dtype  <- match.arg(prefer_dtype)
  prefer_device <- match.arg(prefer_device)
  
  torch <- reticulate::import("torch", convert = FALSE)
  
  # --- device pick ---
  has_mps <- reticulate::py_has_attr(torch, "backends") &&
    reticulate::py_has_attr(torch$backends, "mps") &&
    as.logical(toupper(torch$backends$mps$is_available()))
  
  device <- switch(
    prefer_device,
    "cuda" = if (as.logical(toupper(torch$cuda$is_available()))) "cuda" else "cpu",
    "mps"  = if (has_mps) "mps" else "cpu",
    "cpu"  = "cpu",
    "auto" = if (as.logical(toupper(torch$cuda$is_available()))) "cuda" else if (has_mps) "mps" else "cpu"
  )
  dev <- torch$device(device)
  
  # --- dtype pick (best-effort) ---
  dtype <- switch(
    prefer_dtype,
    "float16"  = torch$float16,
    "bfloat16" = if (reticulate::py_has_attr(torch, "bfloat16")) torch$bfloat16 else torch$float32,
    "float32"  = torch$float32
  )
  
  # Move model to same device
  model$to(device = dev)
  model$eval()
  
  # Helpers
  is_py_dict <- function(x) reticulate::py_has_attr(x, "keys")
  to_device <- function(batch, dev) {
    torch    <- reticulate::import("torch", convert = FALSE)
    builtins <- reticulate::import_builtins(convert = FALSE)
    
    # 1) Get keys robustly
    if (reticulate::py_has_attr(batch, "keys")) {
      # batch is a Python dict-like (e.g., HF BatchEncoding)
      # Convert dict_keys → list → R char vector
      keys_pylist <- builtins$list(batch$keys())
      keys_r      <- reticulate::py_to_r(keys_pylist)
      keys        <- as.character(unlist(keys_r, use.names = FALSE))
    } else if (is.list(batch)) {
      # batch is an R list
      keys <- names(batch)
    } else {
      stop("Unsupported batch type: expected a Python dict-like or an R list.")
    }
    
    if (is.null(keys) || length(keys) == 0L) return(batch)
    
    # 2) Move tensor-like values to device
    for (k in keys) {
      # Fetch the value by key from either a Python dict or an R list
      v <- if (reticulate::py_has_attr(batch, "__getitem__")) {
        reticulate::py_get_item(batch, k)
      } else {
        batch[[k]]
      }
      
      # If it's a torch tensor (has .to method), move it
      if (!is.null(v) && reticulate::py_has_attr(v, "to")) {
        v2 <- v$to(device = dev)
        if (reticulate::py_has_attr(batch, "__setitem__")) {
          reticulate::py_set_item(batch, k, v2)
        } else {
          batch[[k]] <- v2
        }
      }
    }
    batch
  }
  
  
  # Pooling in torch (stay on device as long as possible)
  pool_torch <- function(lhs, attn, pool_mode) {
    if (pool_mode == "none") return(lhs)               # [B, L, H]
    if (pool_mode == "cls")  return(lhs$select(1L, 0L))# [B, H] take token 0
    # mean pool over valid tokens only
    mask <- attn$unsqueeze(2L)$to(dtype = lhs$dtype)   # [B, L, 1]
    sum_emb <- (lhs * mask)$sum(dim = reticulate::tuple(1L))   # [B, H]
    denom   <- mask$sum(dim = reticulate::tuple(1L))$clamp(min = 1)
    sum_emb / denom
  }
  
  # Chunking logic
  make_chunks <- function(tb, chunk_size) {
    stopifnot(!is.null(tb))
    if (is.null(chunk_size) || chunk_size < 1) stop("chunk_size must be a positive integer")
    
    builtins <- reticulate::import_builtins(convert = FALSE)
    
    is_py_dict <- function(x) reticulate::py_has_attr(x, "keys")
    
    # Already chunked? (list/pylist of dicts)
    if (is.list(tb) && length(tb) > 0 && all(vapply(tb, is_py_dict, TRUE))) return(tb)
    if (reticulate::py_has_attr(tb, "__len__") && !is_py_dict(tb)) return(tb)  # python list
    
    # Single big dict
    if (!is_py_dict(tb)) stop("tokenized.batch must be a Python dict or a list of Python dicts.")
    
    # ---- Key change: get batch size via py_len (size of dim 0) ----
    n <- as.integer(reticulate::py_len(tb$input_ids))
    
    starts <- seq.int(0L, n - 1L, by = as.integer(chunk_size))
    
    lapply(starts, function(s) {
      e <- min(s + as.integer(chunk_size) - 1L, n - 1L)
      slicer <- builtins$slice(s, e + 1L, NULL)
      
      d <- builtins$dict()
      reticulate::py_set_item(d, "input_ids", tb$input_ids[slicer])
      
      if (reticulate::py_has_attr(tb, "attention_mask") && !reticulate::py_is_null_xptr(tb$attention_mask)) {
        reticulate::py_set_item(d, "attention_mask", tb$attention_mask[slicer])
      }
      if (reticulate::py_has_attr(tb, "token_type_ids") && !reticulate::py_is_null_xptr(tb$token_type_ids)) {
        reticulate::py_set_item(d, "token_type_ids", tb$token_type_ids[slicer])
      }
      d
    })
  }
  
  
  chunks <- make_chunks(tokenized.batch, chunk_size)
  
  pooled_list <- list()
  tokenlevel_list <- list()
  
  # Use inference_mode (cheaper than no_grad)
  with(torch$inference_mode(), {
    for (i in seq_along(chunks)) {
      batch <- to_device(chunks[[i]], dev = dev)
      
      # optional autocast (GPU/MPS + reduced dtype); safe no-op elsewhere
      use_autocast <- (device == "cuda") &&
        !identical(dtype, torch$float32) &&
        reticulate::py_has_attr(torch, "autocast")
      if (use_autocast) {
        ctx <- torch$autocast(device_type = device, dtype = dtype)
        ctx$`__enter__`()
      } else {
        ctx <- NULL
      }
      
      out <- model(
        input_ids     = batch$input_ids,
        attention_mask = batch$attention_mask
      )
      
      if (!is.null(ctx)) ctx$`__exit__`(NULL, NULL, NULL)
      
      lhs <- out$last_hidden_state  # [B, L, H] on device
      
      if (pool == "none") {
        arr <- lhs$to(device = torch$device("cpu"))$to(dtype = torch$float32)$numpy()
        tokenlevel_list[[length(tokenlevel_list) + 1L]] <- reticulate::py_to_r(arr)
      } else {
        pooled <- pool_torch(lhs, batch$attention_mask, pool)  # [B, H]
        arr <- pooled$to(device = torch$device("cpu"))$to(dtype = torch$float32)$numpy()
        pooled_list[[length(pooled_list) + 1L]] <- reticulate::py_to_r(arr)
      }
      
      if (device == "cuda") torch$cuda$empty_cache()
    }
  })
  
  if (pool == "none") {
    return(tokenlevel_list)
  } else {
    emb <- do.call(rbind, pooled_list)  # [N, H]
    rownames(emb) <- NULL
    emb
  }
}
