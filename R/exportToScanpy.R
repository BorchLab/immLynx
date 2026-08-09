# Export a SingleCellExperiment / Seurat object to a scanpy/scirpy-compatible
# H5AD or H5MU bundle. Public entry point: exportToScanpy().

#' @keywords internal
# Run a function inside the scanpyExportEnv basilisk env. Mirrors the
# .run_in_basilisk pattern in R/calculate_helpers.R but targets the
# scirpy/muon stack rather than the immLynxEnv stack.
.run_in_scanpy_env <- function(FUN, ...) {
  proc <- basilisk::basiliskStart(scanpyExportEnv)
  on.exit(basilisk::basiliskStop(proc))
  basilisk::basiliskRun(proc, FUN, ...)
}

#' @keywords internal
.pruneColData <- function(sce, obs_columns, drop_ct) {
  cd <- SummarizedExperiment::colData(sce)

  # Always drop list-columns (incompatible with H5AD encoding).
  is_list_col <- vapply(seq_len(ncol(cd)),
                        function(j) is.list(cd[[j]]),
                        logical(1))
  if (any(is_list_col)) {
    cd <- cd[, !is_list_col, drop = FALSE]
  }

  # Drop the exact scRepertoire clonal-metadata columns when writing H5MU
  # (the AIRR side carries that info). Use an explicit list to avoid
  # collateral damage on user columns like CTLA4.
  if (isTRUE(drop_ct)) {
    ct_cols <- intersect(c("CTgene", "CTnt", "CTaa", "CTstrict"),
                         colnames(cd))
    if (length(ct_cols)) {
      cd <- cd[, setdiff(colnames(cd), ct_cols), drop = FALSE]
    }
  }

  # Subset by obs_columns if supplied.
  if (!is.null(obs_columns)) {
    unknown <- setdiff(obs_columns, colnames(cd))
    if (length(unknown)) {
      stop("obs_columns not found in colData: ",
           paste(unknown, collapse = ", "), call. = FALSE)
    }
    cd <- cd[, obs_columns, drop = FALSE]
  }

  SummarizedExperiment::colData(sce) <- cd
  sce
}

#' @keywords internal
# Map AIRR-canonical locus names to the chain strings immApex::getIR()
# accepts. immApex uses "Heavy" / "Light" rather than IGH / IGK / IGL,
# and returns light chains mixed (IGK + IGL); we disambiguate later by
# V-gene prefix when the caller asked for a specific light chain.
.airr_locus_to_immapex <- function(locus) {
  switch(locus,
         "TRA" = "TRA", "TRB" = "TRB",
         "TRG" = "TRG", "TRD" = "TRD",
         "IGH" = "Heavy",
         "IGK" = "Light", "IGL" = "Light",
         locus)
}

#' @keywords internal
.buildAIRR <- function(input, chains) {
  loci <- if (identical(chains, "both")) c("TRA", "TRB") else chains

  parts <- list()
  for (locus in loci) {
    immapex_locus <- .airr_locus_to_immapex(locus)

    # immApex returns the CDR3 sequence in a column named cdr3_aa
    # regardless of sequence.type; call twice to get both AA and NT.
    ir_aa <- tryCatch(
      immApex::getIR(input, chains = immapex_locus, sequence.type = "aa"),
      error = function(e) NULL
    )
    if (is.null(ir_aa) || nrow(ir_aa) == 0) next

    ir_nt <- tryCatch(
      immApex::getIR(input, chains = immapex_locus, sequence.type = "nt"),
      error = function(e) NULL
    )

    keep <- !is.na(ir_aa$cdr3_aa) & nzchar(ir_aa$cdr3_aa)

    # immApex pools IGK + IGL under "Light"; if the caller wanted one
    # specifically, filter rows by V-gene prefix.
    if (locus %in% c("IGK", "IGL")) {
      v_safe <- ifelse(is.na(ir_aa$v), "", ir_aa$v)
      keep <- keep & startsWith(v_safe, locus)
    }

    if (!any(keep)) next

    barcode <- ir_aa$barcode[keep]
    n <- length(barcode)
    junction_aa <- ir_aa$cdr3_aa[keep]
    v_call <- ir_aa$v[keep]
    j_call <- ir_aa$j[keep]

    junction <- if (!is.null(ir_nt) && "cdr3_aa" %in% colnames(ir_nt)) {
      ir_nt$cdr3_aa[keep]
    } else {
      rep(NA_character_, n)
    }

    # Heuristic productivity: a valid CDR3 has no stop codon (*) or
    # frameshift marker (_). scRepertoire upstream-filters most non-
    # productive contigs, so this is mostly TRUE in practice but is at
    # least derived from the data rather than fabricated.
    productive <- !grepl("\\*|_", junction_aa)

    parts[[locus]] <- data.frame(
      cell_id     = barcode,
      sequence_id = paste(barcode, locus, seq_len(n), sep = "_"),
      locus       = locus,
      productive  = productive,
      v_call      = v_call,
      j_call      = j_call,
      junction_aa = junction_aa,
      junction    = junction,
      stringsAsFactors = FALSE
    )
  }

  if (length(parts) == 0) return(NULL)
  do.call(rbind, c(parts, list(make.row.names = FALSE)))
}

#' Export Single-Cell + Immune Receptor Data to scanpy/scirpy Format
#'
#' @description Exports a \code{SingleCellExperiment} or \code{Seurat} object,
#'   optionally combined with scRepertoire immune receptor metadata, to
#'   scanpy/scirpy-compatible H5AD or H5MU files. Writes an AIRR
#'   rearrangement TSV sidecar by default when receptor data is present.
#'
#' @param input A \code{SingleCellExperiment} or \code{Seurat} object.
#' @param output_dir Directory to write outputs into. Created if missing.
#'   Files written with fixed names: \code{adata.h5ad},
#'   \code{airr_rearrangement.tsv.gz}, and (for H5MU) \code{mudata.h5mu}.
#' @param format Either \code{"h5ad"} (gene expression only) or
#'   \code{"h5mu"} (combined GEX + AIRR via scirpy/muon). Default \code{"h5ad"}.
#' @param write_airr Logical. Write the standalone AIRR TSV sidecar when
#'   immune receptor data is present. Default \code{TRUE}.
#' @param obs_columns Character vector of \code{colData} columns to keep.
#'   Default \code{NULL}: keep all non-list columns; additionally drop
#'   scRepertoire \code{CT*} columns when \code{format = "h5mu"}.
#' @param reductions Character vector of \code{reducedDimNames()} to write.
#'   Default \code{NULL}: write all.
#' @param assay Assay name to populate \code{adata.X}. Default \code{"counts"}.
#' @param chains Which receptor loci to extract: \code{"both"} (TRA + TRB),
#'   or any of \code{"TRA"}, \code{"TRB"}, \code{"TRG"}, \code{"TRD"},
#'   \code{"IGH"}, \code{"IGL"}, \code{"IGK"}.
#' @param overwrite Logical. Overwrite existing output files. Default \code{FALSE}.
#' @param verbose Logical. Emit progress messages. Default \code{TRUE}.
#'
#' @return Invisibly, a list with paths: \code{h5ad}, \code{airr}, \code{h5mu}
#'   (any of which may be \code{NULL} if not written).
#'
#' @export
#' @importFrom methods is
#' @importFrom SummarizedExperiment colData colData<- assayNames
#'
#' @examples
#' data(immLynx_example)
#'
#' # Inspect what would be written without invoking Python.
#' airr <- immLynx:::.buildAIRR(immLynx_example, chains = "TRB")
#' head(airr)
#'
#' \donttest{
#' # Gene expression only -> adata.h5ad (+ AIRR sidecar by default).
#' out_h5ad <- exportToScanpy(
#'   immLynx_example,
#'   output_dir = tempfile("scanpy_h5ad_"),
#'   format     = "h5ad"
#' )
#' file.exists(out_h5ad$h5ad)
#' file.exists(out_h5ad$airr)
#'
#' # Subset metadata and embeddings before writing.
#' out_lite <- exportToScanpy(
#'   immLynx_example,
#'   output_dir  = tempfile("scanpy_lite_"),
#'   format      = "h5ad",
#'   obs_columns = c("Patient", "Type", "clusters"),
#'   write_airr  = FALSE
#' )
#'
#' # Combined GEX + AIRR -> mudata.h5mu (requires the scanpy export
#' # basilisk env; the first call may take several minutes to install).
#' out_h5mu <- exportToScanpy(
#'   immLynx_example,
#'   output_dir = tempfile("scanpy_h5mu_"),
#'   format     = "h5mu",
#'   chains     = "TRB"
#' )
#' }
exportToScanpy <- function(input,
                           output_dir,
                           format      = c("h5ad", "h5mu"),
                           write_airr  = TRUE,
                           obs_columns = NULL,
                           reductions  = NULL,
                           assay       = "counts",
                           chains      = c("both", "TRA", "TRB",
                                           "TRG", "TRD",
                                           "IGH", "IGL", "IGK"),
                           overwrite   = FALSE,
                           verbose     = TRUE) {

  format <- match.arg(format)
  chains <- match.arg(chains)

  # --- Validate output_dir up front ------------------------------------
  if (missing(output_dir) || is.null(output_dir) ||
      !is.character(output_dir) || length(output_dir) != 1L ||
      is.na(output_dir) || !nzchar(output_dir)) {
    stop("output_dir must be a single non-empty path string",
         call. = FALSE)
  }
  if (file.exists(output_dir) && !dir.exists(output_dir)) {
    stop("output_dir exists but is not a directory: ", output_dir,
         call. = FALSE)
  }

  # --- Input dispatch ---------------------------------------------------
  sce <- .coerceToSCE(input)

  # --- Validate args ----------------------------------------------------
  if (!assay %in% SummarizedExperiment::assayNames(sce)) {
    stop("assay '", assay, "' not in assayNames(input). Available: ",
         paste(SummarizedExperiment::assayNames(sce), collapse = ", "),
         call. = FALSE)
  }
  rd_names <- SingleCellExperiment::reducedDimNames(sce)
  if (!is.null(reductions)) {
    if (!is.character(reductions) || length(reductions) == 0L) {
      stop("reductions must be NULL or a non-empty character vector",
           call. = FALSE)
    }
    missing_r <- setdiff(reductions, rd_names)
    if (length(missing_r)) {
      stop("reductions not found: ", paste(missing_r, collapse = ", "),
           ". Available: ", paste(rd_names, collapse = ", "),
           call. = FALSE)
    }
  }

  # --- Up-front dependency check for AIRR sidecar ----------------------
  # Catch missing data.table before any I/O so we don't leave behind a
  # half-written H5AD when the user can't satisfy the dependency.
  if (write_airr && !requireNamespace("data.table", quietly = TRUE)) {
    stop("write_airr = TRUE requires the data.table package; ",
         "install.packages('data.table') or pass write_airr = FALSE",
         call. = FALSE)
  }

  # --- AIRR up front (so we can refuse h5mu without TCR) ---------------
  airr_df <- .buildAIRR(sce, chains = chains)
  has_tcr <- !is.null(airr_df) && nrow(airr_df) > 0
  if (format == "h5mu" && !has_tcr) {
    stop("format = 'h5mu' requires immune receptor metadata; ",
         "found none. Use format = 'h5ad' or add scRepertoire data first.",
         call. = FALSE)
  }

  # --- Output dir + file paths ----------------------------------------
  if (!dir.exists(output_dir)) {
    ok <- dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    if (!ok || !dir.exists(output_dir)) {
      stop("Failed to create output_dir: ", output_dir, call. = FALSE)
    }
    if (verbose) message("Created ", output_dir)
  }
  h5ad_path <- file.path(output_dir, "adata.h5ad")
  airr_path <- file.path(output_dir, "airr_rearrangement.tsv.gz")
  h5mu_path <- file.path(output_dir, "mudata.h5mu")

  paths_to_check <- c(
    h5ad_path,
    if (write_airr && has_tcr) airr_path,
    if (format == "h5mu") h5mu_path
  )
  for (p in paths_to_check) {
    if (file.exists(p)) {
      if (!isTRUE(overwrite)) {
        stop("Output exists (use overwrite = TRUE to replace): ", p,
             call. = FALSE)
      }
      unlink(p)
    }
  }

  # --- Prune colData and reducedDims ----------------------------------
  drop_ct <- (format == "h5mu")
  sce_out <- .pruneColData(sce, obs_columns = obs_columns,
                           drop_ct = drop_ct)
  if (!is.null(reductions)) {
    drop_rd <- setdiff(rd_names, reductions)
    for (nm in drop_rd) {
      SingleCellExperiment::reducedDim(sce_out, nm) <- NULL
    }
  }

  # --- Write H5AD ------------------------------------------------------
  if (verbose) message("Writing H5AD: ", h5ad_path)
  zellkonverter::writeH5AD(sce_out, file = h5ad_path, X_name = assay)

  result <- list(h5ad = h5ad_path, airr = NULL, h5mu = NULL)

  # --- AIRR sidecar ----------------------------------------------------
  if (write_airr) {
    if (has_tcr) {
      if (verbose) message("Writing AIRR TSV: ", airr_path)
      # fwrite defaults are AIRR-compliant: na = "", quote = "auto",
      # eol = "\n", bom = FALSE. Do not change without re-checking
      # scirpy.io.read_airr behavior.
      data.table::fwrite(airr_df, airr_path, sep = "\t",
                         compress = "gzip")
      result$airr <- airr_path
    } else {
      warning("write_airr = TRUE but no immune receptor data found; ",
              "skipping AIRR sidecar (H5AD was written).",
              call. = FALSE)
    }
  }

  # --- H5MU (scirpy + muon) --------------------------------------------
  if (format == "h5mu") {
    if (verbose) message("Writing H5MU: ", h5mu_path)

    # scirpy.io.read_airr needs a path on disk. If the user asked for
    # write_airr = FALSE, drop a temp gzipped TSV ourselves.
    airr_for_scirpy <- if (!is.null(result$airr)) {
      result$airr
    } else {
      tmp <- tempfile(fileext = ".tsv.gz")
      data.table::fwrite(airr_df, tmp, sep = "\t", compress = "gzip")
      on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
      tmp
    }

    # Known scirpy quirk: gzipped AIRR can fail; decompress to a temp
    # uncompressed TSV in R first, then hand the plain path to scirpy.
    tmp_tsv <- tempfile(fileext = ".tsv")
    on.exit(if (file.exists(tmp_tsv)) unlink(tmp_tsv), add = TRUE)
    con_in <- gzfile(airr_for_scirpy, "rb")
    on.exit(try(close(con_in), silent = TRUE), add = TRUE)
    con_out <- file(tmp_tsv, "wb")
    on.exit(try(close(con_out), silent = TRUE), add = TRUE)
    repeat {
      buf <- readBin(con_in, "raw", n = 1024L * 1024L)
      if (length(buf) == 0L) break
      writeBin(buf, con_out)
    }
    close(con_in)
    close(con_out)

    # All Python imports below MUST use convert = FALSE so MuData and
    # AnnData stay as Python proxies; auto-conversion silently corrupts
    # column dtypes and breaks the chained $-method calls.
    h5mu_result <- tryCatch(
      .run_in_scanpy_env(function(h5ad, airr_tsv, out) {
        ad <- reticulate::import("anndata", convert = FALSE)
        ir <- reticulate::import("scirpy", convert = FALSE)
        mu <- reticulate::import("muon",   convert = FALSE)

        adata_gex  <- ad$read_h5ad(h5ad)
        adata_airr <- ir$io$read_airr(airr_tsv)

        mdata <- mu$MuData(reticulate::dict(gex = adata_gex,
                                            airr = adata_airr))
        ir$pp$index_chains(mdata)
        ir$tl$chain_qc(mdata)
        mdata$write(out, compression = "gzip")
        invisible(NULL)
      },
      h5ad = h5ad_path, airr_tsv = tmp_tsv, out = h5mu_path),
      error = function(e) e
    )

    if (inherits(h5mu_result, "error")) {
      partial <- c(if (file.exists(h5ad_path)) h5ad_path,
                   if (!is.null(result$airr) &&
                       file.exists(result$airr)) result$airr)
      stop("H5MU assembly failed inside scirpy/muon: ",
           conditionMessage(h5mu_result), "\n",
           "Partial outputs left on disk:\n  ",
           paste(partial, collapse = "\n  "), "\n",
           "Re-run with overwrite = TRUE after addressing the cause.",
           call. = FALSE)
    }

    result$h5mu <- h5mu_path
  }

  invisible(result)
}

#' @keywords internal
.coerceToSCE <- function(input) {
  if (methods::is(input, "SingleCellExperiment")) return(input)
  if (methods::is(input, "Seurat")) {
    if (!requireNamespace("Seurat", quietly = TRUE)) {
      stop("Seurat input requires the Seurat package. ",
           "Install with: install.packages('Seurat')", call. = FALSE)
    }
    # Seurat emits one warning per empty layer ("Layer 'data' is empty") when
    # converting an object that only carries counts. That is the normal state
    # for a freshly created object and is not actionable here, so muffle just
    # those and let every other warning through.
    return(withCallingHandlers(
      Seurat::as.SingleCellExperiment(input),
      warning = function(w) {
        if (grepl("^Layer '.*' is empty", conditionMessage(w))) {
          invokeRestart("muffleWarning")
        }
      }
    ))
  }
  stop("input must be a SingleCellExperiment or Seurat object",
       call. = FALSE)
}
