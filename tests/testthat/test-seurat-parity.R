# Every immLynx wrapper must accept a SingleCellExperiment or a Seurat object
# and return the class it was given.
#
# The extraction layer is already class-agnostic because immApex::getIR()
# handles both, so these tests target the guards and the write-back paths,
# which are where the two diverge.

sce_fixture <- function() {
  data("immLynx_example", package = "immLynx")
  immLynx_example
}

seurat_fixture <- function() {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("Matrix")
  sce <- sce_fixture()
  cts <- SummarizedExperiment::assay(sce, "counts")
  seu <- Seurat::CreateSeuratObject(
    counts = Matrix::Matrix(as.matrix(cts), sparse = TRUE)
  )
  md <- as.data.frame(SummarizedExperiment::colData(sce))
  for (nm in colnames(md)) seu[[nm]] <- md[[nm]]
  seu
}

# Read a per-cell column back out regardless of container.
cell_col <- function(obj, nm) {
  if (methods::is(obj, "SingleCellExperiment")) {
    SummarizedExperiment::colData(obj)[[nm]]
  } else {
    obj[[nm]][, 1]
  }
}

# ===========================================================================
# Guards: no wrapper may reject a Seurat object outright
# ===========================================================================

test_that("no wrapper reports the old SingleCellExperiment-only error", {
  # The historical guard read "Input must be a SingleCellExperiment object"
  # and turned Seurat users away. Assert on behaviour rather than on source
  # text, so this holds under R CMD check where the sources are not present.
  wrappers <- list(runClustTCR, runOLGA, runEmbeddings, runTCRdist,
                   runSymdelNeighbors, runDeepTCR, runSoNNia)

  for (fn in wrappers) {
    msg <- tryCatch(fn(data.frame(x = 1)),
                    error = function(e) conditionMessage(e))
    expect_false(grepl("must be a SingleCellExperiment object", msg,
                       fixed = TRUE))
    # It must still name Seurat as an accepted class.
    expect_match(msg, "Seurat")
  }
})

test_that("wrappers still reject objects that are neither class", {
  for (fn in list(runClustTCR, runOLGA, runEmbeddings, runTCRdist,
                  runSymdelNeighbors, runDeepTCR, runSoNNia)) {
    expect_error(fn(data.frame(x = 1)),
                 "SingleCellExperiment or Seurat")
  }
})

test_that("the shared guard accepts a Seurat object", {
  expect_silent(immLynx:::.assertSCObject(seurat_fixture()))
})

# ===========================================================================
# Round trip: class in, same class out
# ===========================================================================

test_that("runSymdelNeighbors round-trips both classes", {
  skip_if_no_python()

  sce <- runSymdelNeighbors(sce_fixture(), chains = "TRB", max_edits = 1)
  seu <- runSymdelNeighbors(seurat_fixture(), chains = "TRB", max_edits = 1)

  expect_s4_class(sce, "SingleCellExperiment")
  expect_s4_class(seu, "Seurat")

  a <- cell_col(sce, "symdel_degree")
  b <- cell_col(seu, "symdel_degree")
  expect_equal(a, b)
})

test_that("runClustTCR round-trips both classes", {
  skip_if_no_python()

  sce <- runClustTCR(sce_fixture(), chains = "TRB")
  seu <- runClustTCR(seurat_fixture(), chains = "TRB")

  expect_s4_class(sce, "SingleCellExperiment")
  expect_s4_class(seu, "Seurat")
  expect_true("clustcr_TRB" %in% colnames(seu[[]]))
})

test_that("runOLGA round-trips both classes", {
  skip_if_no_python()

  sce <- runOLGA(sce_fixture(), chains = "TRB")
  seu <- runOLGA(seurat_fixture(), chains = "TRB")

  expect_s4_class(sce, "SingleCellExperiment")
  expect_s4_class(seu, "Seurat")
  expect_equal(cell_col(sce, "olga_pgen_TRB"),
               cell_col(seu, "olga_pgen_TRB"))
})

test_that("runEmbeddings writes a usable reduction on both classes", {
  skip_if_no_python()

  sce <- runEmbeddings(sce_fixture(), chains = "TRB")
  seu <- runEmbeddings(seurat_fixture(), chains = "TRB")

  expect_s4_class(sce, "SingleCellExperiment")
  expect_s4_class(seu, "Seurat")

  expect_true("tcr_esm" %in% SingleCellExperiment::reducedDimNames(sce))
  expect_true("tcr_esm" %in% Seurat::Reductions(seu))

  rd <- SingleCellExperiment::reducedDim(sce, "tcr_esm")
  emb <- Seurat::Embeddings(seu, "tcr_esm")
  expect_equal(dim(rd), dim(emb))
  expect_identical(rownames(emb), colnames(seu))
})

test_that("runTCRdist stores results on both classes", {
  skip_if_no_tcrdist()

  seu <- runTCRdist(seurat_fixture(), chains = "beta", add_to_object = TRUE)

  expect_s4_class(seu, "Seurat")
  expect_true("tcrdist" %in% names(methods::slot(seu, "misc")))
})

test_that("runDeepTCR writes a reduction on Seurat", {
  skip_if_no_deeptcr()

  seu <- runDeepTCR(seurat_fixture(), chains = "TRB",
                    latent_dim = 16, verbose = FALSE)

  expect_s4_class(seu, "Seurat")
  expect_true("tcr_deeptcr" %in% Seurat::Reductions(seu))
})
