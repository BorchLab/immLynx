# Tests for the shared single-cell object accessors.
#
# immLynx functions accept either a SingleCellExperiment or a Seurat object
# and return the same class they were given. These helpers are the only place
# that knows the difference, so they carry the burden of that contract.

make_sce <- function() {
  data("immLynx_example", package = "immLynx")
  immLynx_example
}

make_seurat <- function() {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("Matrix")
  sce <- make_sce()
  cts <- SummarizedExperiment::assay(sce, "counts")
  seu <- Seurat::CreateSeuratObject(
    counts = Matrix::Matrix(as.matrix(cts), sparse = TRUE)
  )
  md <- as.data.frame(SummarizedExperiment::colData(sce))
  for (nm in colnames(md)) seu[[nm]] <- md[[nm]]
  seu
}

# ===========================================================================
# .isSCObject / .assertSCObject
# ===========================================================================

test_that(".isSCObject recognizes both supported classes", {
  expect_true(immLynx:::.isSCObject(make_sce()))
  expect_true(immLynx:::.isSCObject(make_seurat()))
})

test_that(".isSCObject rejects everything else", {
  expect_false(immLynx:::.isSCObject(data.frame(x = 1)))
  expect_false(immLynx:::.isSCObject(matrix(1)))
  expect_false(immLynx:::.isSCObject(NULL))
})

test_that(".assertSCObject names both accepted classes in its error", {
  err <- tryCatch(immLynx:::.assertSCObject(data.frame(x = 1)),
                  error = function(e) conditionMessage(e))
  expect_match(err, "SingleCellExperiment")
  expect_match(err, "Seurat")
})

test_that(".assertSCObject passes valid objects through silently", {
  expect_silent(immLynx:::.assertSCObject(make_sce()))
  expect_silent(immLynx:::.assertSCObject(make_seurat()))
})

# ===========================================================================
# .writeCellColumn
# ===========================================================================

test_that(".writeCellColumn writes to colData for SingleCellExperiment", {
  sce <- make_sce()
  cells <- colnames(sce)[1:5]

  out <- immLynx:::.writeCellColumn(sce, "myscore", seq_len(5), cells)

  expect_true("myscore" %in% colnames(SummarizedExperiment::colData(out)))
  vals <- SummarizedExperiment::colData(out)$myscore
  expect_identical(length(vals), ncol(sce))
  expect_identical(vals[1:5], seq_len(5))
  # Cells outside the supplied set stay NA.
  expect_true(all(is.na(vals[6:length(vals)])))
})

test_that(".writeCellColumn writes to meta.data for Seurat", {
  seu <- make_seurat()
  cells <- colnames(seu)[1:5]

  out <- immLynx:::.writeCellColumn(seu, "myscore", seq_len(5), cells)

  expect_true("myscore" %in% colnames(out[[]]))
  vals <- out[["myscore"]][, 1]
  # Seurat's ncol() returns a double where SingleCellExperiment returns an
  # integer, so compare by value rather than by type.
  expect_equal(length(vals), ncol(seu))
  expect_identical(vals[1:5], seq_len(5))
  expect_true(all(is.na(vals[6:length(vals)])))
})

test_that(".writeCellColumn returns the class it was given", {
  sce <- make_sce()
  seu <- make_seurat()

  expect_s4_class(
    immLynx:::.writeCellColumn(sce, "x", 1, colnames(sce)[1]),
    "SingleCellExperiment"
  )
  expect_s4_class(
    immLynx:::.writeCellColumn(seu, "x", 1, colnames(seu)[1]),
    "Seurat"
  )
})

# ===========================================================================
# .writeReduction
# ===========================================================================

test_that(".writeReduction stores a reduction on SingleCellExperiment", {
  sce <- make_sce()
  m <- matrix(rnorm(ncol(sce) * 4), nrow = ncol(sce), ncol = 4,
              dimnames = list(colnames(sce), NULL))

  out <- immLynx:::.writeReduction(sce, "myred", m, "MY_")

  expect_true("myred" %in% SingleCellExperiment::reducedDimNames(out))
  rd <- SingleCellExperiment::reducedDim(out, "myred")
  expect_identical(dim(rd), c(ncol(sce), 4L))
})

test_that(".writeReduction stores a reduction on Seurat", {
  seu <- make_seurat()
  m <- matrix(rnorm(ncol(seu) * 4), nrow = ncol(seu), ncol = 4,
              dimnames = list(colnames(seu), NULL))

  out <- immLynx:::.writeReduction(seu, "myred", m, "MY_")

  expect_true("myred" %in% Seurat::Reductions(out))
  emb <- Seurat::Embeddings(out, "myred")
  expect_equal(dim(emb), c(ncol(seu), 4L))
  expect_identical(rownames(emb), colnames(seu))
})

test_that(".writeReduction round-trips values unchanged", {
  sce <- make_sce()
  seu <- make_seurat()
  m <- matrix(seq_len(ncol(sce) * 3) * 1.0, nrow = ncol(sce), ncol = 3,
              dimnames = list(colnames(sce), NULL))

  rd <- SingleCellExperiment::reducedDim(
    immLynx:::.writeReduction(sce, "r", m, "R_"), "r")
  emb <- Seurat::Embeddings(
    immLynx:::.writeReduction(seu, "r", m, "R_"), "r")

  expect_equal(unname(rd), unname(m))
  expect_equal(unname(emb), unname(m))
})

test_that(".writeReduction tolerates rows of NA for cells without sequences", {
  sce <- make_sce()
  m <- matrix(NA_real_, nrow = ncol(sce), ncol = 2,
              dimnames = list(colnames(sce), NULL))
  m[1:10, ] <- 1

  expect_silent(immLynx:::.writeReduction(sce, "r", m, "R_"))
})

# ===========================================================================
# .writeObjMetadata
# ===========================================================================

test_that(".writeObjMetadata stores under metadata() for SingleCellExperiment", {
  sce <- make_sce()
  out <- immLynx:::.writeObjMetadata(sce, "mykey", list(a = 1))

  expect_identical(S4Vectors::metadata(out)$mykey, list(a = 1))
})

test_that(".writeObjMetadata stores under misc for Seurat", {
  seu <- make_seurat()
  out <- immLynx:::.writeObjMetadata(seu, "mykey", list(a = 1))

  expect_identical(methods::slot(out, "misc")$mykey, list(a = 1))
})
