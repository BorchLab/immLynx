# Regression guards for basilisk environment specifications.
#
# Named "zz" so testthat runs these last, after any test that may have
# triggered a first-use environment build.
#
# Background: basilisk::setupBasiliskEnv() hands the `pip` vector straight to
# system2(), which pastes the arguments into a single string and runs it via
# /bin/sh without quoting. A spec such as "anndata>=0.8" is therefore split by
# the shell into the word "anndata", a ">" redirect, and the target "=0.8".
# Two things go wrong at once: the version floor is silently dropped, and pip's
# stdout is written to a file literally named "=0.8" in the working directory.
#
# The conda `packages` vector takes a different route: it reaches
# reticulate::conda_install(), which applies maybe_shQuote() to every argument
# before calling system2(). Comparison operators are therefore safe there, and
# a floor such as "pytorch-cpu>=2.5" is deliberate — pinning conda specs
# exactly makes the solve brittle across platforms for no safety gain. So the
# two vectors get two different rules, matching the two code paths.

# pip: characters that survive an unquoted /bin/sh word split. Everything in
# the current specs is covered — version pins ("==", "="), dotted versions,
# dashes, underscores, local-version "+", and the git+https://...@tag URL used
# for clusTCR (":", "/", "@", "+").
SHELL_SAFE <- "^[A-Za-z0-9._+/@=:-]+$"

# conda: quoting makes redirects harmless, but a spec should still never carry
# whitespace, quotes, or command-substitution characters. Those signal a typo
# ("numpy >= 1.2") or an injection rather than a version bound.
# Note "<" and ">" are intentionally absent: they are valid, quoted, and
# meaningful in a conda version bound.
CONDA_JUNK <- "[[:space:]\"'`$;&|()]"

# Collect every BasiliskEnvironment object defined in the package namespace,
# so environments added later are covered without editing this test.
basilisk_envs <- function() {
  ns <- asNamespace("immLynx")
  found <- list()
  for (nm in ls(ns, all.names = TRUE)) {
    obj <- tryCatch(get(nm, envir = ns), error = function(e) NULL)
    if (methods::is(obj, "BasiliskEnvironment")) {
      found[[nm]] <- obj
    }
  }
  found
}

test_that("every basilisk environment is discoverable for inspection", {
  envs <- basilisk_envs()

  # Guard the guard: if the collector silently returns nothing, the spec
  # checks below would pass vacuously.
  expect_gt(length(envs), 0)
  expect_true("immLynxEnv" %in% names(envs))
  expect_true("scanpyExportEnv" %in% names(envs))
})

test_that("pip specs contain no shell metacharacters", {
  envs <- basilisk_envs()

  for (nm in names(envs)) {
    specs <- envs[[nm]]@pip
    if (!length(specs)) next

    bad <- specs[!grepl(SHELL_SAFE, specs)]
    expect_identical(
      bad, character(0),
      info = paste0(
        "Unsafe pip spec(s) in ", nm, ": ",
        paste(sprintf("'%s'", bad), collapse = ", "),
        ". basilisk passes the pip vector through an unquoted shell, so '>' ",
        "or '<' becomes a redirect: the version bound is dropped and a stray ",
        "file is created. Use an exact '==' pin instead."
      )
    )
  }
})

test_that("conda package specs are well formed", {
  envs <- basilisk_envs()

  for (nm in names(envs)) {
    specs <- envs[[nm]]@packages
    if (!length(specs)) next

    bad <- specs[grepl(CONDA_JUNK, specs)]
    expect_identical(
      bad, character(0),
      info = paste0(
        "Malformed conda spec(s) in ", nm, ": ",
        paste(sprintf("'%s'", bad), collapse = ", "),
        ". Version bounds such as 'pytorch-cpu>=2.5' are fine here, but a ",
        "spec must not contain whitespace, quotes, or shell command ",
        "characters."
      )
    )
  }
})

test_that("no shell-redirect artefact files were left in the package tree", {
  # An unquoted ">" redirect writes to the working directory of the R process
  # that built the environment. Under devtools::test() that is tests/testthat;
  # under R CMD check it is the corresponding directory inside .Rcheck. Walking
  # up two levels covers both, plus the package root for interactive use.
  here <- normalizePath(getwd(), mustWork = FALSE)
  roots <- unique(c(
    here,
    normalizePath(file.path(here, ".."), mustWork = FALSE),
    normalizePath(file.path(here, "..", ".."), mustWork = FALSE)
  ))
  roots <- roots[dir.exists(roots)]

  stray <- unlist(lapply(roots, function(r) {
    list.files(r, pattern = "^=", all.files = TRUE, full.names = TRUE)
  }), use.names = FALSE)
  stray <- if (is.null(stray)) character(0) else stray

  expect_identical(
    stray, character(0),
    info = paste0(
      "Found file(s) whose name begins with '=': ",
      paste(stray, collapse = ", "),
      ". These are shell redirect artefacts from an unquoted version spec. ",
      "Delete them and fix the offending spec in R/basilisk.R."
    )
  )
})
