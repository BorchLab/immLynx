# Bioconductor-compliant management of Python environments
#
# This script defines basilisk environments for immLynx.
#
# Three separate environments are used so that installation failures in
# heavy or problematic packages do not prevent the rest of the Python
# tooling from working:
#
#   immLynxEnv       -- lightweight, stable dependencies: OLGA, tcrdist3,
#                       metaclonotypist, pyrepseq.  Used by the majority
#                       of the package functions.
#
#   immLynxTorchEnv  -- torch + transformers for ESM-2 embeddings.
#                       Separated because torch is a very large download
#                       (~850 MB) and can cause timeouts or disk pressure
#                       in CI environments.
#
#   immLynxExtraEnv  -- packages that have historically had installation
#                       issues: soNNia and clusTCR.  Used only by
#                       calculate.sonia() and calculate.clustcr().
#
# Note: Version checking is disabled because clusTCR is installed directly
# from its GitHub repository via a git+https:// URL in the `pip` vector,
# which does not use the standard `==` version specifier that basilisk
# expects.  The `paths` parameter cannot be used for git URLs because
# basilisk prepends the package system directory to each path entry.

basilisk::setBasiliskCheckVersions(FALSE)

# ---------------------------------------------------------------------------
# Core environment -- lightweight, stable packages
# ---------------------------------------------------------------------------
immLynxEnv <- basilisk::BasiliskEnvironment(
    envname = "immLynxEnv",
    pkgname = "immLynx",
    packages = c(
        "python=3.9",
        "numpy=1.23.4",
        "scipy=1.8.0",
        "pandas=1.4.4",
        "scikit-learn=1.1.3"
    ),
    pip = c(
        "tcrdist3==0.2.2",
        "olga==1.2.4",
        "metaclonotypist==0.2.0",
        "pyrepseq==1.5.1"
    )
)

# ---------------------------------------------------------------------------
# Torch environment -- torch + transformers for ESM-2 embeddings
# ---------------------------------------------------------------------------
immLynxTorchEnv <- basilisk::BasiliskEnvironment(
    envname = "immLynxTorchEnv",
    pkgname = "immLynx",
    packages = c(
        "python=3.9",
        "numpy=1.23.4"
    ),
    pip = c(
        "torch==2.1.2",
        "transformers==4.36.2"
    )
)

# ---------------------------------------------------------------------------
# Extra environment -- soNNia + clusTCR (may fail to install)
# ---------------------------------------------------------------------------
immLynxExtraEnv <- basilisk::BasiliskEnvironment(
    envname = "immLynxExtraEnv",
    pkgname = "immLynx",
    packages = c(
        "python=3.9",
        "numpy=1.23.4",
        "scipy=1.8.0",
        "pandas=1.4.4",
        "scikit-learn=1.1.3",
        "markov-clustering=0.0.6.dev0"
    ),
    pip = c(
        "olga==1.2.4",
        "sonnia==0.2.5",
        "git+https://github.com/svalkiers/clusTCR.git@1.0.3"
    )
)
