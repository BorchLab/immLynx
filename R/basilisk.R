# Bioconductor-compliant management of Python environments
#
# This script defines basilisk environments for immLynx.
#
# Two separate environments are used so that installation failures in
# problematic packages (soNNia and clusTCR) do not prevent the rest of
# the Python tooling from working.
#
#   immLynxEnv      -- stable dependencies used by the majority of the
#                      package: OLGA, tcrdist3, metaclonotypist,
#                      transformers/torch (ESM-2 embeddings), etc.
#
#   immLynxExtraEnv -- packages that have historically had installation
#                      issues: soNNia and clusTCR.  Used only by
#                      calculate.sonia() and calculate.clustcr().
#
# Note: Version checking is disabled because clusTCR is installed directly
# from its GitHub repository via a git+https:// URL in the `pip` vector,
# which does not use the standard `==` version specifier that basilisk
# expects.  The `paths` parameter cannot be used for git URLs because
# basilisk prepends the package system directory to each path entry.

basilisk::setBasiliskCheckVersions(FALSE)

# ---------------------------------------------------------------------------
# Core environment -- stable packages
# ---------------------------------------------------------------------------
immLynxEnv <- basilisk::BasiliskEnvironment(
    envname = "immLynxEnv",
    pkgname = "immLynx",
    packages = c(
        "python=3.9",
        "numpy=1.23.4",
        "scipy=1.8.0",
        "pandas=1.4.4",
        "matplotlib=3.5.3",
        "scikit-learn=1.1.3",
        "statsmodels=0.13.2",
        "seaborn=0.12.1",
        "faiss-cpu=1.7.4"
    ),
    pip = c(
        "tcrdist3==0.2.2",
        "olga==1.2.4",
        "metaclonotypist==0.2.0",
        "pyrepseq==1.5.1",
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
