# Bioconductor-compliant management of Python environments
#
# This script defines a basilisk environment for immlynx.
#
# The environment is defined using the BasiliskEnvironment function, which
# specifies the name of the environment, the package name, and the list of
# required Python packages.
#
# The defined environment is then used in the wrapper functions to execute
# Python code using basiliskRun.

immLynxEnv <- basilisk::BasiliskEnvironment(
    envname = "immLynxEnv",
    pkgname = "immLynx",
    packages = c(
        "python=3.9",
        "pandas>=2.2.1",
        "tensorflow",
        "matplotlib",
        "svalkiers::clustcr",
        "pytorch::faiss-cpu",
        "conda-forge::markov_clustering",
        "conda-forge::scikit-learn",
        "conda-forge::numpy>=1.26.4",
        "conda-forge::scipy>=1.12.0",
        "conda-forge::statsmodels>=0.14.1",
        "conda-forge::seaborn"
    ),
    pip = c(
        "tcrdist3",
        "DeepTCR",
        "olga",
        "sonnia",
        "metaclonotypist",
        "pyrepseq>=1.5.1"
    )
)
