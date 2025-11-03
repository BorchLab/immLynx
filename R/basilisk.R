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
        "python=3.8",
        "pandas",
        "tensorflow",
        "matplotlib",
        "svalkiers::clustcr",
        "pytorch::faiss-cpu",
        "conda-forge::markov_clustering",
        "conda-forge::scikit-learn"
    ),
    pip = c(
        "tcrdist3",
        "DeepTCR",
        "olga",
        "sonnia"
    )
)
