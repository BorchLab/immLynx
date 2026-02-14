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
        "conda-forge::pandas=1.4.4",
        "conda-forge::matplotlib=3.5.3",
        "bioconda::markov_clustering=0.0.6",
        "svalkiers::clustcr=1.0.3",
        "pytorch::faiss-cpu=1.7.4",
        "conda-forge::scikit-learn=1.1.3",
        "conda-forge::numpy=1.23.4",
        "conda-forge::scipy=1.8.0",
        "conda-forge::statsmodels=0.13.2",
        "conda-forge::seaborn=0.12.1"
    ),
    channels = c("conda-forge", "bioconda", "svalkiers", "pytorch"),
    pip = c(
        "tcrdist3==0.2.2",
        "olga==1.2.4",
        "sonnia==0.1.0",
        "metaclonotypist==0.2.0",
        "pyrepseq==1.5.1"
    )
)
