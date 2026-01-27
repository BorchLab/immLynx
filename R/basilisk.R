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
        "pandas=2.2.1",
        "matplotlib=3.8.0",
        "bioconda::markov_clustering=0.0.6",
        "svalkiers::clustcr=1.0.3",
        "pytorch::faiss-cpu=1.7.4",
        "conda-forge::scikit-learn=1.4.0",
        "conda-forge::numpy=1.26.4",
        "conda-forge::scipy=1.12.0",
        "conda-forge::statsmodels=0.14.1",
        "conda-forge::seaborn=0.13.2"
    ),
    channels = c("conda-forge", "bioconda", "svalkiers", "pytorch"),
    pip = c(
        "tcrdist3==0.2.2",
        "DeepTCR==2.1.1",
        "olga==1.2.4",
        "sonnia==0.1.0",
        "metaclonotypist==0.2.0",
        "pyrepseq==1.5.1",
        "tensorflow==2.15.0"
    )
)
