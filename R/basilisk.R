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
#
# Note: Version checking is disabled because clusTCR is installed directly
# from its GitHub repository via a git+https:// URL in the `pip` vector,
# which does not use the standard `==` version specifier that basilisk
# expects.  The `paths` parameter cannot be used for git URLs because
# basilisk prepends the package system directory to each path entry.

basilisk::setBasiliskCheckVersions(FALSE)

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
        "markov-clustering=0.0.6.dev0",
        "faiss-cpu=1.7.4"
    ),
    pip = c(
        "tcrdist3==0.2.2",
        "olga==1.2.4",
        "sonnia==0.1.0",
        "metaclonotypist==0.2.0",
        "pyrepseq==1.5.1",
        "torch==2.1.2",
        "transformers==4.36.2",
        "git+https://github.com/svalkiers/clusTCR.git@1.0.3"
    )
)
