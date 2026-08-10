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

# Separate environment for the scanpy/scirpy export workflow.
# Kept apart from immLynxEnv to avoid version conflicts with tcrdist3,
# olga, sonnia, and clusTCR (which pin older numpy/pandas/scipy).
scanpyExportEnv <- basilisk::BasiliskEnvironment(
    envname = "scanpyExportEnv",
    pkgname = "immLynx",
    packages = c(
        "python=3.10",
        "numpy",
        "pandas",
        "scipy",
        "h5py"
    ),
    pip = c(
        # Pinned exactly rather than ">=0.8": basilisk hands the pip vector to
        # system2(), which runs it through an unquoted shell, so ">=0.8" is
        # parsed as a redirect to a file named "=0.8".  That both dropped the
        # version floor and littered the working directory.  0.11.4 is the
        # version pip already resolved to under python=3.10.
        "anndata==0.11.4",
        "scanpy",
        "muon",
        "scirpy"
    )
)

# Separate environment for scXpand (clonal expansion prediction from gene
# expression).  Cannot share scanpyExportEnv (python 3.10) or immLynxEnv
# (python 3.9, torch 2.1.2, numpy 1.23): scxpand requires python >= 3.11
# and torch >= 2.5.
#
# pytorch-cpu is pulled from conda-forge rather than letting pip resolve
# `torch`, because the default PyPI torch wheel on linux-x86_64 bundles
# CUDA and runs to several gigabytes.  basilisk's `pip` entries are bare
# specifiers, so there is no way to inject
# --index-url https://download.pytorch.org/whl/cpu.  Installing the CPU
# build through conda first satisfies scxpand's torch>=2.5 requirement.
#
# Only scxpand itself is listed under `pip`; it resolves the rest of the
# stack (scanpy, anndata, scirpy, lightgbm, optuna, pooch, shap, ...).
# Enumerating those here would only create version-pin drift.
scXpandEnv <- basilisk::BasiliskEnvironment(
    envname = "scXpandEnv",
    pkgname = "immLynx",
    packages = c(
        "python=3.11",
        "pytorch-cpu>=2.5",
        "numpy",
        "pandas",
        "scipy",
        "h5py"
    ),
    pip = c(
        "scxpand==0.4.6"
    )
)

# Separate environment for DeepTCR. DeepTCR pins its whole scientific stack
# (numpy 1.23.5, pandas 1.5.3, scipy 1.10.1, TensorFlow 2.12) with "==", so it
# cannot share immLynxEnv.
#
# Three constraints drove this layout, each verified by building the
# environment and training a VAE:
#
#   1. python=3.10, not 3.11.  DeepTCR pins biopython==1.76, which ships no
#      wheel past cp38 and none for macOS arm64, so it compiles from source.
#      Its C extension assigns to Py_TYPE(), which CPython 3.11 made a hard
#      error.  It compiles cleanly on 3.10.
#
#   2. The stack comes from conda rather than pip so that TensorFlow arrives
#      without Apple's metal plugin.  On macOS, DeepTCR's requirements pull
#      tensorflow-metal==0.8.0, which is built against TF 2.11 and aborts the
#      process at import under TF 2.12 with "platform is already registered
#      with name: METAL".  conda-forge's tensorflow has no such plugin.
#
#   3. biopython stays on pip.  conda-forge's oldest osx-arm64 build is 1.78,
#      which removed Bio.Alphabet, and DeepTCR still imports it.
#
# "--no-deps" therefore applies to the whole pip step: every dependency is
# already satisfied by conda, and it also keeps DeepTCR's jupyterlab and
# notebook requirements out of the environment.
deepTCREnv <- basilisk::BasiliskEnvironment(
    envname = "deepTCREnv",
    pkgname = "immLynx",
    packages = c(
        "python=3.10",
        "tensorflow=2.12",
        "numpy=1.23.5",
        "pandas=1.5.3",
        "scipy=1.10.1",
        "h5py=3.8.0",
        "scikit-learn=1.2.2",
        "matplotlib-base=3.7.2",
        "seaborn=0.12.2",
        "umap-learn",
        "networkx",
        "tqdm",
        "psutil"
    ),
    pip = c(
        "--no-deps",
        "DeepTCR==2.1.29",
        "biopython==1.76",
        "logomaker==0.8",
        "distinctipy==1.2.1",
        "python-louvain==0.16"
    )
)
