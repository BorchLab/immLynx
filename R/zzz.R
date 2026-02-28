.onLoad <- function(libname, pkgname) {
    # Disable basilisk version checking so that clusTCR can be installed from
    # its GitHub repository using a git+https:// URL in the pip vector.
    # This must happen before any basilisk environment is created.
    basilisk::setBasiliskCheckVersions(FALSE)
}
