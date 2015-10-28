.onAttach <- function(libname, pkgname) {
    packageStartupMessage("\nRcolombos version ", utils::packageVersion("Rcolombos"), ", ?quick_search to start.")
    # options("REST.version"="http://rest.colombos.fmach.it/")
    options("REST.version"="http://luca.colombos-dev.fmach.it/")
}