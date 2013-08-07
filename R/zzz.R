.onAttach <- function(libname, pkgname) {
  packageStartupMessage("\nRcolombos version ", utils::packageVersion("Rcolombos"), ", ?quick_search to start.")
  ## address bug in Windows RCurl
  if(.Platform$OS.type == "windows")
    options(RCurlOptions = list(cainfo = system.file("CurlSSL", "cacert.pem", package = "RCurl"), verbose = TRUE))
}