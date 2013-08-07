\name{listContrasts}
\alias{listContrasts}
\title{This method takes as parameter a single string, representing an organism,
and returns a character vector corresponding to the currently available organisms.}
\usage{
  listContrasts(organism = "ecoli")
}
\arguments{
  \item{organism}{A string containing the organism of
  interest.}
}
\value{
  A data.frame containing the contrasts and GSM of all the
  contrasts for the selected organism.
}
\description{
  This method takes as parameter a single string,
  representing an organism, and returns a character vector
  corresponding to the currently available organisms.
}
\examples{
\dontrun{
library('Rcolombos')
listContrasts()
}
}
\references{
  http://colombos.net
}

