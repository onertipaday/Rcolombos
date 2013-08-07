\name{listGenes}
\alias{listGenes}
\title{This method takes as parameter a single string, representing an organism,
and returns a character vector corresponding to the currently available organisms.}
\usage{
  listGenes(organism = "ecoli")
}
\arguments{
  \item{organism}{A string containing the organism of
  interest.}
}
\value{
  A data.frame containing the locustag and description of
  all the genes for the selected organism.
}
\description{
  This method takes as parameter a single string,
  representing an organism, and returns a character vector
  corresponding to the currently available organisms.
}
\examples{
\dontrun{
library('Rcolombos')
listGenes()
}
}
\references{
  http://colombos.net
}

