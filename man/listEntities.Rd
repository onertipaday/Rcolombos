\name{listEntities}
\alias{listEntities}
\title{This method takes a string containing the nickname for the selected organism and a vector of string,
representing the genes of interest for the specified organism,
and returns a list containing the locustags, gene_names,
contrasts and M-values for the current selection.}
\usage{
  listEntities(organism = "ecoli", annotation = "Pathway")
}
\arguments{
  \item{organism}{A character containing the organism id:
  use \code{\link{listOrganisms}} to display the available
  organisms.}

  \item{annotation}{A character containing the selected
  annotation type: use \code{\link{listAnnotationTypes}} to
  display the available types.}
}
\value{
  A vector containing the available entities for the
  selected annotation type.
}
\description{
  This method takes a string containing the nickname for
  the selected organism and a vector of string,
  representing the genes of interest for the specified
  organism, and returns a list containing the locustags,
  gene_names, contrasts and M-values for the current
  selection.
}
\examples{
\dontrun{
 library("Rcolombos")
 bsubt_available_pathway_entities <- listEntities(organism="bsubt", annotation="Pathway")
}
}
\references{
  http://colombos.net
}

