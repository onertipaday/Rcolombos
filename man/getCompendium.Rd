\name{getCompendium}
\alias{getCompendium}
\title{This method allows to download/import the full compendium for the selected organism}
\usage{
  getCompendium(organism = "ecoli", path = NULL)
}
\arguments{
  \item{organism}{A character containing the organism id:
  use \code{\link{listOrganisms}} to display the available
  organisms.}

  \item{path}{A string indicating the path where the file
  will be either downloaded or read, if already retrieved}
}
\value{
  A list containing three data.frame: \item{exprdata}{the
  full compendium for the selected organism}
  \item{condannot}{The condition annotation for the
  selected organism} \item{condontol}{the condition
  ontology for the selected organism}
}
\description{
  This method allows to download/import the full compendium
  for the selected organism
}
\examples{
\dontrun{
library('Rcolombos')
mtube <- getCompendium("mtube")
}
}
\references{
  http://colombos.net
}

