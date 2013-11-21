\name{parseCompendium}
\alias{parseCompendium}
\title{This method allows importing the full compendium for the selected organism from a local file}
\usage{
  parseCompendium(destfile)
}
\arguments{
  \item{destfile}{A character containing the full path of
  the downloaded file}
}
\value{
  A list containing three data.frame: \item{exprdata}{the
  full compendium for the selected organism}
  \item{condannot}{The condition annotation for the
  selected organism} \item{condontol}{the condition
  ontology for the selected organism}
}
\description{
  This method allows importing the full compendium for the
  selected organism from a local file
}
\examples{
\dontrun{
library('Rcolombos')
mtube <- parseCompendium("mtube_compendium_data.zip")
}
}
\references{
  http://colombos.net
}

