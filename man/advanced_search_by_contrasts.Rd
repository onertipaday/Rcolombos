\name{advanced_search_by_contrasts}
\alias{advanced_search_by_contrasts}
\title{Accessory function allowing the advanced_search by contrast_ids, go, experiment, condition}
\usage{
  advanced_search_by_contrasts(organism = NULL, ids = NULL,
    geneNames = FALSE, c_search_type = NULL)
}
\arguments{
  \item{organism}{A character containing the organism id:
  use \code{\link{listOrganisms}} to display the available
  organisms.}

  \item{ids}{A vector of strings representing contrast_id,
  go terms, experiment id or condition id according the
  search type.}

  \item{geneNames}{boolean if == FALSE (default) return the
  locustag otherwise the gene_name for the selected genes.}

  \item{c_search_type}{A string either contrast_names.
  experiment, go, condition.}
}
\value{
  A data.frame containing locustag (gene_names), contrasts
  and M-values for the current organism and genes.
}
\description{
  Accessory function allowing the advanced_search by
  contrast_ids, go, experiment, condition
}
\examples{
\dontrun{
c.cn <- advanced_search_by_contrasts(organism="bsubt",
ids=c("GSM27217.ch2-vs-GSM27217.ch1","GSM27218.ch1-vs-GSM27218.ch2"), search_type="contrast_names")
}
}
\references{
  http://colombos.net
}

