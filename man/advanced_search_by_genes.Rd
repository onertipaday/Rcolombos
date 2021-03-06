% Generated by roxygen2 (4.1.0): do not edit by hand
% Please edit documentation in R/Rcolombos.R
\name{advanced_search_by_genes}
\alias{advanced_search_by_genes}
\title{Accessory function allowing the advanced_search by gene_ids, go, annotation}
\usage{
advanced_search_by_genes(organism = "bsubt", ids = NULL,
  geneNames = FALSE, g_search_type = "genes", ann_type)
}
\arguments{
\item{organism}{A character containing the organism id: use \code{\link{listOrganisms}} to display
the available organisms.}

\item{ids}{A vector of strings representing gene_id, go terms or annotation entities according the search type.}

\item{geneNames}{boolean if == FALSE (default) return the locustag otherwise the gene_name for the selected genes.}

\item{g_search_type}{A string either genes, go or annotation.}

\item{ann_type}{A string containing the selected gene_annotation_type: use \code{\link{listEntities}} to display the available entities.}
}
\value{
A data.frame containing locustag (gene_names),
contrasts and M-values for the current organism and genes.
}
\description{
Accessory function allowing the advanced_search by gene_ids, go, annotation
}
\references{
http://colombos.net
}

