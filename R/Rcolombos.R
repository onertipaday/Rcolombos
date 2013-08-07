#' This method mimics the quick_search functionality of Colombos.
#' It takes a string containg the nickname for the selected organism and a vector of string, 
#' representing the genes of interest for the specified organism, 
#' and returns a list containing the locustags, gene_names, 
#' contrasts and M-values for the current selection.
#'
#' @param  organism A character containing the organism id: use \code{\link{listOrganisms}} to display
#' the available organisms.
#' @param genes A vector of strings representing the genes of interest.
#' @param geneNames boolean if == FALSE (default) return the locustag otherwise the gene_name for the selected genes.
#' 
#' 
#' @return A data.frame containing locustag (gene_names), 
#' contrasts and M-values for the current organism and genes.
#'
#' @references http://colombos.net
#'
#' @export
#'
#' @examples
#' \dontrun{
#'  library("Rcolombos")
#'  my_module <- quick_search(organism="ecoli", genes=c("b0400","b2805","b0567"), geneNames=FALSE)
#'  heatmap(as.matrix(my_module), col=terrain.colors(15))
#' }
#'
quick_search <- function(organism="ecoli", genes, geneNames=FALSE){
  if(is.null(genes)) stop("Insert a character vector with the genes to be imputed.") else {}
  #
  t <- basicTextGatherer()
  h <- basicHeaderGatherer()
  curlPerform( url = paste("http://rest.colombos.net/rest/quick_search/", 
                           organism, paste(genes, collapse=","), sep="/"),
               .opts = list(httpheader = c('Content-Type' = "application/json"), verbose = FALSE),
               curl = getCurlHandle(),
               writefunction = t$update,
               headerfunction = h$update
  ) 
  
  output <- list(data = t$value(),
                 status = h$value()[['status']],
                 status.message = h$value()[['statusMessage']])
  httpstatus <- as.numeric(output$status)
  if (httpstatus != 200) {
    return(output$status.message)
  }  else {
      tmp <- fromJSON(output$data, nullValue = NA)$data;
      if (geneNames){
        response = suppressWarnings( data.frame(matrix(as.numeric(tmp$Mvalues), 
                   nrow=length(tmp$geneNames), ncol=length(tmp$contrasts), byrow=T)))
        colnames(response) <- tmp$contrasts; rownames(response) <- tmp$geneNames
      } else {
        response = suppressWarnings( data.frame(matrix(as.numeric(tmp$Mvalues), 
                   nrow=length(tmp$geneLocustags), ncol=length(tmp$contrasts), byrow=T)) )
        colnames(response) <- tmp$contrasts; rownames(response) <- tmp$geneLocustags
      }
      return(response)
  }
}
