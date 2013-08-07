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

#' This method mimics the advanced_search functionality of Colombos.
#' It takes a string containg the nickname for the selected organism and a vector of string, 
#' representing the genes of interest for the specified organism, 
#' and returns a list containing the locustags, gene_names, 
#' contrasts and M-values for the current selection.
#'
#' @param  organism A character containing the organism id: use \code{\link{listOrganisms}} to display
#' the available organisms.
#' @param genes A vector of strings representing the genes of interest.
#' @param geneNames boolean if == FALSE (default) return the locustag otherwise the gene_name for the selected genes.
#' @param
#' @param
#' @param by A string, either "genes", "contrasts" or "both" (default is genes)
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
#'  g.gn <- advanced_search(organism="bsubt", gene_ids=c("cgeB","yfnG"), by="genes", search_type="genes")
#'  g.go <- advanced_search(organism="bsubt", gene_ids="response to antibiotic, transcription", by="genes", search_type="go")
#'  g.anno <- advanced_search(organism="bsubt", gene_ids="biotin-carboxyl carrier protein assembly", 
#'  by="genes", search_type="annotation", ann_type="Pathway")
#'  
#'  b.gn.cn <- advanced_search(organism="bsubt", gene_ids=c("cgeB","yfnG"), 
#'  geneNames=FALSE, contrast_ids=c("GSM27217.ch2-vs-GSM27217.ch1","GSM27218.ch1-vs-GSM27218.ch2", by="both")
#'  heatmap(as.matrix(my_module), col=terrain.colors(15))
#' }
#'
advanced_search <- function(organism="bsubt", gene_ids=c("cgeB","yfnG"), geneNames=FALSE, contrast_ids, by="genes", search_type="genes", ann_type){
    switch(by,
           genes = advanced_search_by_genes(organism, gene_ids, geneNames, search_type, ann_type),
           contrasts = "contrasts",
           both = "both"
    )
}
#' Accessory function (not exposed) allowing the advanced_search by gene_ids, go, annotation
#'
#' @param  organism A character containing the organism id: use \code{\link{listOrganisms}} to display
#' the available organisms.
#' @param ids A vector of strings representing gene_id, go terms or annotation entities according the search type.
#' @param geneNames boolean if == FALSE (default) return the locustag otherwise the gene_name for the selected genes.
#' @param search_type A string either genes, go or annotation
#' @param entity 
#' 
#' @return A data.frame containing locustag (gene_names), 
#' contrasts and M-values for the current organism and genes.
#'
#' @references http://colombos.net
#'
#' @export
#'
advanced_search_by_genes <- function(organism="bsubt", ids=c("cgeB","yfnG"), geneNames=FALSE, search_type="genes", ann_type){
    #         if(is.null(genes)) stop("Insert a character vector with the genes to be imputed.") else {}
    t <- basicTextGatherer()
    h <- basicHeaderGatherer()
    if(search_type=="genes"){
        url_string <- paste("http://rest.colombos.net/advanced_search_by_genes", organism, "genes", paste(ids, collapse=","), sep="/")
    }
    else if(search_type=="go"){
        url_string <- paste("http://rest.colombos.net/advanced_search_by_genes", organism, "go", gsub(" ","%20", ids), sep="/")
    }
    else if(search_type=="annotation"){
        url_string <- paste("http://rest.colombos.net/advanced_search_by_genes", organism, "annotation", gsub(" ","%20", ann_type), 
                            gsub(" ","%20", ids), sep="/")
        
    } else {
        stop("Wrong search_type: it should be either genes, go or annotation!")
    }
    curlPerform( url = url_string,
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
#' Accessory function (not exposed) allowing the advanced_search by contrast_ids, go, experiment, condition
#'
#' @param  organism A character containing the organism id: use \code{\link{listOrganisms}} to display
#' the available organisms.
#' @param ids A vector of strings representing gene_id, go terms or annotation entities according the search type.
#' @param geneNames boolean if == FALSE (default) return the locustag otherwise the gene_name for the selected genes.
#' @param search_type A string either genes, go or annotation
#' @param entity 
#' 
#' @return A data.frame containing locustag (gene_names), 
#' contrasts and M-values for the current organism and genes.
#'
#' @references http://colombos.net
#'
#' @export
#'
advanced_search_by_contrasts <- function(organism="bsubt", ids=c("cgeB","yfnG"), geneNames=FALSE, search_type="genes", ann_type){
    #         if(is.null(genes)) stop("Insert a character vector with the genes to be imputed.") else {}
    t <- basicTextGatherer()
    h <- basicHeaderGatherer()
    if(search_type=="genes"){
        url_string <- paste("http://rest.colombos.net/advanced_search_by_genes", organism, "genes", paste(ids, collapse=","), sep="/")
    }
    else if(search_type=="go"){
        url_string <- paste("http://rest.colombos.net/advanced_search_by_genes", organism, "go", gsub(" ","%20", ids), sep="/")
    }
    else if(search_type=="annotation"){
        url_string <- paste("http://rest.colombos.net/advanced_search_by_genes", organism, "annotation", gsub(" ","%20", ann_type), 
                            gsub(" ","%20", ids), sep="/")
        
    } else {
        stop("Wrong search_type: it should be either genes, go or annotation!")
    }
    curlPerform( url = url_string,
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

advanced_search_by_both <- function(organism="bsubt", gene_ids=c("cgeB","yfnG"), geneNames=FALSE, contrast_ids, by="genes"){
    
}
