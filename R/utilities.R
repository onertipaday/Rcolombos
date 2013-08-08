#' Returns a character vector corresponding to the currently available organisms.
#'
#' @return A list containing the currently available organisms.
#'
#' @references http://colombos.net
#'
#' @export
#'
#' @examples
#' \dontrun{
#' library('Rcolombos')
#' listOrganisms()
#' }
#'
listOrganisms <- function(){
  header.field = c('Content-Type' = "application/json")
  curl <- getCurlHandle() 
  curlSetOpt(.opts = list(httpheader = header.field, verbose = FALSE), curl = curl) 
  t <- basicTextGatherer()
  h <- basicHeaderGatherer()
  body = curlPerform(url = "http://rest.colombos.net/get_organisms",
                     curl = curl,
                     writefunction = t$update,
                     headerfunction = h$update)
  output <- list(data = t$value(),
                 status = h$value()[['status']],
                 status.message = h$value()[['statusMessage']])
  httpstatus <- as.numeric(output$status)
  if (httpstatus != 200) {
    return(output$status.message)
  } else {
    tmp <- fromJSON(output$data, nullValue = NA)$data;
    response <- data.frame(matrix(unlist(tmp, recursive=F), length(tmp), 3, byrow=T))
    for (i in 1:3) {
      response[,i] <- sapply(response[,i], as.character)
    }
    colnames(response) <- c("name", "description", "nickname")
    return(response)
  }
}

#' This method takes as parameter a single string, representing an organism,
#' and returns a character vector corresponding to the currently available organisms.
#'
#' @param organism A character containing the organism id: use \code{\link{listOrganisms}} to display
#' the available organisms.
#'
#' @return A data.frame containing the locustag and description
#' of all the genes for the selected organism.
#'
#' @references http://colombos.net
#'
#' @export
#'
#' @examples
#' \dontrun{
#' library('Rcolombos')
#' listGenes()
#' }
#'
listGenes <- function(organism="ecoli"){
  header.field = c('Content-Type' = "application/json")
  curl <- getCurlHandle() 
  curlSetOpt(.opts = list(httpheader = header.field, verbose = FALSE), curl = curl) 
  t <- basicTextGatherer()
  h <- basicHeaderGatherer()
  body = curlPerform(url=paste("http://rest.colombos.net/get_genes/", organism, sep=""),
                     curl = curl,
                     writefunction = t$update,
                     headerfunction = h$update)
  output <- list(data = t$value(),
                 status = h$value()[['status']],
                 status.message = h$value()[['statusMessage']])
  httpstatus <- as.numeric(output$status)
  if (httpstatus != 200) {
    return(output$status.message)
  } else {
    tmp <- fromJSON(output$data, nullValue = NA)$data;
    response <- data.frame(matrix(unlist(tmp, recursive=F), length(tmp), 2, byrow=T))
    for (i in 1:2) {
      response[,i] <- sapply(response[,i], as.character) # sanitize the list in the df
    }
    colnames(response) <- c("locustag", "gene_name")
    return(response)
  }
}

#' This method takes as parameter a single string, representing an organism,
#' and returns a character vector corresponding to the currently available organisms.
#'
#' @param organism A character containing the organism id: use \code{\link{listOrganisms}} to display
#' the available organisms.
#'
#' @return A data.frame containing the contrasts and GSM
#' of all the contrasts for the selected organism.
#'
#' @references http://colombos.net
#'
#' @export
#'
#' @examples
#' \dontrun{
#' library('Rcolombos')
#' listContrasts()
#' }
#'
listContrasts <- function(organism="ecoli"){
  header.field = c('Content-Type' = "application/json")
  curl <- getCurlHandle() 
  curlSetOpt(.opts = list(httpheader = header.field, verbose = FALSE), curl = curl) 
  t <- basicTextGatherer()
  h <- basicHeaderGatherer()
  body = curlPerform(url=paste("http://rest.colombos.net/get_contrasts/", organism, sep=""),
                     curl = curl,
                     writefunction = t$update,
                     headerfunction = h$update)
  output <- list(data = t$value(),
                 status = h$value()[['status']],
                 status.message = h$value()[['statusMessage']])
  httpstatus <- as.numeric(output$status)
  if (httpstatus != 200) {
    return(output$status.message)
  } else {
    tmp <- fromJSON(output$data, nullValue = NA)$data;
    response <- data.frame(matrix(unlist(tmp, recursive=F), length(tmp), 2, byrow=T))
    for (i in 1:2) {
      response[,i] <- sapply(response[,i], as.character) # sanitize the list in the df
    }
    colnames(response) <- c("name", "description")
    return(response)
  }
}

#' This method allows to download/import the full compendium for the selected organism
#'
#' @param organism A character containing the organism id: use \code{\link{listOrganisms}} to display
#' the available organisms.
#'
#' @param path A string indicating the path where the file will be either downloaded or read,
#' if already retrieved
#' 
#'
#' @return A data.frame containing the full compendium for the selected organism.
#'
#' @references http://colombos.net
#'
#' @export
#'
#' @examples
#' \dontrun{
#' library('Rcolombos')
#' getCompendium()
#' }
#'
getCompendium <- function(organism="ecoli", path=NULL){
  if(is.null(path)) path <- getwd() else {}
    destfile <- paste(path,"/",organism, "_compendium_data.txt.gz",sep="")
    header.field = c('Content-Type' = "application/json")
    curl <- getCurlHandle() 
    curlSetOpt(.opts = list(httpheader = header.field, verbose = FALSE), curl = curl) 
    t <- basicTextGatherer()
    h <- basicHeaderGatherer()
    body = curlPerform(url=paste("http://rest.colombos.net/get_organism_data/", organism, "/",sep=""),
                     curl = curl,
                     writefunction = t$update,
                     headerfunction = h$update)
    output <- list(data = t$value(),
                 status = h$value()[['status']],
                 status.message = h$value()[['statusMessage']])
    httpstatus <- as.numeric(output$status)
    if (httpstatus != 200) {
      return(output$status.message)
    } else {
      if(!file.exists(destfile)){
        tmp <- fromJSON(output$data, nullValue = NA)$data;
        download.file( tmp, destfile )
#       return(tmp)
      } else {}
    ## TO.DO create a ExpressionSet object take care of all contrasts and genes info
    my_cols <- scan(destfile, what="c", nlines=1)
    out <- read.csv(destfile, row.names=1,
                    skip=7, stringsAsFactors=FALSE, sep="\t")
    out <- out[,c(2:dim(out)[[2]])] 
    colnames(out) = my_cols; out <- out[,c(2:dim(out)[[2]])]
  }    
}

#' This method takes as parameter a single string, representing an organism,
#' and returns a character vector corresponding to the currently available organisms.
#'
#' @param organism A character containing the organism id: use \code{\link{listOrganisms}} to display
#' the available organisms.
#'
#' @return A data.frame containing the locustag and description
#' of all the genes for the selected organism.
#'
#' @references http://colombos.net
#'
#' @export
#'
#' @examples
#' \dontrun{
#' library('Rcolombos')
#' listGenes()
#' }
#'
listGenes <- function(organism="ecoli"){
    header.field = c('Content-Type' = "application/json")
    curl <- getCurlHandle() 
    curlSetOpt(.opts = list(httpheader = header.field, verbose = FALSE), curl = curl) 
    t <- basicTextGatherer()
    h <- basicHeaderGatherer()
    body = curlPerform(url=paste("http://rest.colombos.net/get_genes/", organism, sep=""),
                       curl = curl,
                       writefunction = t$update,
                       headerfunction = h$update)
    output <- list(data = t$value(),
                   status = h$value()[['status']],
                   status.message = h$value()[['statusMessage']])
    httpstatus <- as.numeric(output$status)
    if (httpstatus != 200) {
        return(output$status.message)
    } else {
        tmp <- fromJSON(output$data, nullValue = NA)$data;
        response <- data.frame(matrix(unlist(tmp, recursive=F), length(tmp), 2, byrow=T))
        for (i in 1:2) {
            response[,i] <- sapply(response[,i], as.character) # sanitize the list in the df
        }
        colnames(response) <- c("locustag", "gene_name")
        return(response)
    }
}

#' This method takes as parameter a string (the nickname of an organism) and returns a character vector 
#' corresponding to the currently available annotation type for the selected organism.
#'
#' @param organism A character containing the organism id: use \code{\link{listOrganisms}} to display
#' the available organisms.
#'
#' @return A data.frame containing the name and description of the annotation
#' for the selected organism.
#'
#' @references http://colombos.net
#'
#' @export
#'
#' @examples
#' \dontrun{
#' library('Rcolombos')
#' listAnnotationTypes()
#' }
#'
listAnnotationTypes <- function(organism="ecoli"){
    header.field = c('Content-Type' = "application/json")
    curl <- getCurlHandle() 
    curlSetOpt(.opts = list(httpheader = header.field, verbose = FALSE), curl = curl) 
    t <- basicTextGatherer()
    h <- basicHeaderGatherer()
    body = curlPerform(url=paste("http://rest.colombos.net/get_annotation_types/", organism, sep=""),
                       curl = curl,
                       writefunction = t$update,
                       headerfunction = h$update)
    output <- list(data = t$value(),
                   status = h$value()[['status']],
                   status.message = h$value()[['statusMessage']])
    httpstatus <- as.numeric(output$status)
    if (httpstatus != 200) {
        return(output$status.message)
    } else {
        tmp <- fromJSON(output$data, nullValue = NA)$data;
        response <- data.frame(matrix(unlist(tmp, recursive=F), length(tmp), 2, byrow=T))
        for (i in 1:2) {
            response[,i] <- sapply(response[,i], as.character) # sanitize the list in the df
        }
        colnames(response) <- c("name", "description")
        return(response)
    }
}

#' This method takes a string containing the nickname for the selected organism and a string containing the annotation type
#' and return the available entities
#'
#' @param  organism A character containing the organism id: use \code{\link{listOrganisms}} to display
#' the available organisms.
#' @param annotation A character containing the selected annotation type: use \code{\link{listAnnotationTypes}} to display
#' the available types.
#' 
#' @return A vector containing the available entities for the selected annotation type.
#'
#' @references http://colombos.net
#'
#' @export
#'
#' @examples
#' \dontrun{
#'  library("Rcolombos")
#'  pathway_entities <- listEntities(organism="bsubt", annotation="Pathway")
#'  Tr_entities <- listEntities("bsubt","Transcriptional regulation")
#' }
#'
listEntities <- function(organism="ecoli", annotation="Pathway"){
    if(is.null(annotation)) stop("Insert a string with the annotation type.\n See listAnnotationTypes for the available types.") else {}
    #
    t <- basicTextGatherer()
    h <- basicHeaderGatherer()
    curlPerform( url = paste("http://rest.colombos.net/get_entities", organism, gsub(" ","%20", annotation), sep="/"),
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
            response = suppressWarnings( unlist(tmp) )

        return(response)
    }
}
