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
#' @return A list containing three data.frame: 
#' \item{exprdata}{the full compendium for the selected organism}
#' \item{condannot}{The condition annotation for the selected organism}
#' \item{condontol}{the condition ontology for the selected organism}
#'
#' @references http://colombos.net
#'
#' @export
#'
#' @examples
#' \dontrun{
#' library('Rcolombos')
#' hpylo <- getCompendium("hpylo")
#' }
#'
getCompendium <- function(organism="hpylo", path=NULL){
    if(is.null(path)) path <- getwd() else {}
    destfile <- paste(path,"/",organism, "_compendium_data.zip",sep="")
    if(!file.exists(destfile)){
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
        } else {
            tmp <- fromJSON(output$data, nullValue = NA)$data;
            download.file( tmp, destfile )
            return(parseCompendium(destfile))
        }
    } else {
        return(parseCompendium(destfile))
    }
}
#' This method allows importing the full compendium for the selected organism from a local file
#'
#' @param destfile A character containing the full path of the downloaded file
#'
#' @return A list containing three data.frame: 
#' \item{exprdata}{the full compendium for the selected organism}
#' \item{condannot}{The condition annotation for the selected organism}
#' \item{condontol}{the condition ontology for the selected organism}
#'
#' @references http://colombos.net
#'
#' @export
#'
#' @examples
#' \dontrun{
#' library('Rcolombos')
#' mtube <- parseCompendium("mtube_compendium_data.zip")
#' }
#'
parseCompendium <- function(destfile){
    out_dir <- strsplit(destfile, "\\.")[[1]][1]
    unzip(destfile, exdir=out_dir) # unzip the files in the proper directory 
    files <- dir(path=out_dir,pattern="colombos_[a-z]+_[a-z]+_[0-9]+.txt") # reg exp for matching only the downloaded files
    temp <- paste(out_dir, files[grep("colombos_[a-z]+_exprdata_[0-9]+.txt", files)], sep="/")
    my_cols <- na.omit(scan(temp, nlines=1, sep="\t", what="c", na.strings="", quiet=TRUE))
    exprdata <- read.csv(temp, row.names=1, skip=7, stringsAsFactors=FALSE, sep="\t", header=FALSE)
    exprdata <- exprdata[,c(2:dim(exprdata)[[2]])] 
    colnames(exprdata) = my_cols; exprdata <- exprdata[,c(2:dim(exprdata)[[2]])]
    ## condition annotations 
    temp <- paste(out_dir, files[grep("colombos_[a-z]+_condannot_[0-9]+.txt", files)], sep="/")
    condannot <- read.csv(temp, stringsAsFactors=FALSE, sep="\t", header=T)
    ## condition ontology
    temp <- paste(out_dir, files[grep("colombos_[a-z]+_condontol_[0-9]+.txt", files)], sep="/")
    condontol <- read.csv(temp, stringsAsFactors=FALSE, sep="\t", header=T)
    ## return a list with three data.frame
    return( list(exprdata=exprdata, condannot=condannot, condontol=condontol) )
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
