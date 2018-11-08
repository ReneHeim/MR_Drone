#' Reformat wide to long data for ggplot2
#'
#' @param datawide A wide data frame
#' @param agg If FALSE (default) all spectra will be included in the long df. If
#' TRUE, spectra will be aggregated by their mean.
#' @return A long data frame
#' @examples
#' DropClass(ori.data, ori.data$Type, "Healthy")
#' Note: This function is custom made for spectral data and will rename the
#' columns of the return object ('Wavelength', 'Reflectance'). It is important 
#' the column names are numbers for each waveband only and we must prevent R 
#' checking the names when loading a spectral dataset. R would add an X in front 
#' of each number as this is what R thinks is good name for a column.

prep_gg <- function(datawide, agg=FALSE){
  require(reshape2)
  
  if(agg == TRUE) {
    datawide <- aggregate(.~Type, data=datawide, mean)
  }
  
  datamelt<-melt(datawide, id=c("Type"))
  names(datamelt)[2:length(names(datamelt))] <-  c('Wavelength', 'Reflectance')
  #datamelt$Wavelength<-gsub("X", "", paste(datamelt$Wavelength))
  datamelt$Wavelength <- as.character(datamelt$Wavelength)
  datamelt$Wavelength <- as.numeric(datamelt$Wavelength)
        
        datamelt
}

