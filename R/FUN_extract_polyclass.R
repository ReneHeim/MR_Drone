#' Extract digital values from a raster file (e.g. geotiff).
#'
#' @param shp A GIS shapefile shapefile format can spatially describe vector 
#'            features: points, lines, and polygons, representing specific 
#'            landscape regions. Each feature usually has descriptive attributes. 
#'            In this case: id (1,2,3,...) and Class (Water, Soil, Forest,...).
#' @param rsta Image file which can be imported bz using the raster pkg
#' @return A data frame containing a class column (Class: Factor) and  one or 
#'         multiple predictor columns (Band1: num, Band2: num, ...).
#' @examples
#' extract_polyclass(shapepoly, image)
#' @details !This function only processes shapefiles where the attribute table 
#'          contains the attributes written as described above (i.e. id, Class)! 


extract_polyclass <- function(shp, rsta){
  
  require("raster")
  
  vec <- unique(shp@data$Class)
  
  res <- list()
  
  for (i in vec){
    
    Class <- subset(shp, shp$Class == i)
    
    dat <- raster::extract(rsta, Class)
    
    findf <- as.data.frame(do.call("rbind", dat))
    
    Class <- rep(i, each = length(findf[,1]))
    
    res[[i]] <- cbind(Class, findf)
    
  }
  
  fin <- do.call(rbind.data.frame, res)
  
  return(fin)
}
