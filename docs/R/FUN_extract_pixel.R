#' Extract pixel values from a .tif file
#'
#' @param rsta Multiband raster file (e.g. raster stack or brick)
#' @param shp SpatialPolygonsDataFrame containing sampling polygons (ID and 
#' Treatment columns required)
#' @return SpatialPolygonsDataFrame conatining all values that were captured by
#' the polygons and also contains an ID and Treatment column
#' @examples
#' img <- brick("data/anymultilayerimage.tif")
#' poly <- shapefile("data/samplepolygons.shp") #were sampled on img using GIS
#' df <- extract_pixel(img, poly)
extract_pixel <- function(rsta, shp){
  
  responseCol <- "id"
  
  dfAll<-data.frame(matrix(vector(), nrow = 0, ncol = length(names(rsta)) + 1))
  #container
  
  for (i in 1:length(unique(shp[[responseCol]]))){                          
    category <- unique(shp[[responseCol]])[i]
    categorymap <- shp[shp[[responseCol]] == category,]
    dataSet <- raster::extract(rsta, categorymap)
    dataSet <- dataSet[!unlist(lapply(dataSet, is.null))]
    dataSet <- lapply(dataSet, function(x){cbind(x, class = as.numeric(rep(category, nrow(x))))})
    df <- do.call("rbind", dataSet)
    dfAll <- rbind(dfAll, df)
  }

  return(dfAll)
}