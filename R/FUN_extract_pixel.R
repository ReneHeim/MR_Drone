####
# Function accepts a multi-layered RasterBrick to extract pixel values from 
# pre-defined polygons saved as a GIS shape file. Then a df is created, 
# containing an ID and classification labels for each pixel.
# The user must define a "Treatment" column and an "ID" column when drawing the 
# polygons by using a GIS software (e.g. QGIS).
####

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