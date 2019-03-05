#' Transforms a VSURF feature selection object into a matrix table
#'
#' @param vobj A VSURF package feature selection output
#' @param classifdat The classification dataframe without factor column
#' @return A matrix that can be exported as a (e.g.) .csv table
#' @examples
#' 
#' Run a VSURF feature selection:
#' 
#'   fs.I <- VSURF(classif.allclasses[,2:10], 
#'   classif.allclasses[,1], 
#'   clusterType = "FORK", 
#'   ntree = 500,mtry = 4) #All feature selection parameters are only examples.
#' 
#' Use presented function:
#' 
#'     varimp_I <- VSURF_table(fs.I, classif.df[,-1])
#'     
#' Export as .csv:
#' 
#'   write.csv(varimp_I, 'output/example.csv', row.names = TRUE)


VSURF_table<- function(vobj,classifdat){
  
  fs.norm.1 <- (vobj$imp.varselect.thres-min(vobj$imp.varselect.thres))/
    (max(vobj$imp.varselect.thres)-min(vobj$imp.varselect.thres)) #44
  
  out <- classifdat
  varimp.canopy <- rbind(names(out[,as.numeric(vobj$varselect.thres)]),
                         round(vobj$imp.varselect.thres,2), #45
                         round(fs.norm.1, 2)) #46
  
  row.names(varimp.canopy) <- c('Predictor', 'Abs. Imp.', 'Rel. Imp.') #47
  
  
   return(varimp.canopy)
  
}

