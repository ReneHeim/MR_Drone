#' Dropping factor levels
#'
#' @param caretmodel A prediction model created with the caret pkg
#' @return A data frame containing the confusion matrix incl UA and PA
#' @examples
#' caretout(carmodel.1)
#' 
caretout <- function(caretmodel) {
  con <- round(caretmodel$confusion$byClass[, c(1, 3)],3)*100
  colnames(con) <- c('Producer Accuracy [%]', 'User Accuracy [%]')
  
  matri <- caretmodel$confusion$table
  
  colnames(matri) <- paste("REF", colnames(matri), sep = ".")
  rownames(matri) <- paste("PRED", rownames(matri), sep = ".")
  maco <- cbind(matri, con)
  
  maco <- as.data.frame(maco)
  return(maco)
}



