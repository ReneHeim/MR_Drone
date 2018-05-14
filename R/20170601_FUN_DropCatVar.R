####
# This function is dropping a categorical response variable in a df and is also 
# dropping the according factor level.
####

DropClass <-  function(data, column, classtodrop){
        
        require('gdata')
        
        data2 <- subset(data, column != classtodrop)
        
        #Check if levels are structured correctly and adapt if necessary (gdata) #
        
        data2 <- drop.levels(data2)
        
        return(data2)
}

