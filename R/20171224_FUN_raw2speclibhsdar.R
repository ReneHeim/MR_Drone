raw2speclib <- function(dat){
    #fct converts spectral data into a hsdar spectral lib if first col contains
    #ID variable (here:Type) and all other colums are numeric, named with the 
    #number of the wavelength only (e.g. 350, 355, 365,...) and contain refl
    #values
    require(hsdar)
    
    data.new <- t(dat[,-1])
    colnames(data.new) <- dat[,1]
    rownames(data.new) <-  c()
    
    Wavelength <- colnames(dat[,2:length(dat)])
    Wavelength <- as.numeric(Wavelength)
    
    # NB -- have investigated warning here and it is fine
    spectra <- suppressWarnings(speclib(data.new,Wavelength))
    
    
    mat <- as.matrix(dat[,1])
    colnames(mat) <- colnames(dat[1])
    
    #attribute(spectra) <- mat 
    
    return(spectra)
    
    
}