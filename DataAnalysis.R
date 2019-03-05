####
# The following code reproduces results of the article:
# 
# Heim, Wright, Carnegie, Taylor, Scarth and Oldeland, 2019. Using drones and 
# multispectral imagery to detect myrtle rust on a lemon myrtle plantation. Drones
#
# The code is split in two major sections according to our research questions.
#
# A)  Is it possible to accurately discriminate fungicide-treated and untreated 
#     sunlit plants at canopy-level
#
# B)  What are important spectral regions for the classification?
####


# Before we start to analyse the data we set up our working environment.

dir.create('output', FALSE, FALSE)

# install.packages(c("rgdal", 
#                    "raster",
#                    "reshape2",
#                    "ggplot2", 
#                    "caret",
#                    "e1071", 
#                    "gdata",
#                    "rasterVis",
#                    "corrplot"))

library(rgdal)#extract reflectance values
library(corrplot)#test if predictor vars are correlated
library(raster)#extract reflectance values
library(ggplot2)#figure 3
library(caret)#machine learning library
library(gdata)#required for drop_class()
library(reshape2)#for prep_gg function, reshaping spectra to use ggplot2
library(rasterVis)#modify risk map (Part D)
library(VSURF)#feature selection



source("R/FUN_drop_cat_var.R")#drops factor and according factor level
source("R/FUN_extract_polyclass.R")#extract reflectance and index values
source("R/FUN_prepggwide2long.R")#reshape spectral data for ggplot
source("R/FUN_VSURF_table.R")#create feature selection output

# PART A------------------------Cleaning Data-----------------------------------

# Loading image brick.

  img <-
    brick("data/20180528_ortho_ground.tif")

# Convert from DN to reflectance.

  img <- img/65535

# Calculate indices and rename bands for clarity (MicaSense RedEdge Camera)

  names(img) <- c("blue", "green", "red", "re", 
                "nir", "alpha")
  
    ndvi <- (img$nir-img$red)/(img$nir+img$red)
    rg <- img$red/img$green
    sipi <- (img$nir-img$blue)/(img$nir-img$red)
    ari <- (1/img$green)-(1/img$re)

      img <- addLayer(img, c(ndvi, rg, sipi, ari))

  names(img) <- c("blue", "green", "red", "re", 
                "nir", "alpha", "ndvi", "rg", "sipi", "ari")

# Loading QGIS shape file where sample polygons have been defined.

  allclasses <- shapefile("data/samplepolygons.shp")
  
# Check what classes are contained

  unique(allclasses@data$Class)
    
# Extract reflectance values from pixels based on sample polygons.

  allclasses.df<- extract_polyclass(allclasses, img)
  
# Check if all classes were imported

  unique(allclasses.df$Class) 

# Check for NAs

  apply(allclasses.df, 2, function(x) any(is.na(x)))

  #full.df.wona <- na.omit(full.df) If NAs contained, they can be removed.

  summary(allclasses.df$Class)
  #summary(full.df.wona$Class)
  #summary(full.df$Class)-summary(full.df.wona$Class)

  allclasses.noalpha.df <- 
    allclasses.df[ , -which(names(allclasses.df) %in% c("alpha"))] #remove alpha

# Write full data to .csv and reload with new name for classification.

  write.csv(allclasses.noalpha.df, 
            'output/droneclassif_plusSHD.csv', 
            row.names = FALSE)
  
  classif.allclasses <- read.csv("output/droneclassif_plusSHD.csv", 
                                 check.names = FALSE)

# Preparing data for a second, mixed shadow, classification

SHDmixdata <- shapefile("data/samplepolygons_noSHD.shp")

  unique(SHDmixdata@data$Class) # To check what classes are contained

  SHDmixdata.df<- extract_polyclass(SHDmixdata, img)

unique(SHDmixdata.df$Class)

# Check for NAs

  apply(SHDmixdata.df, 2, function(x) any(is.na(x)))

  #full.df.wona <- na.omit(full.df)

  summary(SHDmixdata.df$Class)
  #summary(full.df.wona$Class)
  #summary(full.df$Class)-summary(full.df.wona$Class)

  SHDmixdata.noalpha.df <- 
    SHDmixdata.df[ , -which(names(SHDmixdata.df) %in% c("alpha"))]

# Write mix shadow data to .csv and reload with new name for classification.
  
  write.csv(SHDmixdata.noalpha.df, 
          'output/droneclassif_mixSHD.csv', 
          row.names = FALSE)

  classif.mixSHD <- read.csv("output/droneclassif_mixSHD.csv", 
                     check.names = FALSE)


# Plot Figure 3 (Multispectral signatures - TR, UN, SHD)

  fig3 <- classif.allclasses[,1:6]

    names(fig3) <- c("Type","475", "560", "668", "717", "840")

      spectragg <- prep_gg(fig3, agg = TRUE)

C <- ggplot(spectragg, aes(Wavelength, Reflectance*100, colour = Type)) +
  geom_line(aes(linetype=Type), size = 1)+
  geom_point(aes(shape=Type), size = 2)+
  scale_color_manual(values=c("#525252", "#00CD00", "#cd5c00"))+
  labs(x = "", y = "")+
  theme_minimal(base_size=20)+
  theme(legend.title = element_blank(),
        legend.key.width = unit(1.5, "cm"))+
  annotate("rect",
    xmin = 455,
    xmax = 495,
    ymin = -Inf,
    ymax = Inf,
    alpha = .2,
    fill = 'blue'
  ) +
  annotate(
    "rect",
    xmin = 540,
    xmax = 580,
    ymin = -Inf,
    ymax = Inf,
    alpha = .2,
    fill = 'green'
  ) +
  annotate(
    "rect",
    xmin = 658,
    xmax = 678,
    ymin = -Inf,
    ymax = Inf,
    alpha = .2,
    fill = c("#FF2626")
  ) +
  annotate(
    "rect",
    xmin = 707,
    xmax = 727,
    ymin = -Inf,
    ymax = Inf,
    alpha = .2,
    fill = c("#7A0A0A")
  ) +
  annotate(
    "rect",
    xmin = 800,
    xmax = 880,
    ymin = -Inf,
    ymax = Inf,
    alpha = .2,
    fill = c("lightgrey")
  )  +
  annotate(
    "text",
    x = 475,
    y = 80,
    label = "B",
    fontface = "bold",
    size = 5
  ) +
  annotate(
    "text",
    x = 560,
    y = 80,
    label = "G",
    fontface = "bold",
    size = 5
  ) +
  annotate(
    "text",
    x = 668,
    y = 80,
    label = "R",
    fontface = "bold",
    size = 5
  ) +
  annotate(
    "text",
    x = 717,
    y = 80,
    label = "RE",
    fontface = "bold",
    size = 5
  ) +
  annotate(
    "text",
    x = 840,
    y = 80,
    label = "NIR",
    fontface = "bold",
    size = 5
  )+
  labs(x = "Wavelength [nm]", y = "Reflectance [%]")
  

ggsave("output/Figure3.allspectra.png",
       plot = C,
       width = 40,
       height = 20,
       units = "cm",
       dpi = 400
)


# PART B------------------Random Forest Classification--------------------------

# All Classes (TR, UN, SHD) ----------------------------------------------------

set.seed(2019)

# Create test and training data

  inTraining <- createDataPartition(classif.allclasses$Class, 
                                    p = .75, 
                                    list = FALSE)
  
    train <- classif.allclasses[ inTraining,]
    test  <- classif.allclasses[-inTraining,]
    
# Define random forest model training parameters

  rfControl <- trainControl(
      method = "boot",
      number = 100)

  rfGrid <- expand.grid(mtry = seq(1, ncol(train)-1, 1))
    
# Train random forest model

  rfFit <- train(Class ~ ., data = train,
               method = "rf",
               importance = TRUE, ntree=500,
               trControl = rfControl, tuneGrid = rfGrid,
               metric = "Accuracy", maximize = TRUE)

# Validate random forest model on test data

  rfPred <- 
    predict.train(rfFit, test[, !names(test) %in% c("Class")], type = "raw")

# Store relevant random forest results in a list() container for proper output

  RFdataall <- 
    list(fit = rfFit,
      pred = predict.train(rfFit, 
                           test[, !names(test) %in% c("Class")], 
                           type = "raw"),
      confusion = confusionMatrix(rfPred, test$Class),
      varImp = varImp(rfFit, scale = TRUE))

# Use sink() to export list() container contents as text file
  
  sink(file = 'output/I_allclasses_RFclassif_report.txt')
    RFdataall
  sink()

# Save list() container as R data object to not have to rerun analysis

  saveRDS(RFdataall, 'output/I_allclasses_RFclassif_object.rds')
  #RFdataall <- readRDS("output/I_allclasses_RFclassif_object.rds")


# Mix Shadow (TR_S, UN_S) ------------------------------------------------------

# Create test and training data

  inTraining <- createDataPartition(classif.mixSHD$Class, p = .75, list = FALSE)

  train <- classif.mixSHD[ inTraining,]
  test  <- classif.mixSHD[-inTraining,]

# Train random forest model

  rfFit <- train(Class ~ ., data = train,
               method = "rf",
               importance = TRUE, ntree=500,
               trControl = rfControl, tuneGrid = rfGrid,
               metric = "Accuracy", maximize = TRUE)

# Validate random forest model on test data

  rfPred <- 
  predict.train(rfFit, test[, !names(test) %in% c("Class")], type = "raw")

# Store relevant random forest results in a list() container for proper output

  RFdatamix <- 
    list(fit = rfFit,
       pred = predict.train(rfFit, 
                            test[, !names(test) %in% c("Class")],
                            type = "raw"),
       confusion = confusionMatrix(rfPred, test$Class),
       varImp = varImp(rfFit, scale = TRUE))
  
# Use sink() to export list() container contents as text file

  sink(file = 'output/II_mix_RFclassif_report.txt')
    RFdatamix
  sink()

# Save list() container as R data object to not have to rerun analysis

  saveRDS(RFdata, 'output/II_mix_RFclassif_object.rds')
  #RFdatamix <- readRDS("output/II_mix_RFclassif_object.rds")


# Only UN and TR ---------------------------------------------------------------

# Drop shadow class

  classif.TRUN <- drop_class(classif.allclasses, 
                             classif.allclasses$Class, 
                             "SHD")
  
# Create test and training data

  inTraining <- createDataPartition(classif.TRUN$Class, p = .75, list = FALSE)

    train <- classif.TRUN[ inTraining,]
    test  <- classif.TRUN[-inTraining,]

# Train random forest model
  
  rfFit <- train(Class ~ ., data = train,
               method = "rf",
               importance = TRUE, ntree=500,
               trControl = rfControl, tuneGrid = rfGrid,
               metric = "Accuracy", maximize = TRUE)
  
# Validate random forest model on test data

  rfPred <- 
    predict.train(rfFit, test[, !names(test) %in% c("Class")], type = "raw")

# Store relevant random forest results in a list() container for proper output

  RFdata_TRUN <- 
    list(fit = rfFit,
       pred = predict.train(rfFit, 
                            test[, !names(test) %in% c("Class")], 
                            type = "raw"),
       confusion = confusionMatrix(rfPred, test$Class),
       varImp = varImp(rfFit, scale = TRUE)) #38

# Use sink() to export list() container contents as text file

  sink(file = 'output/III_onlyTRUN_RFclassif_report.txt')
    RFdata_TRUN
  sink()

# Save list() container as R data object to not have to rerun analysis

  saveRDS(RFdata_TRUN, 'output/III_onlyTRUN_RFclassif_object.rds') #40
  #RFdata_TRUN <- readRDS("output/III_onlyTRUN_RFclassif_object.rds") #41

# PART C----------------------Feature Selection---------------------------------

# Test for all data sets whether predictor variables are correlated.

corr <- round(cor(classif.allclasses[,2:10]), 2)
corr2 <- round(cor(classif.mixSHD[,2:10]), 2)
corr3 <- round(cor(classif.TRUN[,2:10]), 2)

  corrplot.mixed(corr)
  corrplot.mixed(corr2)
  corrplot.mixed(corr3)
  
  # NOTE: As they are correlated, the VSURF package is used for feature 
  # selection to confirm random forest feature selection.

set.seed(201903)

fs.I <- VSURF(classif.allclasses[,2:10], 
             classif.allclasses[,1], 
             clusterType = "FORK", 
             ntree = 500,mtry = 4)

fs.II <- VSURF(classif.mixSHD[,2:10], 
               classif.mixSHD[,1], 
              clusterType = "FORK", 
              ntree = 500,mtry = 4)

fs.III <- VSURF(classif.TRUN[,2:10], 
                classif.TRUN[,1], 
              clusterType = "FORK", 
              ntree = 500,mtry = 4)

  # NOTE: Warning message can be ignored. set mtry = best model reported in 
  # RFdata_report.txt.

# Save feature selection results as R data object to not have to rerun if requ.

  saveRDS(fs.I, 'output/I_fs_allRFdata.rds')
  saveRDS(fs.II, 'output/II_fs_mixRFdata.rds')
  saveRDS(fs.III, 'output/III_fs_TRUNRFdata.rds')
  
# Use VSURF_table() to create nice feature selection output.

  varimp_I<- VSURF_table(fs.I, classif.allclasses[,-1])
  varimp_II<- VSURF_table(fs.II, classif.mixSHD[,-1])
  varimp_III<- VSURF_table(fs.III, classif.TRUN[,-1])
  
# Write VSURF_table() objects as .csv file for manuscript use.

  write.csv(varimp_I, 'output/varimp.all_I_RFdata.csv', row.names = TRUE)
  write.csv(varimp_II, 'output/varimp.mix_II_RFdata.csv', row.names = TRUE)
  write.csv(varimp_III, 'output/varimp.TRUN_III_RFdata.csv', row.names = TRUE)

# PART D------------------Creating a Disease Risk Map --------------------------

# Import raster file similar to initial raster but without any ground pixel.
  
  imgpred <-  brick("data/20180528_ortho_no_ground.tif")

# Convert DN into reflectance
  
  imgpred <- imgpred/65535
  
# Rename raster layers for clarity (must be identical with initial layer names)

  names(imgpred) <- c("blue", "green", "red", "re", "nir", "alpha")

    ndvi.p <- (imgpred$nir-imgpred$red)/(imgpred$nir+imgpred$red)
    rg.p <- imgpred$red/imgpred$green
    sipi.p <- (imgpred$nir-imgpred$blue)/(imgpred$nir-imgpred$red)
    ari.p <- (1/imgpred$green)-(1/imgpred$re)

      imgpred <- addLayer(imgpred, c(ndvi.p, rg.p, sipi.p, ari.p))

  names(imgpred) <- c("blue", "green", "red", "re", 
                    "nir", "alpha", "ndvi", "rg", "sipi", "ari")
  
# Remove alpha channel
  
  riskpre <- imgpred[[c(1,2,3,4,5,7,8,9,10)]]
    names(riskpre)
  
# Start predicting pixel values (risk map) based on random forest model
  
  riskpred <- predict(riskpre, RFdata$fit)

    plot(riskpred)

  currentDate <- Sys.Date()
  rstFileName <- paste("output/riskmap",currentDate,".tif",sep="")
  writeRaster(riskpred, 
            file=rstFileName, 
            format = "GTiff", 
            bylayer=TRUE, 
            overwrite=TRUE)

# Note: The exported risk map was further processed in QGIS to change class 
# colours and have it ready for publication





#END
