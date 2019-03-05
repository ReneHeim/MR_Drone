####
# The following code reproduces results of the article:
# 
# Heim, Wright, Carnegie, Taylor, Scarth and Oldeland, 2018. Using drones and 
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
#                    "roxygen2",
#                    "tictoc", 
#                    "tidyverse", 
#                    "caret",
#                    "e1071", 
#                    "gdata", 
#                    "hsdar", 
#                    "utils", 
#                    "magrittr", 
#                    "rasterVis",  
#                    "rmarkdown"))

library(rgdal) #load installed pkgs
library(corrplot)
library(raster)
library(ggplot2)
library(tictoc)
library(caret)
#library(gdata)
library(reshape2)
library(rasterVis)
#library(magrittr)
#library(knitr)
library(cowplot)
library(VSURF)
#library(dplyr)


source("R/FUN_drop_cat_var.R")#drops factor and factor level
source("R/FUN_extract_polyclass.R")
source("R/FUN_prepggwide2long.R")
source("R/FUN_VSURF_table.R")

# PART A------------------------Cleaning Data-----------------------------------

# Loading image brick (1), convert to reflectance (2) and calculate indices (3).
# Also the band names were renamed for clarity according to camera manufacturer
# specifications (4).

img <-
  brick("data/20180528_ortho_ground.tif") #(1)

img <- img/65535 #(2)

  names(img) <- c("blue", "green", "red", "re", 
                "nir", "alpha")
  
    ndvi <- (img$nir-img$red)/(img$nir+img$red)
    rg <- img$red/img$green
    sipi <- (img$nir-img$blue)/(img$nir-img$red)
    ari <- (1/img$green)-(1/img$re)

      img <- addLayer(img, c(ndvi, rg, sipi, ari)) #(3)

  names(img) <- c("blue", "green", "red", "re", 
                "nir", "alpha", "ndvi", "rg", "sipi", "ari") #(4)

# Loading QGIS shape file (5) where sample polygons have been defined and
# extracting pixel from sample polygons (6). Then check data if everything 
# worked (7)

allclasses <- shapefile("data/samplepolygons.shp") #5

  unique(allclasses@data$Class) # To check what classes are contained

    allclasses.df<- extract_polyclass(allclasses, img)

  unique(allclasses.df$Class) #Check if all classes were imported

# Check for NAs

apply(allclasses.df, 2, function(x) any(is.na(x)))

  #full.df.wona <- na.omit(full.df) If NAs contained, they can be removed.

summary(allclasses.df$Class)
  #summary(full.df.wona$Class)
  #summary(full.df$Class)-summary(full.df.wona$Class)

allclasses.noalpha.df <- 
  allclasses.df[ , -which(names(allclasses.df) %in% c("alpha"))] #remove alpha



  write.csv(allclasses.noalpha.df, 
            'output/droneclassif_plusSHD.csv', 
            row.names = FALSE) #21
  
  classif.allclasses <- read.csv("output/droneclassif_plusSHD.csv", 
                                 check.names = FALSE) #22

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
  SHDmixdata.df[ , -which(names(SHDmixdata.df) %in% c("alpha"))] #19

write.csv(SHDmixdata.noalpha.df, 
          'output/droneclassif_mixSHD.csv', 
          row.names = FALSE) #21

classif.mixSHD <- read.csv("output/droneclassif_mixSHD.csv", 
                     check.names = FALSE) #22

# It is tested whether the bands in the classification data are correlated (23) 
# to later select a suitable selection method for relevant classification bands.

corr <- cor(classif.allclasses[,2:10]) #23
corr2 <- cor(classif.mixSHD[,2:10])

  round(corr,2)
  round(corr2,2)


corrplot(corr, type = "upper", 
         order = "hclust", 
         tl.col = "black", 
         tl.srt = 45, 
         sig.level = 0.01, 
         insig = "p-value")

corrplot.mixed(corr)
corrplot.mixed(corr2)

# Finally, re-assign (24) and rename (25) column names of the just exported df 
# to prepare (26) the data for a ggplot2 output (wide to long format). Plot (27)
# spectra to build figure 2 in the manuscript.

fig2 <- classif.allclasses[,1:6] #24

names(fig2) <- c("Type","475", "560", "668", "717", "840") #25

spectragg <- prep_gg(fig2, agg = TRUE) #26

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
  

ggsave("output/Figure2.allspectra.png",
       plot = C,
       width = 40,
       height = 20,
       units = "cm",
       dpi = 400
)


# PART B------------------Random Forest Classification--------------------------

#########################

#All classes (TR, UN, SHD)

# First we set a seed (30) to consistently reproduce the results of the 
# classification. Then we partition the extracted pixel data (classif) into a 
# training and test subset (31). Finally, we tune the settings for the random forest 
# training process (32, 33).

#set.seed(20180524) #30

inTraining <- createDataPartition(classif.allclasses$Class, p = .75, list = FALSE)
train <- classif.allclasses[ inTraining,]
test  <- classif.allclasses[-inTraining,] #31

rfControl <- trainControl(
  method = "boot",
  number = 100
) #32

rfGrid <- expand.grid(mtry = seq(1, ncol(train)-1, 1)) #33

# Then we initialize a timer (34) to measure the computation time of the 
# training (35) and prediction (36) process. Then, we stop timing (37). We can 
# assign (38) and export (39) a report of our classification. We are also 
# saving (40) the model object for later re-use (41) which will save time when 
# re-running the analysis.

tic("RF.I") #34

rfFit <- train(Class ~ ., data = train,
               method = "rf",
               importance = TRUE, ntree=500,
               trControl = rfControl, tuneGrid = rfGrid,
               metric = "Accuracy", maximize = TRUE) #35

rfPred <- 
  predict.train(rfFit, test[, !names(test) %in% c("Class")], type = "raw") #36

toc() #37

RFdataall <- 
  list(fit = rfFit,
       pred = predict.train(rfFit, test[, !names(test) %in% c("Class")], type = "raw"),
       confusion = confusionMatrix(rfPred, test$Class),
       varImp = varImp(rfFit, scale = TRUE)) #38

sink(file = 'output/I_allclasses_RFclassif_report.txt')
RFdataall #39
sink()

saveRDS(RFdataall, 'output/I_allclasses_RFclassif_object.rds') #40
RFdataall <- readRDS("output/I_allclasses_RFclassif_object.rds") #41
#########################

#Mix shadow (TR_S, UN_S)

# First we set a seed (30) to consistently reproduce the results of the 
# classification. Then we partition the extracted pixel data (classif) into a 
# training and test subset (31). Finally, we tune the settings for the random forest 
# training process (32, 33).

#set.seed(20180524) #30

inTraining <- createDataPartition(classif.mixSHD$Class, p = .75, list = FALSE)
train <- classif.mixSHD[ inTraining,]
test  <- classif.mixSHD[-inTraining,] #31

rfControl <- trainControl(
  method = "boot",
  number = 100
) #32

rfGrid <- expand.grid(mtry = seq(1, ncol(train)-1, 1)) #33

# Then we initialize a timer (34) to measure the computation time of the 
# training (35) and prediction (36) process. Then, we stop timing (37). We can 
# assign (38) and export (39) a report of our classification. We are also 
# saving (40) the model object for later re-use (41) which will save time when 
# re-running the analysis.

tic("RF.I") #34

rfFit <- train(Class ~ ., data = train,
               method = "rf",
               importance = TRUE, ntree=500,
               trControl = rfControl, tuneGrid = rfGrid,
               metric = "Accuracy", maximize = TRUE) #35

rfPred <- 
  predict.train(rfFit, test[, !names(test) %in% c("Class")], type = "raw") #36

toc() #37

RFdatamix <- 
  list(fit = rfFit,
       pred = predict.train(rfFit, test[, !names(test) %in% c("Class")], type = "raw"),
       confusion = confusionMatrix(rfPred, test$Class),
       varImp = varImp(rfFit, scale = TRUE)) #38

sink(file = 'output/II_mix_RFclassif_report.txt')
RFdatamix #39
sink()

saveRDS(RFdata, 'output/II_mix_RFclassif_object.rds') #40
RFdatamix <- readRDS("output/II_mix_RFclassif_object.rds") #41

#########################

#Only UN and TR

classif.TRUN <- drop_class(classif.allclasses, 
                           classif.allclasses$Class, 
                           "SHD")

inTraining <- createDataPartition(classif.TRUN$Class, p = .75, list = FALSE)
train <- classif.TRUN[ inTraining,]
test  <- classif.TRUN[-inTraining,] #31

rfControl <- trainControl(
  method = "boot",
  number = 100
) #32

rfGrid <- expand.grid(mtry = seq(1, ncol(train)-1, 1)) #33

# Then we initialize a timer (34) to measure the computation time of the 
# training (35) and prediction (36) process. Then, we stop timing (37). We can 
# assign (38) and export (39) a report of our classification. We are also 
# saving (40) the model object for later re-use (41) which will save time when 
# re-running the analysis.

tic("RF.I") #34

rfFit <- train(Class ~ ., data = train,
               method = "rf",
               importance = TRUE, ntree=500,
               trControl = rfControl, tuneGrid = rfGrid,
               metric = "Accuracy", maximize = TRUE) #35

rfPred <- 
  predict.train(rfFit, test[, !names(test) %in% c("Class")], type = "raw") #36

toc() #37

RFdata_TRUN <- 
  list(fit = rfFit,
       pred = predict.train(rfFit, test[, !names(test) %in% c("Class")], type = "raw"),
       confusion = confusionMatrix(rfPred, test$Class),
       varImp = varImp(rfFit, scale = TRUE)) #38

sink(file = 'output/III_onlyTRUN_RFclassif_report.txt')
RFdata_TRUN #39
sink()

saveRDS(RFdata_TRUN, 'output/III_onlyTRUN_RFclassif_object.rds') #40
RFdata_TRUN <- readRDS("output/III_onlyTRUN_RFclassif_object.rds") #41

# PART C----------------------Feature Selection---------------------------------

# Here (42) we select relevant classification features for our dataset.
# Also, we export the object containing the features (43), create a relative version 
# of ranked features (44), an absolute version (45), combine them (46), 
# rename these objects (47) and export them as a table (48). VSURF can select 
# relevant features from datasets that contain correlated predictor variables.

#set.seed(20190301)

fs.I <- VSURF(classif.allclasses[,2:10], 
             classif.allclasses[,1], 
             clusterType = "FORK", 
             ntree = 500,mtry = 4) #42

fs.II <- VSURF(classif.mixSHD[,2:10], 
               classif.mixSHD[,1], 
              clusterType = "FORK", 
              ntree = 500,mtry = 4) #42

fs.III <- VSURF(classif.TRUN[,2:10], 
                classif.TRUN[,1], 
              clusterType = "FORK", 
              ntree = 500,mtry = 4) #42

# Warning message can be ignored. set mtry = best model reported in RFdata_report.txt


saveRDS(fs.I, 'output/I_fs_allRFdata.rds') #43
saveRDS(fs.II, 'output/II_fs_mixRFdata.rds')
saveRDS(fs.III, 'output/III_fs_TRUNRFdata.rds')


  varimp_I<- VSURF_table(fs.I, classif.allclasses[,-1])
  varimp_II<- VSURF_table(fs.II, classif.mixSHD[,-1])
  varimp_III<- VSURF_table(fs.III, classif.TRUN[,-1])

write.csv(varimp_I, 'output/varimp.all_I_RFdata.csv', row.names = TRUE)
write.csv(varimp_II, 'output/varimp.mix_II_RFdata.csv', row.names = TRUE)
write.csv(varimp_III, 'output/varimp.TRUN_III_RFdata.csv', row.names = TRUE)

# PART D------------------Creating a Disease Risk Map --------------------------

# Import aerial scene as brick image (49), similar to the one imported initially
# (see #1). Any ground/soil pixel has been removed (Agisoft Photoscan Pro) to 
# have only the lemon myrtle trees remaining. The imported image can be used to 
# test the model by predicting the class of each pixel in this images. To display 
# the image, we set background pixel as being tranparent (50) by removing the 
# alpha channel which was not used as a predictor variable. Then we divide by 
# the max DN value to yield reflectance values as before (51) and rename (52) 
# the bands to match the data on which we trained the random forest model. Now 
# we can predict the pixel-classes (53), plot (54) and export the result (55).

imgpred <-  brick("data/20180528_ortho_no_ground.tif") #49

imgpred <- imgpred/65535 #51

ndvi.p <- 
  (imgpred$X20180528_ortho_no_ground.5-imgpred$X20180528_ortho_no_ground.3)/
  (imgpred$X20180528_ortho_no_ground.5+imgpred$X20180528_ortho_no_ground.3)
rg.p <- 
  imgpred$X20180528_ortho_no_ground.2/imgpred$X20180528_ortho_no_ground.3
sipi.p <- 
  (imgpred$X20180528_ortho_no_ground.5-imgpred$X20180528_ortho_no_ground.1)/
  (imgpred$X20180528_ortho_no_ground.5-imgpred$X20180528_ortho_no_ground.3)
ari.p <- 
  (1/imgpred$X20180528_ortho_no_ground.2)-(1/imgpred$X20180528_ortho_no_ground.4)

imgpred <- addLayer(imgpred, c(ndvi.p, cal.p, sipi.p, ari.p))

names(imgpred) <- c("blue", "green", "red", "re", 
                    "nir", "alpha", "ndvi", "rg", "sipi", "ari") #52

riskpre <- imgpred[[c(1,2,3,4,5,7,8,9,10)]] #50 remove alpha channel

riskpred <- predict(riskpre, RFdata$fit) #53

plot(riskpred) #54

currentDate <- Sys.Date()
rstFileName <- paste("output/riskmap",currentDate,".tif",sep="")
writeRaster(riskpred, 
            file=rstFileName, 
            format = "GTiff", 
            bylayer=TRUE, 
            overwrite=TRUE) #55

# Note: The exported risk map was further processed in QGIS to change class 
# colours and have it ready for publication





#END
