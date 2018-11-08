####
# The following code reproduces results of the article:
# 
# Heim, Wright and Oldeland, 2018. Using drones and multispectral imagery to 
#detect myrtle rust on a lemon myrtle plantation. JOURNAL
#
# The code is split in two major sections according to our research questions.
#
# A)  Can we discriminate fungicide treated and untreated lemon myrtle trees 
#     based on multispectral imagery captured by a drone?
#       - Load GIS files to extract raw pixel data for all available spectral 
#         bands and pre-defined classes (TReated, UNtreated, SHAdow)
#       - Raw pixel data will be classified using a random forest classifier and 
#         all spectral bands.
#       - Accuracy metrics will be reported in an error matrix
#
# B)  Is a leaf-scale classification equally accurate as a canopy-scale 
#     classification?
#       - Resample hyperspectral leaf data to multispectral (Micasense) 
#         resolution to be able to compare two classifications based on the same 
#         sensor.
#       - Perform classifcation on this resampled data.
#       - Report accuracies and compare to results found in A) (DN vs refl.)
#
# C)  Are the important spectral regions found on the leaf scale similar to 
#     those found from the air?
#       - Compare relevant spectral regions between A) and B)
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
library(raster)
library(tictoc)
library(caret)
library(gdata)
library(hsdar)
library(utils)
library(rasterVis)
library(rmarkdown)
library(roxygen2)
library(magrittr)
library(knitr)
library(reshape2)
library(cowplot)
library(VSURF)
library(dplyr)

source("R/FUN_raw2speclibhsdar.R")#coverts spec data to hsdar lib
source("R/FUN_drop_cat_var.R")#drops factor and factor level
source("R/FUN_extract_pixel.R")
source("R/FUN_prepggwide2long.R")

# SECTION A --------------------------------------------------------------------

# Loading image brick (1) and 

img <-
  brick("data/20180528_ortho_ground.tif") 

img <- img/65535
  
ndvi <- (img$X20180528_ortho_ground.5-img$X20180528_ortho_ground.3)/(img$X20180528_ortho_ground.5+img$X20180528_ortho_ground.3)
rg <- img$X20180528_ortho_ground.2/img$X20180528_ortho_ground.3
sipi <- (img$X20180528_ortho_ground.5-img$X20180528_ortho_ground.1)/(img$X20180528_ortho_ground.5-img$X20180528_ortho_ground.3)
ari <- (1/img$X20180528_ortho_ground.2)-(1/img$X20180528_ortho_ground.4)
img <- addLayer(img, c(ndvi, rg, sipi, ari))

names(img) <- c("blue", "green", "red", "re", "nir", "alpha", "ndvi", "rg", "sipi", "ari")

# Loading QGIS shape file where sample polygons have been defined (3)
alldata <- shapefile("data/samplepolygons.shp") #3

# Extracting pixel from sample polygons (4) and check data (5)

dfAll <- extract_pixel(img, alldata) #4
head(dfAll) #5

# Use "Treatment" (5) and "ID" (6) which should have been defined in QGIS to 
# construct a new df (7).

Treatment <- alldata$Type #5
ID <- alldata$id #6

dfclass <- cbind(Treatment, ID) #7

# Look up ID and the assigned class (TR, UN, SHD) from the "Type" column (8) and
# then save them in a list (9) while looping through each ID (10).

classvec <- unique(as.character(dfclass[,1])) #8
li <- list() #9

for(i in classvec){
  x <- subset(dfclass, dfclass[,1] == i)
  li[[i]] <- assign(paste0(i,"_","num"), as.numeric(x[,2]))
} #10

# Writes a new df that has written the class labels (11, 12, 14) instead of ID 
# numbers. New df (14) has the "Type" column added (15).

for(i in UN_num){
  dfAll$class[dfAll$class == i] <- "UN" #11
}

for(i in TR_num){
  dfAll$class[dfAll$class == i] <- "TR" #12
}

for(i in SHD_num){
  dfAll$class[dfAll$class == i] <- "SHD" #13
}

dfAll$class <- as.factor(dfAll$class) #14
names(dfAll)[11] <- c("Type") #15

# Divide by 65535 (DN max value) to yield refl between 0 and 1 (16).rename band names according to Micasense RedEdge
# specifications (2). #Band1(blue=475nm), Band2 (green=560nm), Band3 (Red=668nm), 
# Band4 (840nm), Band5 (rededge=717nm), Band6 (alphaband)


dfAll <- dfAll[,c(11,1,2,3,4,5,7,8,9,10)]
names(dfAll) <- c("Type","blue", "green", "red", "re", "nir", "ndvi", "rg", "sipi", "ari")

# Save (17) and reload (18) df containing reflectance values and classes stored 
# in the "Type" column. Also remove alpha (transparency) channel (19) as it is 
# not used for classification. Reorder columns (20) and check if all classes are 
# present (21).

write.csv(dfAll, 'output/2018MyrtleRust_Refl.csv', row.names = FALSE) #17
classif <- read.csv("output/2018MyrtleRust_Refl.csv", check.names = FALSE) #18

# mycor <- cor(classif[,-6]) #test correlation - all bands correlated



# Finally, re-assign (22) and rename (23) column names of the just exported df 
# to prepare (24) the data for a ggplot2 output (wide to long format). Plot (25)
# spectra to build Figure 2B(C).

gplot <- classif[,1:6] #22

names(gplot) <- c("Type","475", "560", "668", "717", "840") #23

spectragg <- prep_gg(gplot, agg = TRUE) #24

C <- ggplot(spectragg, aes(Wavelength, Reflectance*100, colour = Type)) +
  geom_line(aes(linetype=Type), size = 1)+
  geom_point(aes(shape=Type), size = 2)+
  scale_color_manual(values=c("#525252", "#00CD00", "#cd5c00"))+
  labs(x = "", y = "")+
  theme_minimal(base_size=20)+
  theme(legend.title = element_blank(),
        legend.key.width = unit(1.5, "cm"))#25

# Now we drop the "shadow" class (26), rename the columns (27), convert wide to 
# long format (28) and create Figure 2B(C-S) (29).

dropSHD <- drop_class(classif, classif$Type, "SHD")
dropSHD <- dropSHD[,1:6] #26
names(dropSHD) <- c("Type","475", "560", "668", "717", "840") #27

spectraggdropSHD <- prep_gg(dropSHD, agg = TRUE) #28

C_S <- ggplot(spectraggdropSHD, aes(Wavelength, Reflectance*100, colour = Type)) +
  geom_line(aes(linetype=Type), size = 1)+
  geom_point(aes(shape=Type), size = 2)+
  scale_color_manual(values=c("#00CD00", "#cd5c00"))+
  labs(x = "", y = "")+
  theme_minimal(base_size=20)+
  theme(legend.title = element_blank(),
        legend.key.width = unit(1.5, "cm"))#29

# RF (all classes)--------------------------------------------------------------

# First we set a seed (30) to avoid random number generation in vulnerable 
# processes. Then we partition the extracted pixel data in a training and 
# test subset (31). Finally, we tune the settings for the random forest 
# training process (32, 33).

#set.seed(20180524) #30

inTraining <- createDataPartition(classif$Type, p = .75, list = FALSE)
train <- classif[ inTraining,]
test  <- classif[-inTraining,] #31

rfControl <- trainControl(
  method = "boot",
  number = 100
) #32

rfGrid <- expand.grid(mtry = seq(1, ncol(train)-1, 1)) #33

# Then we initialize a timer (34) to measure the computation time of the 
# training (35) and prediction (36) process. We stop timing (37). Then, we can 
# assign (38) and export (39) a report of our classification. We are also 
# saving (40) the model object to reload (41) and save time when re-running 
# the code.

tic("RF.I") #34

rfFit <- train(Type ~ ., data = train,
               method = "rf",
               importance = TRUE, ntree=500,
               trControl = rfControl, tuneGrid = rfGrid,
               metric = "Accuracy", maximize = TRUE) #35

rfPred <- 
  predict.train(rfFit, test[, !names(test) %in% c("Type")], type = "raw") #36

toc() #37

canopydata <- 
  list(fit = rfFit,
       pred = predict.train(rfFit, test[, !names(test) %in% c("Type")], type = "raw"),
       confusion = confusionMatrix(rfPred, test$Type),
       varImp = varImp(rfFit, scale = TRUE)) #38

sink(file = 'output/I_canopydata_report.txt')
canopydata #39
sink()

saveRDS(canopydata, 'output/I_canopydata_object.rds') #40
canopydata <- readRDS("output/I_canopydata_object.rds") #41

# Feature selection canopydata

# Here (42) we select relevant classification features for the canopy dataset.
# Also, we export the object containing the features (43), create a relative version 
# of ranked features (44), an absolute version (45), combine them (46), 
# rename these objects (47) and export them as a table (48).

fs.1 <- VSURF(classif[,2:10], 
             classif[,1], 
             clusterType = "FORK", 
             ntree = 500,mtry = 4) #42

#warning can be ignored. set mtry = best model reported in canopydata_report.txt


saveRDS(fs.1, 'output/I_fs_canopy.rds') #43

fs.norm.1 <- (fs.1$imp.varselect.thres-min(fs.1$imp.varselect.thres))/
  (max(fs.1$imp.varselect.thres)-min(fs.1$imp.varselect.thres)) #44

out <- classif[,-1]
varimp.canopy <- rbind(names(out[,as.numeric(fs.1$varselect.thres)]),
              round(fs.1$imp.varselect.thres,2), #45
              round(fs.norm.1, 2)) #46

row.names(varimp.canopy) <- c('Predictor', 'Abs. Imp.', 'Rel. Imp.') #47


write.csv(varimp.canopy, 'output/varimp.canopy.csv', row.names = TRUE) #48

# RiskMap ----------------------------------------------------------------------

# Import a brick image (49), similar to the one imported initially (see #1). Any
# ground pixel has been removed (Agisoft Photoscan Pro) to have only the lemon
# myrtle trees available. The image can be used to test the model by predicted 
# the class of each pixel in this images. To display the image, we set 
# background pixel as being tranparent (50) by removing the alpha channel
# which was not used as a predictor variable, divide by the max DN value
# (51) and rename (52) the bands to match the data on which we trained the 
# random forest model. Now we can predict the pixel-classes (53) and plot (54)
# and export the result (55).

imgpred <-  brick("data/20180528_ortho_no_ground.tif") #49

imgpred <- imgpred/65535 #51

ndvi.p <- (imgpred$X20180528_ortho_no_ground.5-imgpred$X20180528_ortho_no_ground.3)/(imgpred$X20180528_ortho_no_ground.5+imgpred$X20180528_ortho_no_ground.3)
rg.p <- imgpred$X20180528_ortho_no_ground.2/imgpred$X20180528_ortho_no_ground.3
sipi.p <- (imgpred$X20180528_ortho_no_ground.5-imgpred$X20180528_ortho_no_ground.1)/(imgpred$X20180528_ortho_no_ground.5-imgpred$X20180528_ortho_no_ground.3)
ari.p <- (1/imgpred$X20180528_ortho_no_ground.2)-(1/imgpred$X20180528_ortho_no_ground.4)

imgpred <- addLayer(imgpred, c(ndvi.p, cal.p, sipi.p, ari.p))

names(imgpred) <- c("blue", "green", "red", "re", "nir", "alpha", "ndvi", "rg", "sipi", "ari") #52

riskpre <- imgpred[[c(1,2,3,4,5,7,8,9,10)]] #50 remove alpha channel for transparancy



riskpred <- predict(riskpre, canopydata$fit) #53

plot(riskpred) #54

currentDate <- Sys.Date()
rstFileName <- paste("output/riskmap",currentDate,".tif",sep="")
writeRaster(riskpred, file=rstFileName, format = "GTiff", bylayer=TRUE, overwrite=TRUE) #55

# Note: The exported risk map was further processed in QGIS to change class colours
# and have it publication ready

# SECTION B --------------------------------------------------------------------

# In this section, we code the analysis to compare leaf-scale vs canopy scale 
# classification results. Therefore, we import leaf spectral reflectance
# data (56) that was recorded on the same plants but in a previous study. Then,
# we drop (57) the class "Healthy" as this was not recorded on the plantation 
# where the aerial imagery was captured. We create a spectral library (58) to 
# make use of the hsdar package and then resample the hyperspectral data to 
# the Micasense RedEdge band specifications (59-62).

hypdata <- 
  read.csv('data/data.wo.out.binned.cut.csv', check.names = FALSE) %>% #56
  drop_class(., .$Type, "Healthy") %>% 
  .[,1:36]#57

speclib2 <- raw2speclib(hypdata) #58

center <-  c(475, 560, 668, 717, 840) #59
fwhm <- c(20, 20, 10, 10, 40) #60
micasense <- as.data.frame(cbind(center, fwhm)) #61

data_mica <- spectralResampling(speclib2, micasense) #62

# Now we would like to use the resampled data (leaf scale = L) for classification to 
# then compare it with the data that was captured on the canopy level (Micasense
# RedEdge). We extract the data from our spectral library (65) and add a "Type"
# column to assign the classes (66). We also rename the columns according to the 
# Micasense RedEdge bands (67).

micadata <- as.data.frame(data_mica@spectra@spectra_ma) #65
micadata <- cbind('Type'=hypdata$Type, micadata) #66

names(micadata) <- c("Type", "blue", "green", "red", "re", "nir") #67

# We plot the hyperspectral leaf data (68) to check if the export worked and
# then plot the resampled version (69) to include it in Figure 2B (L).

spectraggII <- prep_gg(hypdata, agg = TRUE)

spectraggII <- spectraggII %>% 
  mutate(Type = recode(Type, `Treated` = "TR", `Untreated` = "UN"))

hyperplot <- ggplot(spectraggII, aes(Wavelength, Reflectance, colour = Type)) +
  geom_line(aes(linetype=Type), size = 1)+
  geom_point(aes(shape=Type), size = 2)+
  scale_color_manual(values=c("#00CD00", "#cd5c00"))+
  labs(x = "Wavelength [nm]", y = "Reflectance [%]")#68

mplot <- micadata
names(mplot) <- c("Type", "475", "560", "668", "717", "840")

spectraggIII <- prep_gg(mplot, agg = TRUE)

spectraggIII <- spectraggIII %>% 
  mutate(Type = recode(Type, `Treated` = "TR", `Untreated` = "UN"))

L <- ggplot(spectraggIII, aes(Wavelength, Reflectance, colour = Type)) +
  geom_line(aes(linetype=Type), size = 1)+
  geom_point(aes(shape=Type), size = 2)+
  scale_color_manual(values=c("#00CD00", "#cd5c00"))+
  labs(x = "", y = "")+
  theme_minimal(base_size=20)+
  theme(legend.title = element_blank(),
        legend.key.width = unit(1.5, "cm"))#69

# Create spectral vegetation indices (SVI) to be included in the classification

rgsvi <- '560/668'
sipisvi <- '(840-475)/(840-668)'
arisvi <- '(1/560)-(1/717)'

SIPI2 <- vegindex(speclib2, sipisvi)
NDVI2 <- vegindex(speclib2, "NDVI")
RG <- vegindex(speclib2, rgsvi)
ARI2 <- vegindex(speclib2, arisvi)

micadata.ind <- cbind(micadata, NDVI2, SIPI2, RG, ARI2)


# RF (Leaf) --------------------------------------------------------------------

# The seed from our previous classification is still active (#30). We partition 
# the resampled data (70) in a training (71) and test (72) subset. We tune the 
# settings for the random forest training process (73, 73).


inTrainingM <- createDataPartition(micadata.ind$Type, p = .75, list = FALSE) #70
trainM <- micadata.ind[ inTrainingM,] #71
testM  <- micadata.ind[-inTrainingM,] #72

rfControl <- trainControl(
  method = "boot",
  number = 100
) #73

rfGrid <- expand.grid(mtry = seq(1, ncol(trainM)-1, 1)) #74

# Again, we initialize a timer (74) to measure the computation time of the 
# training (75) and prediction (76) process. We stop timing (77). Then, we can 
# assign (78) and export (79) a report of our classification. We are also 
# saving (80) the model object to reload (81) and save time when re-running 
# the code.

tic("RF.II") #74

rfFit.M <- train(Type ~ ., data = trainM,
                 method = "rf",
                 importance = TRUE, ntree=500,
                 trControl = rfControl, tuneGrid = rfGrid,
                 metric = "Accuracy", maximize = TRUE) #75

rfPred.M <- #76
  predict.train(rfFit.M, testM[, !names(testM) %in% c("Type")], type = "raw")

toc() #77

leafdata <- 
  list(fit = rfFit.M,
       pred = predict.train(rfFit.M, testM[, !names(testM) %in% c("Type")], type = "raw"),
       confusion = confusionMatrix(rfPred.M, testM$Type),
       varImp = varImp(rfFit.M, scale = TRUE)) #78

sink(file = 'output/II_leafdata_report.txt')
leafdata #79
sink()

saveRDS(leafdata, 'output/II_leafdata_object.rds') #80
#leafdata <- readRDS("output/II_leafdata_object.rds") #81

# Feature selection leafdata

fs.2 <- VSURF(micadata.ind[,2:10], 
              micadata[,1], 
              clusterType = "FORK", 
              ntree = 500,mtry = 4)

saveRDS(fs.2, 'output/II_fs_leafdata.rds')

fs.norm.2 <- (fs.2$imp.varselect.thres-min(fs.2$imp.varselect.thres))/
  (max(fs.2$imp.varselect.thres)-min(fs.2$imp.varselect.thres))

out2 <- micadata.ind[,-1]
varimp.leaf <- rbind(names(out2[,as.numeric(fs.2$varselect.thres)]),
              round(fs.2$imp.varselect.thres,2), 
              round(fs.norm.2, 2))

row.names(varimp.leaf) <- c('Predictor', 'Abs. Imp.', 'Rel. Imp.')


write.csv(varimp.leaf, 'output/varimp.leaf.csv', row.names = TRUE)

# RF (Canopy without class "Shadow") -------------------------------------------

# Now, we have to repeat the classification using the aerial imagery again. This
# time, we drop the class 'SHADOW' (82). This class was not collected for the 
# leaf spectral data and we want to have a fair comparison between leaf-scale 
# and canopy-scale. We can test if the factor levels were dropped (83).

classif2 <- drop_class(classif, classif$Type, 'SHD') #82
unique(classif2$Type) #83


# We can repeat the data partition (84-86), tune the random forest settings 
# (87-89) and then start the training (90) and prediction (91) process. We 
# export (92) a report of our classification. We are also saving (93) the model 
# object to reload (94) if necessary.

inTrainingD <- createDataPartition(classif2$Type, p = .75, list = FALSE) #84
trainD <- classif2[ inTrainingD,] #85
testD  <- classif2[-inTrainingD,] #86

rfControl <- trainControl(
  method = "boot",
  number = 100
) #87

rfGrid <- expand.grid(mtry = seq(1, ncol(trainD)-1, 1)) #88

tic("RF.III")

rfFit.D <- train(Type ~ ., data = trainD,
                 method = "rf",
                 importance = TRUE, ntree=500,
                 trControl = rfControl, tuneGrid = rfGrid,
                 metric = "Accuracy", maximize = TRUE) #89

rfPred.D <- #90
  predict.train(rfFit.D, testD[, !names(testD) %in% c("Type")], type = "raw")

toc()

canopywoshddata <- #91
  list(fit = rfFit.D,
       pred = predict.train(rfFit.D, testD[, !names(testD) %in% c("Type")], type = "raw"),
       confusion = confusionMatrix(rfPred.D, testD$Type),
       varImp = varImp(rfFit.D, scale = TRUE))

sink(file = 'output/III_canopywoshddata_report.txt')
canopywoshddata #92
sink()

saveRDS(canopywoshddata, 'output/III_canopywoshddata_object.rds') #93
#canopywoshddata <- readRDS("output/III_canopywoshddata_object.rds") #94

# Feature selection canopy data without shadow class

fs.3 <- VSURF(classif2[,2:10], 
              classif2[,1], 
              clusterType = "FORK", 
              ntree = 500,mtry = 4)

saveRDS(fs.3, 'output/III_fs_canopy2.rds')

fs.norm.3 <- (fs.3$imp.varselect.thres-min(fs.3$imp.varselect.thres))/
  (max(fs.3$imp.varselect.thres)-min(fs.3$imp.varselect.thres))

out3 <- classif2[,-1]
varimp.canopywoshd <- rbind(names(out3[,as.numeric(fs.3$varselect.thres)]),
               round(fs.3$imp.varselect.thres,2), 
               round(fs.norm.3, 2))

row.names(varimp.canopywoshd) <- c('Predictor', 'Abs. Imp.', 'Rel. Imp.')


write.csv(varimp.canopywoshd, 'output/varimp.canopywoshd.csv', row.names = TRUE)

# Plot all multispectral signatures to yield Figure 2B

res <- ggdraw() +
  draw_plot(C, x = 0, y = .66, width = 1, height = .33)+
  draw_plot(C_S, x = 0, y = .33, width = 1, height = .33) +
  draw_plot(L, x = 0, y = 0, width = 1, height = .33) +
  draw_plot_label(label = c("", "", ""), size = 18,
                  x = c(.12,.12,.12), y = c(0.99,0.66,0.33))

ggsave("output/Figure2.allspectra.png",
       plot = res,
       width = 40,
       height = 20,
       units = "cm",
       dpi = 400
)

#END
