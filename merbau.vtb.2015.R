# Last run June 9, 2015

rm(list = ls(all=T))
setwd("D:/Github/BGE-SDM")
# save(list=ls(all=TRUE), file="D:/Papers.Projects/Merbau.VTB/R/merbau.vtb.2015.RData") # save RDATA for later use
load("D:/Papers.Projects/Merbau.VTB/R/merbau.vtb.2015.RData")

library(raster)
library(rgdal)
library(dismo)
library(maptools)
library(SDMTools)
library(rgbif)
library(rgeos)
library(adehabitatHS)
#library(rJava); # Sys.setenv(JAVA_HOME='C:\\Program Files\\Java\\jre7') # for 64-bit version
source("D:/R/Scripts/VIF.R")
source("D:/R/Scripts/null_model_function.R")

countries <- rgdal::readOGR("D:/Github/BGE-SDM/GISDATA/Study Area SHP")
proj4string(countries) <- "+proj=longlat +datum=WGS84"
plot(countries)
str(countries)

native.extent <- extent(-43, 109, 25, 81)

### 1. Merbau - Intsia bijuga ####

?name_backbone
key <- name_backbone(name='Intsia bijuga', rank='species')$usageKey
# key.genus <- name_backbone(name='Intsia', rank='genus')$usageKey
key
?occ_search
dat <- occ_search(taxonKey=key, return='data', limit=10000) # limit=300, hasCoordinate = T, spatialIssues=F, basisOfRecord='PRESERVED_SPECIMEN'
# dat.genus <- occ_search(taxonKey=key.genus, return='data', limit=10000) # limit=300, hasCoordinate = T, spatialIssues=F, basisOfRecord='PRESERVED_SPECIMEN'
head(dat); dim(dat) # 1901  110 - 9 June 2015
# head(dat.genus); dim(dat.genus) # 1873  110


unique(dat$basisOfRecord)
dat <- dat[dat$basisOfRecord == "PRESERVED_SPECIMEN",]
head(dat); dim(dat) # 1583  110

dat[,c('decimalLongitude', 'decimalLatitude')]
dat <- dat[!dat$decimalLongitude == 0 & !dat$decimalLatitude == 0, ]
head(dat); dim(dat) # 1161  110
dat <- dat[!(is.na(dat$decimalLongitude)) & !(is.na(dat$decimalLatitude)),]
head(dat); dim(dat) # 529  110

gbifmap(input=dat)
# gbifmap(input=dat.genus)

plot(countries); points(dat$decimalLongitude, dat$decimalLatitude, pch=19, col="red")

merbau.gbif <- dat[, c('species', 'decimalLongitude', 'decimalLatitude')]
head(merbau.gbif); dim(merbau.gbif) # 529 3
duplicates <- duplicated(merbau.gbif)
merbau.gbif <- merbau.gbif[!duplicates,]
head(merbau.gbif); dim(merbau.gbif) # 330 3
# merbau.gbif <- na.omit(merbau.gbif)
write.csv(merbau.gbif, 'D:/Papers.Projects/Merbau.VTB/R/merbau.gbif.unique.csv', row.names=F)

# Select native range

merbau.gbif.native <- subset(merbau.gbif, merbau.gbif$decimalLatitude >= native.extent@ymin & merbau.gbif$decimalLatitude <= native.extent@ymax & merbau.gbif$decimalLongitude >= native.extent@xmin & merbau.gbif$decimalLongitude <= native.extent@xmax)
head(merbau.gbif.native); dim(merbau.gbif.native) # 173 3
write.csv(merbau.gbif.native, 'D:/Papers.Projects/Merbau.VTB/R/merbau.gbif.native.unique.csv', row.names=F)

plot(countries); points(merbau.gbif.native$decimalLongitude, merbau.gbif.native$decimalLatitude, pch=19, col="purple")

coordinates(merbau.gbif.native) <- ~decimalLongitude+decimalLatitude
str(merbau.gbif.native)
merbau.gbif.native@proj4string <- P4S.latlon
plot(countries); plot(merbau.gbif.native, col='purple', add=T)

### 2. Compile climate data ####

files.present <- list.files('D:/GIS/Worldclim/Present/5arcmin/bio/', pattern="[.]bil$", full.names=T) # alternatives for pattern (c|C)(e|E)(l|L)$
# files <- list.files('Z:/World/Climate/Worldclim/05arcmin/Present/bio/', pattern="[.]bil$", full.names=T) # alternatives for pattern (c|C)(e|E)(l|L)$
files.present
present.stack <- stack(files.present[2:20]) # remove alt
head(present.stack)
present.df <- as.data.frame(present.stack, xy=T)
coordinates(present.df) <- ~x+y
gridded(present.df) <- T
present.df@proj4string <- P4S.latlon
present.df$grid.index <- present.df@grid.index # Add grid.index value
head(present.df)
image(present.df, 'bio01')

### 3. Get abiotic bioclim data and remove duplicates - merbau.native.unique ####

str(merbau.gbif.native); str(present.df)
merbau.gbif.native.abiotic <- over(merbau.gbif.native, present.df) # Get climate variables + grid.index
str(merbau.gbif.native.abiotic) # 173  20
head(merbau.gbif.native.abiotic)
dim(merbau.gbif.native.abiotic); dim(merbau.gbif.native) # 173  20  1
head(merbau.gbif.native); str(merbau.gbif.native)
merbau.gbif.native <- cbind(merbau.gbif.native, merbau.gbif.native.abiotic) # Link species col and climate data
head(merbau.gbif.native)
duplicates <- duplicated(merbau.gbif.native[,c("species", "grid.index")]) # Duplicates on grid.index
table(duplicates) # 114 F 59 T
merbau.gbif.native.unique <- merbau.gbif.native[!duplicates,] # remove duplicates
str(merbau.gbif.native.unique) # data.frame
head(merbau.gbif.native.unique); dim(merbau.gbif.native.unique) # 114 24
summary(merbau.gbif.native.unique)
merbau.gbif.native.unique <- na.omit(merbau.gbif.native.unique)
head(merbau.gbif.native.unique); dim(merbau.gbif.native.unique) # 103  24
summary(merbau.gbif.native.unique)
plot(raster(present.df, 'bio01')); points(merbau.gbif.native.unique$decimalLongitude, merbau.gbif.native.unique$decimalLatitude)

boxplot(merbau.gbif.native.unique$bio01, main = "Mean annual Temperature", ylab="Temperature x 10 (in °C)")
boxplot(merbau.gbif.native.unique$bio12, main = "Annual precipitation", ylab="Precipitation (in mm)")

?boxplot
library(ggplot2)

qplot(species, bio01, data = merbau.gbif.native.unique, geom="boxplot") + geom_point()
qplot(species, bio01, data = merbau.gbif.native.unique, geom="boxplot") + geom_jitter(position=position_jitter(w=0.05, h=0))
qplot(species, bio12, data = merbau.gbif.native.unique, geom="boxplot") + geom_jitter(position=position_jitter(w=0.05, h=0))

### 4. Run Maxent model with bioclim data in 500 km buffered area around presences to balance prevalence ####

# climate has priority over soil - hierarchical model, 1 climate, 2 soil
# Boucher-Lalonde, V., A. Morin and D. J. Currie (2012). "How are tree species distributed in climatic space? A simple and general pattern." Global Ecology and Biogeography 21(12): 1157-1166.

### Create empty mask layer

mask <- raster(files.present[1])
plot(mask)
mask <- !is.na(mask) # all values to 1
mask[mask == 0] <- NA # zero values to NA
plot(mask)
summary(mask)
writeRaster(mask, filename  = "D:/Papers.Projects/Merbau.VTB/R/mask.asc", format = 'ascii', NAflag = -9999, overwrite = T)
# mask <- raster('D:/Papers.Projects/Merbau.VTB/R/mask.asc')
plot(mask, col='red')

# add mask to present.df
present.df@data$mask <- as.data.frame(stack('D:/Papers.Projects/Merbau.VTB/R/mask.asc')) # Add mask layer to spdf
head(present.df); dim(present.df) # 7776000      21

head(merbau.gbif.native.unique)
merbau.gbif.native.unique.maxent <- merbau.gbif.native.unique[, c("species", "decimalLongitude", "decimalLatitude")]
names(merbau.gbif.native.unique.maxent) <- c("species", "lon", "lat")
head(merbau.gbif.native.unique.maxent); dim(merbau.gbif.native.unique.maxent) # 103   3
plot(countries); points(merbau.gbif.native.unique.maxent$lon, merbau.gbif.native.unique.maxent$lat, pch=19, col='red')
write.csv(merbau.gbif.native.unique.maxent, 'D:/Papers.Projects/Merbau.VTB/R/merbau.gbif.native.unique.maxent.csv', row.names=F) # write species points file

# Convert to coll.locs spatial Points Data Frame and create 500 km buffer
coordinates(merbau.gbif.native.unique.maxent) <- ~lon+lat
proj4string(merbau.gbif.native.unique.maxent) <- P4S.latlon
plot(countries); points(merbau.gbif.native.unique.maxent, pch=19, cex=0.5, col='red')
x <- circles(merbau.gbif.native.unique.maxent, d=500000, lonlat=TRUE) # 500 km
pol <- gUnaryUnion(x@polygons) # dissolve polygons
extent(pol)
plot(pol, col='blue', add=T); points(merbau.gbif.native.unique.maxent, pch=19, cex=0.5, col='red')

# extract cell numbers for the circles
v <- extract(mask, x@polygons, cellnumbers=T)
str(v)
# use rbind to combine the elements in list v
v <- do.call(rbind, v)
head(v); dim(v) # 211089  2

# remove ocean cells
v <- unique(na.omit(v))
head(v); dim(v) # 56783  2

# to display the results
m <- mask
m[] <- NA # empty mask
m[as.vector(v[,1])] <- 1
plot(m, col='purple')
extent(m)
str(m); summary(m)
plot(m, ext=extent(x@polygons)+1, col='blue')
plot(x@polygons, add=T)
points(merbau.gbif.native.unique.maxent, pch=19, cex=0.5, col='red')
plot(countries, add=T)
str(m) # rasterlayer

# Write mask for buffered areas
writeRaster(m, filename  = "D:/Papers.Projects/Merbau.VTB/R/mask.500km.asc", format = 'ascii', NAflag = -9999, overwrite = T)

### Add mask.buffer to present.df
head(present.df); dim(present.df) # 7776000      21
# present.df <- subset(present.df, select = -c(mask)) # remove mask column

present.species.df <- present.df # copy present.df
mask.buffer.df <- as.data.frame(stack('D:/Papers.Projects/Merbau.VTB/R/mask.500km.asc'), xy=T) # read mask layer
head(mask.buffer.df); dim(mask.buffer.df); colSums(mask.buffer.df, na.rm=T, dims=1)
head(present.species.df); str(present.species.df)
present.species.df$mask <- mask.buffer.df[,'mask.500km'] # replace mask with mask.500km
head(present.species.df); dim(present.species.df) # 7776000      21
str(present.species.df) # SpatialPixelsDataFrame
present.species.df <- na.omit(present.species.df@data)
head(present.species.df); dim(present.species.df) # 56783    21

### 5. Select uncorrelated variables using VIF ####

### Variance Inflation Factor within buffered area
# A VIF for a single explanatory variable is obtained using the r-squared value of the regression of that variable against all other explanatory variables (http://www.r-bloggers.com/collinearity-and-stepwise-vif-selection/)

# x <- sample(1:(dim(present.species.df)[1]), 10000, replace=F) # sample 10k background points for the VIF
# sample.df <- present.species.df[x,]
sample.df <- present.species.df
head(sample.df); dim(sample.df) # 56783    21
# plot(countries); points(sample.df$x, sample.df$y, col='green'); plot(countries, add=T)
sample.matrix <- as.matrix(sample.df)
head(sample.matrix); dim(sample.matrix) # 56783 21

### VIF ###

# keep.dat <- vif_func(in_frame = sample.matrix[,1:19], thresh=10, trace=T) # thresh=5
# keep.dat # "bio02" "bio03" "bio08" "bio13" "bio14" "bio15" "bio18" "bio19"

keep.dat <- colnames(present.species.df[,1:19]) # To use all variables

keep.dat <- c(keep.dat, 'mask')
str(keep.dat)
sample.matrix.keep <- sample.matrix[, (colnames(sample.matrix) %in% keep.dat)]
head(sample.matrix.keep); dim(sample.matrix.keep) # 56783    20 - BACKGROUND SAMPLE
summary(sample.matrix.keep)
sample.df.keep <- data.frame(sample.matrix.keep) 
dim(sample.df.keep) # 56783    20
names(sample.df.keep)

# Species dataframe for keep.dat
head(merbau.gbif.native.unique); dim(merbau.gbif.native.unique); str(merbau.gbif.native.unique) # 103  24
merbau.gbif.native.unique.df <- merbau.gbif.native.unique[, (colnames(merbau.gbif.native.unique) %in% keep.dat)] # retrieve species records
head(merbau.gbif.native.unique.df); dim(merbau.gbif.native.unique.df) # 103  19
merbau.gbif.native.unique.df$mask <- 1 # Add mask column
names(merbau.gbif.native.unique.df)

### Create directory

mainDirMaxent <- "D:/Papers.Projects/Merbau.VTB/R/maxentOutput"

### 6. MAXENT bioclim ####

### CHECK FOLDER NAMES !!!
### Logistic
swd <- rbind(merbau.gbif.native.unique.df, sample.df.keep); dim(swd) # 56886    20; swd dataframe
pa <- c(rep(1, nrow(merbau.gbif.native.unique.df)), rep(0, nrow(sample.df.keep))); length(pa) # presence/absence vector
me <- maxent(swd, pa, args = c("noproduct", "nothreshold", "nohinge", "noextrapolate", "outputformat=logistic", "jackknife", "applyThresholdRule=10 percentile training presence", "projectionlayers=D:/GIS/Worldclim/Present/5arcmin/bio.asc", "redoifexists"), path=file.path(mainDirMaxent))

me
str(me)
me@lambdas
me@results
str(me@results)
me@results['X10.percentile.training.presence.logistic.threshold',] # 0.2309
me@results['Training.AUC',] # 0.8198
plot(me)
response(me, 1:12)
response(me, 13:20)
response(me, var='bio02') # response curves use median values for all other values. Bio02 is (probably) correlated with other vars. Therefor for a single response the maxent model on 1 variable should be run!!!

eval <- evaluate(me, p=merbau.gbif.native.unique.df, a=sample.df.keep)
eval
str(eval)
AUC <- eval@auc
AUC # 0.8198168
threshold(eval)
plot(eval, 'ROC')

# single response bio02
head(swd); str(swd)
swd.bio02 <- as.data.frame(swd[,c('bio02', 'mask')])
head(swd.bio02)
me.bio02 <- maxent(swd.bio02, pa, args = c("noproduct", "nothreshold", "nohinge", "noextrapolate", "outputformat=logistic", "jackknife", "applyThresholdRule=10 percentile training presence", "projectionlayers=D:/GIS/Worldclim/Present/5arcmin/bio.asc", "redoifexists"), path=paste(file.path(mainDirMaxent), '/single.response', sep=""))

me.bio02
response(me.bio02, 'bio02')

# single response bio16
swd.bio16 <- as.data.frame(swd[,c('bio16', 'mask')])
head(swd.bio16)
me.bio16 <- maxent(swd.bio16, pa, args = c("noproduct", "nothreshold", "nohinge", "noextrapolate", "outputformat=logistic", "jackknife", "applyThresholdRule=10 percentile training presence", "projectionlayers=D:/GIS/Worldclim/Present/5arcmin/bio.asc", "redoifexists"), path=paste(file.path(mainDirMaxent), '/single.response', sep=""))

me.bio16
response(me.bio16, 'bio16')

### Null-model ###

maxentResults <- read.csv(paste(file.path(mainDirMaxent), '/', 'maxentResults.csv', sep=""))
# maxentResults
vector <- maxentResults$X.Training.samples
# str(vector)
vector <- as.vector(sort(vector))
vector # 103

head(sample.df.keep); dim(sample.df.keep)
x <- sample.df.keep
head(x)

## Run null-model from source
nm <- nullModel(x, n = vector, rep = 100)
nm # shows the evaluations of the 'rep' null models created
auc <- sapply(nm, function(x){slot(x,'auc')})# get just the auc values of  the null models
auc <- auc[order(auc, decreasing = TRUE)]
hist(auc) #make a histogram
write.csv(auc, paste(file.path(mainDirMaxent), '/', 'nm_auc.csv', sep=""))

auc[5] # 0.6560218
maxentResults$nm <- auc[5] # Add null-model value to maxentResults
write.csv(maxentResults, file=paste(file.path(mainDirMaxent), '/', 'maxentResults.csv', sep=""))

### 7. SDM on ISRIC Soil within climate niche ####

isric.soil.files <- list.files('D:/GIS/ISRIC/wise5by5min_v1b/ascii/', pattern="[.]asc$", full.names=T)
isric.soil.files
isric.soil.df <- as.data.frame(stack(isric.soil.files), xy=T)
head(isric.soil.df); dim(isric.soil.df) # 7240320      21
coordinates(isric.soil.df) <- ~x+y
gridded(isric.soil.df) <- T
head(isric.soil.df)
image(isric.soil.df, 'PHAQ')
ext.isric <- extent(isric.soil.df) # ISRIC has different extent as Worldclim
ext.isric # -180, 179.9986, -56, 83.66611

mask.merbau.present.bioclim <- 'D:/Papers.Projects/Merbau.VTB/R/maxentOutput/species_bio.asc_thresholded.asc'
plot(raster(mask.merbau.present.bioclim))
mask.merbau.present.bioclim.df <- as.data.frame(stack(mask.merbau.present.bioclim), xy=T) # dataframe from merbau.present.bioclim
head(mask.merbau.present.bioclim.df); dim(mask.merbau.present.bioclim.df) # 7776000       3
coordinates(mask.merbau.present.bioclim.df) <- ~x+y
gridded(mask.merbau.present.bioclim.df) <- T
extent(mask.merbau.present.bioclim.df) # -180, 180, -60, 90
mask.merbau.present.bioclim.df <- crop(mask.merbau.present.bioclim.df, ext.isric) # crop to isric extent
image(mask.merbau.present.bioclim.df)
extent(mask.merbau.present.bioclim.df); extent(isric.soil.df)
str(mask.merbau.present.bioclim.df)
mask.merbau.present.bioclim.df@bbox; isric.soil.df@bbox
isric.soil.df@bbox <- mask.merbau.present.bioclim.df@bbox
dim(isric.soil.df); dim(mask.merbau.present.bioclim.df) # 7240320      19; 7240320   1
summary(mask.merbau.present.bioclim.df)
head(mask.merbau.present.bioclim.df@data)

# convert spdf to df to replace 0 for NA
mask.merbau.present.bioclim.df <- cbind(mask.merbau.present.bioclim.df@coords, mask.merbau.present.bioclim.df@data)
head(mask.merbau.present.bioclim.df); str(mask.merbau.present.bioclim.df); dim(mask.merbau.present.bioclim.df) # 7240320       3
mask.merbau.present.bioclim.df[which(mask.merbau.present.bioclim.df$species_bio.asc_thresholded == 0), "species_bio.asc_thresholded"] <- NA # replace 0 to NA
summary(mask.merbau.present.bioclim.df); dim(mask.merbau.present.bioclim.df) # 7240320       3

head(isric.soil.df); str(isric.soil.df); dim(isric.soil.df) # 7240320      19
isric.soil.df@data$mask <- mask.merbau.present.bioclim.df$species_bio.asc_thresholded
head(isric.soil.df); str(isric.soil.df); dim(isric.soil.df) # 7240320      20
isric.soil.df <- cbind(isric.soil.df@coords, isric.soil.df@data) # convert from spdf to df
isric.soil.df <- na.omit(isric.soil.df) # remove missing values --> omit all records outside Merbau bioclim range
dim(isric.soil.df) # 61461    22
summary(isric.soil.df)

x <- isric.soil.df # backup
# isric.soil.df <- x
# rm(x)

# minus indicates no useful substitution (From ISRIC report)

isric.soil.df[which(isric.soil.df$PHAQ < 0), ] <- NA
summary(isric.soil.df)
hist(isric.soil.df$CECc)
isric.soil.df[which(isric.soil.df$CECc < 0), ] <- NA
summary(isric.soil.df)

isric.soil.df <- na.omit(isric.soil.df)
summary(isric.soil.df); dim(isric.soil.df) # 58063    22

coordinates(isric.soil.df) <- ~x+y; gridded(isric.soil.df) <- T # convert df to spdf
image(isric.soil.df, 'PHAQ')
plot(raster(isric.soil.df, 'PHAQ')); plot(countries, add=T)
summary(isric.soil.df); dim(isric.soil.df) # 58063    20
str(isric.soil.df)
isric.soil.df@proj4string <- P4S.latlon
isric.soil.df@data$grid.index <- isric.soil.df@grid.index # Add grid.index value
head(isric.soil.df); dim(isric.soil.df) # 58063    21

### 8. Get ISRIC soil data for Merbau gbif collections within bioclim model prediction ####

merbau.gbif.native.isric <- subset(merbau.gbif, merbau.gbif$decimalLatitude >= native.extent@ymin & merbau.gbif$decimalLatitude <= native.extent@ymax & merbau.gbif$decimalLongitude >= native.extent@xmin & merbau.gbif$decimalLongitude <= native.extent@xmax)
head(merbau.gbif.native.isric); dim(merbau.gbif.native.isric) # 173 3

plot(countries); points(merbau.gbif.native.isric$decimalLongitude, merbau.gbif.native.isric$decimalLatitude, pch=19, col="red")

coordinates(merbau.gbif.native.isric) <- ~decimalLongitude+decimalLatitude
str(merbau.gbif.native.isric)
merbau.gbif.native.isric@proj4string <- P4S.latlon
str(merbau.gbif.native.isric); str(isric.soil.df)
dim(merbau.gbif.native.isric) # 173 1
plot(countries); plot(merbau.gbif.native.isric, col='red', add=T)

merbau.gbif.native.isric.vars <- over(merbau.gbif.native.isric, isric.soil.df) # Get climate variables + grid.index
str(merbau.gbif.native.isric.vars); dim(merbau.gbif.native.isric.vars) # 173 21
merbau.gbif.native.isric.vars <- cbind(merbau.gbif.native.isric@coords, merbau.gbif.native.isric.vars) # add coordinates
head(merbau.gbif.native.isric.vars); dim(merbau.gbif.native.isric.vars) # 173  23
summary(merbau.gbif.native.isric.vars) # 127 NAs, ISRIC is conservative along the coast
merbau.gbif.native.isric.vars <- na.omit(merbau.gbif.native.isric.vars)
head(merbau.gbif.native.isric.vars); dim(merbau.gbif.native.isric.vars) # 46 23
summary(merbau.gbif.native.isric.vars)

duplicates.isric <- duplicated(merbau.gbif.native.isric.vars[,c("grid.index")]) # Duplicates on grid.index
table(duplicates.isric) # 35 F 11 T
merbau.gbif.native.isric.vars <- merbau.gbif.native.isric.vars[!duplicates.isric,] # remove duplicates
str(merbau.gbif.native.isric.vars) # data.frame
head(merbau.gbif.native.isric.vars); dim(merbau.gbif.native.isric.vars) # 35 23; records with ISRIC data

plot(raster(isric.soil.df, 'PHAQ')); points(merbau.gbif.native.isric.vars$decimalLongitude, merbau.gbif.native.isric.vars$decimalLatitude)

boxplot(merbau.gbif.native.isric.vars$PHAQ, main = "pH", ylab="pH")

### 9. Restrict ISRIC background to 500km buffer original dataset ####

head(isric.soil.df); dim(isric.soil.df) # 62369    21; ISRIC data within bioclim model extent
str(isric.soil.df)
image(isric.soil.df, 'PHAQ')
plot(raster("D:/Papers.Projects/Merbau.VTB/R/mask.500km.asc")) # Range within 500km buffer
head(mask.buffer.df); dim(mask.buffer.df); colSums(mask.buffer.df, na.rm=T, dims=1)
mask.buffer.isric.df <- mask.buffer.df
head(mask.buffer.isric.df)
coordinates(mask.buffer.isric.df) <- ~x+y; gridded(mask.buffer.isric.df) <- T
proj4string(mask.buffer.isric.df) <- P4S.latlon # convert to spdf
str(mask.buffer.isric.df)

head(isric.soil.df); dim(isric.soil.df) # 58063    21
str(isric.soil.df)
isric.soil.native.df <- over(isric.soil.df, mask.buffer.isric.df)
head(isric.soil.native.df); dim(isric.soil.native.df)
isric.soil.native.df <- cbind(isric.soil.df@coords, isric.soil.df@data, isric.soil.native.df)
head(isric.soil.native.df); dim(isric.soil.native.df) # 58063    24
isric.soil.native.df <- na.omit(isric.soil.native.df)
head(isric.soil.native.df); dim(isric.soil.native.df) # 18433    24; ISRIC native range background sample
isric.native.background <- isric.soil.native.df[, 3:22]
head(isric.native.background); dim(isric.native.background) # 18433    20
head(merbau.gbif.native.isric.vars); dim(merbau.gbif.native.isric.vars) # 35  23
merbau.gbif.native.isric.swd <- merbau.gbif.native.isric.vars[,3:22]
head(merbau.gbif.native.isric.swd); dim(merbau.gbif.native.isric.swd) # 35 20

### 10. MAXENT ISRIC ####

# mask layer in ISRIC folder, correct header from ISRIC data header
mask.isric <- raster("D:/Papers.Projects/Merbau.VTB/R/maxentOutput/species_bio.asc_thresholded.asc")
plot(mask.isric)
mask.isric <- crop(mask.isric, ext.isric)
mask.isric[mask.isric == 0] <- NA
writeRaster(mask.isric, 'D:/GIS/ISRIC/wise5by5min_v1b/ascii/mask.asc', format='ascii', NAflag = -9999, overwrite = T)

mainDirMaxentISRIC <- "D:/Papers.Projects/Merbau.VTB/R/maxentOutputISRIC/"

### CHECK FOLDER NAMES !!!
### Logistic
swd.isric <- rbind(merbau.gbif.native.isric.swd, isric.native.background); dim(swd.isric); str(swd.isric) # 18468     20 swd dataframe
pa.isric <- c(rep(1, nrow(merbau.gbif.native.isric.swd)), rep(0, nrow(isric.native.background))); length(pa.isric) # presence/absence vector
me.isric <- maxent(swd.isric, pa.isric, args = c("noproduct", "nothreshold", "nohinge", "noextrapolate", "outputformat=logistic", "jackknife", "applyThresholdRule=10 percentile training presence", "projectionlayers=D:/GIS/ISRIC/wise5by5min_v1b/ascii/", "redoifexists"), path=file.path(mainDirMaxentISRIC))

me.isric
str(me.isric)
me.isric@results['X10.percentile.training.presence.logistic.threshold',] #  0.3793
me.isric@results['Training.AUC',] # 0.7036
plot(me.isric)
response(me.isric, 1:9)
response(me.isric, 10:20)

eval.isric <- evaluate(me.isric, p=merbau.gbif.native.isric.swd, a=isric.native.background)
eval.isric
str(eval.isric)
AUC.isric <- eval.isric@auc
AUC.isric # 0.7331161
plot(eval.isric, 'ROC')

### single reponses ISRIC ####

# single response CECc
head(swd.isric); str(swd.isric)
swd.isric.cecc <- as.data.frame(swd.isric[,c('CECc', 'mask')])
head(swd.isric.cecc)
me.swd.isric.cecc <- maxent(swd.isric.cecc, pa.isric, args = c("noproduct", "nothreshold", "nohinge", "noextrapolate", "outputformat=logistic", "jackknife", "applyThresholdRule=10 percentile training presence", "projectionlayers=D:/GIS/ISRIC/wise5by5min_v1b/ascii/", "redoifexists"), path=paste(file.path(mainDirMaxentISRIC), '/single.response', sep=""))

me.swd.isric.cecc
response(me.swd.isric.cecc, 'CECc')

# single response ECEC
head(swd.isric); str(swd.isric)
swd.isric.ecec <- as.data.frame(swd.isric[,c('ECEC', 'mask')])
head(swd.isric.ecec)
me.swd.isric.ecec <- maxent(swd.isric.ecec, pa.isric, args = c("noproduct", "nothreshold", "nohinge", "noextrapolate", "outputformat=logistic", "jackknife", "applyThresholdRule=10 percentile training presence", "projectionlayers=D:/GIS/ISRIC/wise5by5min_v1b/ascii/", "redoifexists"), path=paste(file.path(mainDirMaxentISRIC), '/single.response', sep=""))

me.swd.isric.ecec
response(me.swd.isric.ecec, 'ECEC')

# single response BULK
head(swd.isric); str(swd.isric)
swd.isric.bulk <- as.data.frame(swd.isric[,c('BULK', 'mask')])
head(swd.isric.bulk)
me.swd.isric.bulk <- maxent(swd.isric.bulk, pa.isric, args = c("noproduct", "nothreshold", "nohinge", "noextrapolate", "outputformat=logistic", "jackknife", "applyThresholdRule=10 percentile training presence", "projectionlayers=D:/GIS/ISRIC/wise5by5min_v1b/ascii/", "redoifexists"), path=paste(file.path(mainDirMaxentISRIC), '/single.response', sep=""))

me.swd.isric.bulk
response(me.swd.isric.bulk, 'BULK')

# single response ESP
head(swd.isric); str(swd.isric)
swd.isric.esp <- as.data.frame(swd.isric[,c('ESP', 'mask')])
head(swd.isric.esp)
me.swd.isric.esp <- maxent(swd.isric.esp, pa.isric, args = c("noproduct", "nothreshold", "nohinge", "noextrapolate", "outputformat=logistic", "jackknife", "applyThresholdRule=10 percentile training presence", "projectionlayers=D:/GIS/ISRIC/wise5by5min_v1b/ascii/", "redoifexists"), path=paste(file.path(mainDirMaxentISRIC), '/single.response', sep=""))

me.swd.isric.esp
response(me.swd.isric.esp, 'ESP')

# single response CNrt
head(swd.isric); str(swd.isric)
swd.isric.cnrt <- as.data.frame(swd.isric[,c('CNrt', 'mask')])
head(swd.isric.cnrt)
me.swd.isric.cnrt <- maxent(swd.isric.cnrt, pa.isric, args = c("noproduct", "nothreshold", "nohinge", "noextrapolate", "outputformat=logistic", "jackknife", "applyThresholdRule=10 percentile training presence", "projectionlayers=D:/GIS/ISRIC/wise5by5min_v1b/ascii/", "redoifexists"), path=paste(file.path(mainDirMaxentISRIC), '/single.rcnrtonse', sep=""))

me.swd.isric.cnrt
response(me.swd.isric.cnrt, 'CNrt')

### Null-model ISRIC ####

maxentResultsISRIC <- read.csv(paste(file.path(mainDirMaxentISRIC), '/', 'maxentResults.csv', sep=""))
# maxentResultsISRIC
vector.isric <- maxentResultsISRIC$X.Training.samples
# str(vector)
vector.isric <- as.vector(sort(vector.isric))
vector.isric # 35

head(isric.native.background); dim(isric.native.background)
x <- isric.native.background
head(x)

## Run null-model from source
nm <- nullModel(x, n = vector.isric, rep = 100)
nm # shows the evaluations of the 'rep' null models created
auc <- sapply(nm, function(x){slot(x,'auc')})# get just the auc values of  the null models
auc <- auc[order(auc, decreasing = TRUE)]
hist(auc) #make a histogram
write.csv(auc, paste(file.path(mainDirMaxentISRIC), '/', 'nm_auc.csv', sep=""))

auc[5] #  0.7164486
auc[6] # 0.7122571
auc[7] #  0.71024
auc[8] #  0.7093586
auc[9] #  0.7013586
auc[10] # 0.6981786
maxentResultsISRIC$nm5 <- auc[5] # Add null-model value to maxentResults
maxentResultsISRIC$nm10 <- auc[10] # Add null-model value to maxentResults
write.csv(maxentResultsISRIC, file=paste(file.path(mainDirMaxentISRIC), '/', 'maxentResultsISRIC.csv', sep=""))

### 11. Project to future bioclim ####

# rcp2.6.2050.list <- list.files('D:/GIS/Worldclim/Future/5arcmin/2050/rcp2.6/rcp2.6.2050.worldclim.15gcm.mean'); rcp2.6.2050.list
rcp2.6.2070.list <- list.files('D:/GIS/Worldclim/Future/5arcmin/2070/rcp2.6/rcp2.6.2070.worldclim.15gcm.mean', pattern='[.]asc$'); rcp2.6.2070.list
# rcp8.5.2050.list <- list.files('D:/GIS/Worldclim/Future/5arcmin/2050/rcp8.5/rcp8.5.2050.worldclim.17gcm.mean'); rcp8.5.2050.list
rcp8.5.2070.list <- list.files('D:/GIS/Worldclim/Future/5arcmin/2070/rcp8.5/rcp8.5.2070.worldclim.17gcm.mean', pattern='[.]asc$'); rcp8.5.2070.list

# density.Project rcp2.6 - 2070

mainDirMaxent

# rcp 2.6 - 2070

system(command=paste('java -cp D:/Programs/Maxent.3.3.3k/maxent.jar density.Project ', file.path(mainDirMaxent), '/species.lambdas D:/GIS/Worldclim/Future/5arcmin/2070/rcp2.6/rcp2.6.2070.worldclim.15gcm.mean ', file.path(mainDirMaxent), '/merbau_worldclim_rcp2.6.2070/merbau_worldclim_rcp2.6.2070.noextrapolate noextrapolate', sep=""))

system(command=paste('java -cp D:/Programs/Maxent.3.3.3k/maxent.jar density.Project ', file.path(mainDirMaxent), '/species.lambdas D:/GIS/Worldclim/Future/5arcmin/2070/rcp2.6/rcp2.6.2070.worldclim.15gcm.mean ', file.path(mainDirMaxent), '/merbau_worldclim_rcp2.6.2070/merbau_worldclim_rcp2.6.2070', sep=""))

# rcp8.5 - 2070

system(command=paste('java -cp D:/Programs/Maxent.3.3.3k/maxent.jar density.Project ', file.path(mainDirMaxent), '/species.lambdas D:/GIS/Worldclim/Future/5arcmin/2070/rcp8.5/rcp8.5.2070.worldclim.17gcm.mean ', file.path(mainDirMaxent), '/merbau_worldclim_rcp8.5.2070/merbau_worldclim_rcp8.5.2070.noextrapolate noextrapolate', sep=""))

system(command=paste('java -cp D:/Programs/Maxent.3.3.3k/maxent.jar density.Project ', file.path(mainDirMaxent), '/species.lambdas D:/GIS/Worldclim/Future/5arcmin/2070/rcp8.5/rcp8.5.2070.worldclim.17gcm.mean ', file.path(mainDirMaxent), '/merbau_worldclim_rcp8.5.2070/merbau_worldclim_rcp8.5.2070', sep=""))

### Threshold by 10 percentile values ###

merbau.rcp2.6.2070.files <- list.files('D:/Papers.Projects/Merbau.VTB/R/maxentOutput/merbau_worldclim_rcp2.6.2070', full.name=T, pattern='[.]asc')
merbau.rcp2.6.2070.files
merbau.rcp8.5.2070.files <- list.files('D:/Papers.Projects/Merbau.VTB/R/maxentOutput/merbau_worldclim_rcp8.5.2070', full.name=T, pattern='[.]asc')
merbau.rcp8.5.2070.files

clamping <- grep("clamping", merbau.rcp2.6.2070.files, ignore.case=T, value=T); clamping
merbau.rcp2.6.2070.files <- merbau.rcp2.6.2070.files[! merbau.rcp2.6.2070.files %in% clamping]
merbau.rcp2.6.2070.files
clamping <- grep("clamping", merbau.rcp8.5.2070.files, ignore.case=T, value=T); clamping
merbau.rcp8.5.2070.files <- merbau.rcp8.5.2070.files[! merbau.rcp8.5.2070.files %in% clamping]
merbau.rcp8.5.2070.files

merbau.2070.files <- c(merbau.rcp2.6.2070.files, merbau.rcp8.5.2070.files)
merbau.2070.files

merbau.2070.df <- as.data.frame(stack(merbau.2070.files), xy=T)
head(merbau.2070.df); dim(merbau.2070.df) # 7776000       6

str(me@results)
me@results[which(row.names(me@results) == "X10.percentile.training.presence.logistic.threshold"),] # 0.2309
threshold.10perc <- as.numeric(me@results[which(row.names(me@results)== "X10.percentile.training.presence.logistic.threshold"),])
threshold.10perc # 0.2309

merbau.2070.threshold.df <- merbau.2070.df
head(merbau.2070.threshold.df)
merbau.2070.threshold.df <- na.omit(merbau.2070.threshold.df)
dim(merbau.2070.threshold.df) # 2287025       6

merbau.2070.threshold.df[merbau.2070.threshold.df[,3] >= threshold.10perc, 3] <- 1 # rcp2.6.2070.extrapolate
merbau.2070.threshold.df[merbau.2070.threshold.df[,3] < threshold.10perc, 3] <- 0
merbau.2070.threshold.df[merbau.2070.threshold.df[,4] >= threshold.10perc, 4] <- 1 # rcp2.6.2070.noextrapolate
merbau.2070.threshold.df[merbau.2070.threshold.df[,4] < threshold.10perc, 4] <- 0
merbau.2070.threshold.df[merbau.2070.threshold.df[,5] >= threshold.10perc, 5] <- 1 # rcp8.5.2070.noextrapolate
merbau.2070.threshold.df[merbau.2070.threshold.df[,5] < threshold.10perc, 5] <- 0
merbau.2070.threshold.df[merbau.2070.threshold.df[,6] >= threshold.10perc, 6] <- 1 # rcp8.5.2070.noextrapolate
merbau.2070.threshold.df[merbau.2070.threshold.df[,6] < threshold.10perc, 6] <- 0

head(merbau.2070.threshold.df)
summary(merbau.2070.threshold.df)
colSums(merbau.2070.threshold.df[,3:6])

###############################

coordinates(merbau.2070.threshold.df) <- ~x+y
gridded(merbau.2070.threshold.df) = T
head(merbau.2070.threshold.df); dim(merbau.2070.threshold.df)
# mimage(merbau.2070.threshold.df)

names(merbau.2070.threshold.df)

for(i in 1:4){
  r <- raster(merbau.2070.threshold.df[i])
  # plot(r)
  writeRaster(r,
              filename  = paste('D:/Papers.Projects/Merbau.VTB/R/maxentOutput/', names(merbau.2070.threshold.df)[i], '.thresholded.asc', sep=""),
              format = 'ascii',
              NAflag = -9999,
              overwrite = T)
}

### 11. Project to future ISRIC ####

# Script: Replace mask in ISRIC folder by future thresholded bioclim distributions

ISRIC.list <- list.files('D:/GIS/ISRIC/wise5by5min_v1b/ascii', pattern='[.]asc$'); ISRIC.list

# density.Project ISRIC

mainDirMaxentISRIC

# rcp 2.6 - 2070 - noextrapolate

merbau.bioclim.thesholded  <- list.files('D:/Papers.Projects/Merbau.VTB/R/maxentOutput/', pattern='[.]thresholded[.]asc$', full.name=T)
merbau.bioclim.thesholded # 4 files

r <- crop(raster(merbau.bioclim.thesholded[1]), ext.isric)
r[r == 0] <- NA
plot(r)

writeRaster(r, filename  = "D:/GIS/ISRIC/wise5by5min_v1b/ascii/mask.asc", format = 'ascii', NAflag = -9999, overwrite = T)

### Manually correct ascii header mask !!!

system(command=paste('java -cp D:/Programs/Maxent.3.3.3k/maxent.jar density.Project ', file.path(mainDirMaxentISRIC), '/species.lambdas D:/GIS/ISRIC/wise5by5min_v1b/ascii/ ', file.path(mainDirMaxentISRIC), '/merbau_worldclim_isric_rcp2.6.2070/merbau_worldclim_isric_rcp2.6.2070.noextrapolate noextrapolate', sep=""))

# rcp 2.6 - 2070 - extrapolate

merbau.bioclim.thesholded

r <- crop(raster(merbau.bioclim.thesholded[2]), ext.isric)
r[r == 0] <- NA
plot(r)

writeRaster(r, filename  = "D:/GIS/ISRIC/wise5by5min_v1b/ascii/mask.asc", format = 'ascii', NAflag = -9999, overwrite = T)

### Manually correct ascii header mask !!!

system(command=paste('java -cp D:/Programs/Maxent.3.3.3k/maxent.jar density.Project ', file.path(mainDirMaxentISRIC), '/species.lambdas D:/GIS/ISRIC/wise5by5min_v1b/ascii/ ', file.path(mainDirMaxentISRIC), '/merbau_worldclim_isric_rcp2.6.2070/merbau_worldclim_isric_rcp2.6.2070', sep=""))

# rcp 8.5 - 2070 - noextrapolate

merbau.bioclim.thesholded

r <- crop(raster(merbau.bioclim.thesholded[3]), ext.isric)
r[r == 0] <- NA
plot(r)

writeRaster(r, filename  = "D:/GIS/ISRIC/wise5by5min_v1b/ascii/mask.asc", format = 'ascii', NAflag = -9999, overwrite = T)

### Manually correct ascii header mask !!!

system(command=paste('java -cp D:/Programs/Maxent.3.3.3k/maxent.jar density.Project ', file.path(mainDirMaxentISRIC), '/species.lambdas D:/GIS/ISRIC/wise5by5min_v1b/ascii/ ', file.path(mainDirMaxentISRIC), '/merbau_worldclim_isric_rcp8.5.2070/merbau_worldclim_isric_rcp8.5.2070.noextrapolate noextrapolate', sep=""))

# rcp 8.5 - 2070 - extrapolate

merbau.bioclim.thesholded

r <- crop(raster(merbau.bioclim.thesholded[4]), ext.isric)
r[r == 0] <- NA
plot(r)

writeRaster(r, filename  = "D:/GIS/ISRIC/wise5by5min_v1b/ascii/mask.asc", format = 'ascii', NAflag = -9999, overwrite = T)

### Manually correct ascii header mask !!!

system(command=paste('java -cp D:/Programs/Maxent.3.3.3k/maxent.jar density.Project ', file.path(mainDirMaxentISRIC), '/species.lambdas D:/GIS/ISRIC/wise5by5min_v1b/ascii/ ', file.path(mainDirMaxentISRIC), '/merbau_worldclim_isric_rcp8.5.2070/merbau_worldclim_isric_rcp8.5.2070', sep=""))


### Threshold by 10 percentile values ###

merbau.worldclim.isric.present <- 'D:/Papers.Projects/Merbau.VTB/R/maxentOutputISRIC/species_ascii_thresholded.asc'
plot(raster(merbau.worldclim.isric.present))

merbau.rcp2.6.2070.files <- list.files('D:/Papers.Projects/Merbau.VTB/R/maxentOutputISRIC/merbau_worldclim_isric_rcp2.6.2070', full.name=T, pattern='[.]asc')
merbau.rcp2.6.2070.files
merbau.rcp8.5.2070.files <- list.files('D:/Papers.Projects/Merbau.VTB/R/maxentOutputISRIC/merbau_worldclim_isric_rcp8.5.2070', full.name=T, pattern='[.]asc')
merbau.rcp8.5.2070.files

clamping <- grep("clamping", merbau.rcp2.6.2070.files, ignore.case=T, value=T); clamping
merbau.rcp2.6.2070.files <- merbau.rcp2.6.2070.files[! merbau.rcp2.6.2070.files %in% clamping]
merbau.rcp2.6.2070.files
clamping <- grep("clamping", merbau.rcp8.5.2070.files, ignore.case=T, value=T); clamping
merbau.rcp8.5.2070.files <- merbau.rcp8.5.2070.files[! merbau.rcp8.5.2070.files %in% clamping]
merbau.rcp8.5.2070.files

merbau.2070.files <- c(merbau.worldclim.isric.present, merbau.rcp2.6.2070.files, merbau.rcp8.5.2070.files)
merbau.2070.files

merbau.2070.df <- as.data.frame(stack(merbau.2070.files), xy=T)
head(merbau.2070.df); dim(merbau.2070.df) # 7240320       7
summary(merbau.2070.df)

str(me.isric@results)
me.isric@results[which(row.names(me.isric@results) == "X10.percentile.training.presence.logistic.threshold"),] # 0.3793
threshold.10perc <- as.numeric(me.isric@results[which(row.names(me.isric@results)== "X10.percentile.training.presence.logistic.threshold"),])
threshold.10perc # 0.3793

merbau.2070.threshold.df <- merbau.2070.df
colSums(merbau.2070.threshold.df, na.rm=T)
summary(merbau.2070.threshold.df)
head(merbau.2070.threshold.df)
dim(merbau.2070.threshold.df) # 7240320       7
merbau.2070.threshold.df[is.na(merbau.2070.threshold.df)] <- 0 # replace cannot deal with NA

merbau.2070.threshold.df[merbau.2070.threshold.df[,4] >= threshold.10perc, 4] <- 1 # rcp2.6.2070.extrapolate
merbau.2070.threshold.df[merbau.2070.threshold.df[,4] < threshold.10perc, 4] <- 0
merbau.2070.threshold.df[merbau.2070.threshold.df[,5] >= threshold.10perc, 5] <- 1 # rcp2.6.2070.noextrapolate
merbau.2070.threshold.df[merbau.2070.threshold.df[,5] < threshold.10perc, 5] <- 0
merbau.2070.threshold.df[merbau.2070.threshold.df[,6] >= threshold.10perc, 6] <- 1 # rcp8.5.2070.noextrapolate
merbau.2070.threshold.df[merbau.2070.threshold.df[,6] < threshold.10perc, 6] <- 0
merbau.2070.threshold.df[merbau.2070.threshold.df[,7] >= threshold.10perc, 7] <- 1 # rcp8.5.2070.noextrapolate
merbau.2070.threshold.df[merbau.2070.threshold.df[,7] < threshold.10perc, 7] <- 0

### Gain, Loss, Stable ####
### Future minus present

head(merbau.2070.threshold.df)

merbau.2070.threshold.df$gls.2.6.2070 <- merbau.2070.threshold.df$merbau_worldclim_isric_rcp2.6.2070 - merbau.2070.threshold.df$species_ascii_thresholded

merbau.2070.threshold.df$gls.2.6.2070[merbau.2070.threshold.df$gls.2.6.2070 == 0] <- NA # set both 1-1 and 0-0 to NA
merbau.2070.threshold.df$gls.2.6.2070[merbau.2070.threshold.df$species_ascii_thresholded == 1 & merbau.2070.threshold.df$merbau_worldclim_isric_rcp2.6.2070 == 1] <- 0 # set present = 1 and future = 1 to stable = 0

merbau.2070.threshold.df$gls.2.6.2070.noextr <- merbau.2070.threshold.df$merbau_worldclim_isric_rcp2.6.2070.noextrapolate - merbau.2070.threshold.df$species_ascii_thresholded

merbau.2070.threshold.df$gls.2.6.2070.noextr[merbau.2070.threshold.df$gls.2.6.2070.noextr == 0] <- NA # set both 1-1 and 0-0 to NA
merbau.2070.threshold.df$gls.2.6.2070.noextr[merbau.2070.threshold.df$species_ascii_thresholded == 1 & merbau.2070.threshold.df$merbau_worldclim_isric_rcp2.6.2070.noextrapolate == 1] <- 0 # set present = 1 and future = 1 to stable = 0

merbau.2070.threshold.df$gls.8.5.2070 <- merbau.2070.threshold.df$merbau_worldclim_isric_rcp8.5.2070 - merbau.2070.threshold.df$species_ascii_thresholded

merbau.2070.threshold.df$gls.8.5.2070[merbau.2070.threshold.df$gls.8.5.2070 == 0] <- NA # set both 1-1 and 0-0 to NA
merbau.2070.threshold.df$gls.8.5.2070[merbau.2070.threshold.df$species_ascii_thresholded == 1 & merbau.2070.threshold.df$merbau_worldclim_isric_rcp8.5.2070 == 1] <- 0 # set present = 1 and future = 1 to stable = 0

merbau.2070.threshold.df$gls.8.5.2070.noextr <- merbau.2070.threshold.df$merbau_worldclim_isric_rcp8.5.2070.noextrapolate - merbau.2070.threshold.df$species_ascii_thresholded

merbau.2070.threshold.df$gls.8.5.2070.noextr[merbau.2070.threshold.df$gls.8.5.2070.noextr == 0] <- NA # set both 1-1 and 0-0 to NA
merbau.2070.threshold.df$gls.8.5.2070.noextr[merbau.2070.threshold.df$species_ascii_thresholded == 1 & merbau.2070.threshold.df$merbau_worldclim_isric_rcp8.5.2070.noextrapolate == 1] <- 0 # set present = 1 and future = 1 to stable = 0

summary(merbau.2070.threshold.df)
dim(merbau.2070.threshold.df)

# merbau.2070.threshold.df[merbau.2070.threshold.df[,3] == 0, 3] <- NA
# merbau.2070.threshold.df[merbau.2070.threshold.df[,4] == 0, 4] <- NA
# merbau.2070.threshold.df[merbau.2070.threshold.df[,5] == 0, 5] <- NA
# merbau.2070.threshold.df[merbau.2070.threshold.df[,6] == 0, 6] <- NA
# merbau.2070.threshold.df[merbau.2070.threshold.df[,7] == 0, 7] <- NA

head(merbau.2070.threshold.df)
summary(merbau.2070.threshold.df)
colSums(merbau.2070.threshold.df[,3:11], na.rm=T)

coordinates(merbau.2070.threshold.df) <- ~x+y
gridded(merbau.2070.threshold.df) = T
head(merbau.2070.threshold.df); dim(merbau.2070.threshold.df)
mimage(merbau.2070.threshold.df, 6:9)
image(merbau.2070.threshold.df, 6)

names(merbau.2070.threshold.df)

for(i in 6:9){
  r <- raster(merbau.2070.threshold.df[i])
  # plot(r)
  writeRaster(r,
              filename  = paste('D:/Papers.Projects/Merbau.VTB/R/Future/', names(merbau.2070.threshold.df)[i], '.thresholded.asc', sep=""),
              format = 'ascii',
              NAflag = -9999,
              overwrite = T)
}


### GlobCover 2009 ####

glob.cover.2009 <- raster('D:/GIS/Globcover2009_V2.3_Global/GLOBCOVER_L4_200901_200912_V2.3.tif')
str(glob.cover.2009)
dim(glob.cover.2009)
legend <- read.xlsx('D:/GIS/Globcover2009_V2.3_Global/Globcover2009_Legend.xls', 1)
head(legend)
glob.cover.2009@crs <- P4S.latlon
glob.cover.2009.africa <- crop(glob.cover.2009, ext.africa)

writeRaster(glob.cover.2009.africa, 
            filename  = "D:/Papers.Projects/Andel.Tinde.van/medicinal.plants/predictors/mask/glob.cover.2009.africa.asc",
            format    = 'ascii',
            NAflag    = -9999,
            overwrite = T)

rm(glob.cover.2009)
plot(glob.cover.2009.africa)
summary(glob.cover.2009.africa)
glob.cover.2009.africa[glob.cover.2009.africa >= 210] <- NA # remove water & ice etc. see legend.xls
glob.cover.2009.africa[glob.cover.2009.africa >= 190] <- 0 # set urban & bare to 0
glob.cover.2009.africa.natural.vegetation <- glob.cover.2009.africa
glob.cover.2009.africa.natural.vegetation[glob.cover.2009.africa.natural.vegetation <= 30] <- 0 # set cropland to 0
glob.cover.2009.africa.natural.vegetation[glob.cover.2009.africa.natural.vegetation > 0] <- 1 # set natural forest cover to 1

writeRaster(glob.cover.2009.africa.natural.vegetation, 
            filename  = "D:/Papers.Projects/Andel.Tinde.van/medicinal.plants/predictors/mask/gc2009nvhighres.asc",
            format    = 'ascii',
            NAflag    = -9999,
            overwrite = T)

str(glob.cover.2009.africa.natural.vegetation)
plot(glob.cover.2009.africa.natural.vegetation)
glob.cover.2009.africa.natural.vegetation@ncols #18720
bio01.TA@ncols # 624
glob.cover.2009.africa.natural.vegetation@ncols/bio01.TA@ncols # 30
gc.2009.africa.nv.mean.5arcmin <- aggregate(glob.cover.2009.africa.natural.vegetation, fact = 30, fun=mean) # % natural  vegetation at 5 arcmin
plot(gc.2009.africa.nv.mean.5arcmin)
str(gc.2009.africa.nv.mean.5arcmin)
gc.2009.africa.nv.05.5arcmin <- gc.2009.africa.nv.mean.5arcmin
gc.2009.africa.nv.05.5arcmin[gc.2009.africa.nv.05.5arcmin > 0.5] <- 1 # above 50% present
gc.2009.africa.nv.05.5arcmin[gc.2009.africa.nv.05.5arcmin <= 0.5] <- 0 # below 50% absent
gc.2009.africa.nv.05.5arcmin[gc.2009.africa.nv.05.5arcmin <= 0.5] <- NA # below 50% absent
plot(gc.2009.africa.nv.05.5arcmin, col=rainbow(2))
summary(gc.2009.africa.nv.05.5arcmin)

writeRaster(gc.2009.africa.nv.05.5arcmin, 
            filename  = "D:/Papers.Projects/Andel.Tinde.van/medicinal.plants/predictors/mask/gc.2009.africa.nv.05.5arcmin.asc",
            format    = 'ascii',
            NAflag    = -9999,
            overwrite = T)

gc.2009.africa.nv.min.5arcmin <- aggregate(glob.cover.2009.africa.natural.vegetation, fact = 30, fun=min) # % natural  vegetation at 5 arcmin
plot(gc.2009.africa.nv.min.5arcmin)

writeRaster(gc.2009.africa.nv.mean.5arcmin, 
            filename  = "D:/Papers.Projects/Andel.Tinde.van/medicinal.plants/predictors/mask/gc2009nvmean.asc",
            format    = 'ascii',
            NAflag    = -9999,
            overwrite = T)

writeRaster(gc.2009.africa.nv.min.5arcmin, 
            filename  = "D:/Papers.Projects/Andel.Tinde.van/medicinal.plants/predictors/mask/gc2009nvmin.asc",
            format    = 'ascii',
            NAflag    = -9999,
            overwrite = T)

### Range measures ####

sp.dist <- list.files('D:/Papers.Projects/Andel.Tinde.van/medicinal.plants/maxent.output/', pattern="_present_thresholded[.]asc$", full.names=T)
gc.nv <- "D:/Papers.Projects/Andel.Tinde.van/medicinal.plants/predictors/mask/gc.2009.africa.nv.05.5arcmin.asc"
sp.dist <- c(sp.dist, gc.nv)
names.sp.dist <- gsub('_present_thresholded.asc', '', basename(sp.dist))
names.sp.dist <- gsub('gc.2009.africa.nv.05.5arcmin.asc', 'gc.nv', names.sp.dist)
sp.dist.asc2df <- asc2dataframe(sp.dist, varnames = names.sp.dist)
head(sp.dist.asc2df)
dim(sp.dist.asc2df)
colnames(sp.dist.asc2df)
table(sp.dist.asc2df[1:5,3:17])
x <- data.frame(colSums(sp.dist.asc2df[,3:16]))
y <- data.frame(colSums(sp.dist.asc2df[sp.dist.asc2df$gc.nv==1,3:16]))
x$nv <- y[,1]
x$nv.perc <- round((x[,2]/x[,1])*100,0)
write.csv(x, file='../ms/range.sizes.csv')

### From dismo ##########
ecocrop('potato', 5:16, 15:26, runif(12)*100)
getCrop('Acacia brachystachya Benth.')
crop <- getCrop('Hot pepper')
ecocrop(crop, 5:16, 15:26, rainfed=FALSE)
