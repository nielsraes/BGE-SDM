rm(list = ls(all=T))
setwd("D:/Github/BGE-SDM")
#save(list=ls(all=TRUE), file="D:/Github/BGE-SDM/SDM.RData") # save RDATA for later use
#load("D:/Github/BGE-SDM/SDM.RData")

#### Install required libraries ####
library(rgbif)
library(raster) 
library(sp)
library(mapr)
library(dismo)
library(rgeos)


#### Define extent and import study area shapefile ####
#Import shapefile
countries <- rgdal::readOGR("D:/Github/BGE-SDM/GISDATA/Study Area SHP")
#dev.off()
proj4string(countries) <- "+proj=longlat +datum=WGS84"
plot(countries)
str(countries); countries@bbox

#Define extent
extent <- extent(-43, 109, 25, 81) ##NEED TO BE REARRANGE FOR EAST SIDE

#### 1. GBIF DATA ####
#get the key of species by its name
key = name_backbone(name="Oeneis jutta")$speciesKey

#search data according to key. **occ_data vs occ_Search
dat <- occ_data(taxonKey = key,limit=1000) #limit??
names(dat)
names(dat$data)

#assign data to a variable
dat.data <- dat$data
dat.data

#plot occurences
map_plot(dat.data, lon = "decimalLongitude", lat = "decimalLatitude", size = 1, pch = 3)

#plot occurrences data and the study area
plot(countries); points(dat.data$decimalLongitude, dat.data$decimalLatitude, pch=19, col="red")

plot(countries)

#Prepare the data
fe.gbif <- dat.data[, c('species', 'decimalLongitude', 'decimalLatitude')]
head(fe.gbif); dim(fe.gbif)
duplicates <- duplicated(fe.gbif)
fe.gbif <- fe.gbif[!duplicates,]
head(fe.gbif); dim(fe.gbif)
write.csv(fe.gbif, 'D:/Github/BGE-SDM/Output/Oeneis_jutta.csv', row.names=F)
fe.gbif.maxent <- fe.gbif[, c("species", "decimalLongitude", "decimalLatitude")]
names(fe.gbif.maxent) <- c("species", "lon", "lat")
head(fe.gbif.maxent); dim(fe.gbif.maxent)
plot(countries); points(fe.gbif.maxent$lon, fe.gbif.maxent$lat, pch=19, col='red')
write.csv(fe.gbif.maxent, 'D:/Github/BGE-SDM/Output/fe.gbif.maxent.csv', row.names=F) # write species points file

#assign the fe.gbif as SpatialPointsDataFrame
coordinates(fe.gbif) <- ~decimalLongitude+decimalLatitude
str(fe.gbif)
#assign CRS
P4S.latlon <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")


fe.gbif@proj4string <- P4S.latlon
plot(countries); plot(fe.gbif, col='red', add=T, pch = 16)


### 2. Compile climate data ####
# WORLDCLIM DATA #

# Download precipitation/Bio/+++ data from WorldClim
global.clim <- getData("worldclim", var="bio", res=5, download=T, path="RDATA")

files.present.bil <- list.files('D:/Github/BGE-SDM/RDATA/wc5/', pattern="[.]bil$", full.names=T) # alternatives for pattern (c|C)(e|E)(l|L)$
files.present.bil

# Loop for cropping with extent and write as ascii
for(i in files.present.bil)  {
  raster <- raster(i)
  raster <- crop(raster, extent)
  writeRaster(raster,
              filename  = paste("D:/Github/BGE-SDM/RDATA/clipped/", "fe_buffer_", basename(i), sep = ""),
              format    = 'ascii',
              NAflag    = -9999,
              overwrite = T)
}

# read worldclim
files.present <- list.files('D:/Github/BGE-SDM/RDATA/clipped/', pattern="[.]asc$", full.names=T) # alternatives for pattern (c|C)(e|E)(l|L)$
files.present
present.stack <- stack(files.present)
head(present.stack)
present.df <- as.data.frame(present.stack, xy=T)
coordinates(present.df) <- ~x+y
gridded(present.df) <- T
present.df@proj4string <- P4S.latlon
present.df$grid.index <- present.df@grid.index # Add grid.index value
head(present.df)
image(present.df, 'fe_buffer_bio01') #Annual Mean Temperature
str(present.df)
class(present.df)
present.df@coords
present.df@proj4string <- P4S.latlon

### 3. Get abiotic bioclim data ####
str(fe.gbif); str(present.df)
fe.gbif.abiotic <- over(fe.gbif, present.df) # Get climate variables + grid.index
str(fe.gbif.abiotic) 
fe.gbif.abiotic
class(present.df)
head(fe.gbif.abiotic)
dim(fe.gbif.abiotic); dim(fe.gbif)
names(fe.gbif)
head(fe.gbif); str(fe.gbif)
fe.gbif <- cbind(fe.gbif, fe.gbif.abiotic) # Link species col and climate data
names(fe.gbif)
head(fe.gbif)
#error?? = duplicated() applies only to vectors sol: as.data.frame()
duplicates <- duplicated(as.data.frame(fe.gbif)[,c("species", "grid.index")]) # Duplicates on grid.index
names(fe.gbif)
table(duplicates)
summary(fe.gbif)
fe.gbif <- na.omit(fe.gbif)
head(fe.gbif) ; dim(fe.gbif)

plot(raster(present.df, 'fe_buffer_bio01')); points(fe.gbif$decimalLongitude, fe.gbif$decimalLatitude) #Mean annual Temperature

boxplot(fe.gbif$fe_buffer_bio01, main = "Mean annual Temperature", ylab="Temperature x 10 (in °C)")
boxplot(fe.gbif$fe_buffer_bio12, main = "Annual precipitation", ylab="Precipitation (in mm)")


### 4. Run Maxent model with bioclim data in 500 km buffered area around presences to balance prevalence ####

# climate has priority over soil - hierarchical model, 1 climate, 2 soil
# Boucher-Lalonde, V., A. Morin and D. J. Currie (2012). "How are tree species distributed in climatic space? A simple and general pattern." Global Ecology and Biogeography 21(12): 1157-1166.

### Create empty mask layer

mask <- raster(files.present[1])
#dev.off()
plot(mask)
mask <- !is.na(mask) # all values to 1
mask[mask == 0] <- NA # zero values to NA
plot(mask)
summary(mask)
writeRaster(mask, filename  = "D:/Github/BGE-SDM/Output/mask.asc", format = 'ascii', NAflag = -9999, overwrite = T)
mask <- raster('D:/Github/BGE-SDM/Output/mask.asc')
plot(mask, col='red')

# add mask to present.df
present.df@data$mask <- as.data.frame(stack('D:/Github/BGE-SDM/Output/mask.asc')) # Add mask layer to spdf
head(present.df); dim(present.df)
head(fe.gbif)

# Convert to coll.locs spatial Points Data Frame and create 500 km buffer
coordinates(fe.gbif.maxent) <- ~lon+lat
proj4string(fe.gbif.maxent) <- P4S.latlon
plot(countries); points(fe.gbif.maxent, pch=19, cex=0.5, col='red')

x <- circles(fe.gbif.maxent, d=500000, lonlat=TRUE) # 500 km
pol <- gUnaryUnion(x@polygons) # dissolve polygons
extent(pol)
plot(pol, col='blue', add=T); points(fe.gbif.maxent, pch=19, cex=0.5, col='red')


# extract cell numbers for the circles
v <- extract(mask, x@polygons, cellnumbers=T)
str(v)
# use rbind to combine the elements in list v
v <- do.call(rbind, v)
head(v); dim(v) 

# remove ocean cells
v <- unique(na.omit(v))
head(v); dim(v)

# to display the results
m <- mask
m[] <- NA # empty mask
m[as.vector(v[,1])] <- 1
plot(m, col='purple')
extent(m)
str(m); summary(m)
plot(m, ext=extent(x@polygons)+1, col='blue')
plot(x@polygons, add=T)
points(fe.gbif.maxent, pch=19, cex=0.5, col='red')
plot(countries, add=T)
str(m) # rasterlayer

# Write mask for buffered areas
writeRaster(m, filename  = "D:/Github/BGE-SDM/Output/mask.500km.asc", format = 'ascii', NAflag = -9999, overwrite = T)

### Add mask.buffer to present.df
head(present.df); dim(present.df)
# present.df <- subset(present.df, select = -c(mask)) # remove mask column

present.species.df <- present.df # copy present.df
mask.buffer.df <- as.data.frame(stack('D:/Github/BGE-SDM/Output/mask.500km.asc'), xy=T) # read mask layer
head(mask.buffer.df); dim(mask.buffer.df); colSums(mask.buffer.df, na.rm=T, dims=1)
head(present.species.df); str(present.species.df)

#error undefined columns selected
present.species.df$mask <- mask.buffer.df[,'mask.500km'] # replace mask with mask.500km

head(present.species.df); dim(present.species.df)
str(present.species.df) # SpatialPixelsDataFrame
present.species.df <- na.omit(present.species.df@data)
head(present.species.df); dim(present.species.df)


### 5. Select uncorrelated variables using VIF ####
### Variance Inflation Factor within buffered area
# A VIF for a single explanatory variable is obtained using the r-squared value of the regression of that variable against all other explanatory variables (http://www.r-bloggers.com/collinearity-and-stepwise-vif-selection/)

x <- sample(1:(dim(present.species.df)[1]), 10000, replace=F) # sample 10k background points for the VIF
sample.df <- present.species.df[x,]
sample.df <- present.species.df
head(sample.df); dim(sample.df) 
#plot(countries); points(sample.df$x, sample.df$y, col='green'); plot(countries, add=T)
sample.matrix <- as.matrix(sample.df)
head(sample.matrix); dim(sample.matrix)

### VIF ###

keep.dat <- colnames(present.species.df[,1:19]) # To use all variables

keep.dat <- c(keep.dat, 'mask')
str(keep.dat)
sample.matrix.keep <- sample.matrix[, (colnames(sample.matrix) %in% keep.dat)]
head(sample.matrix.keep); dim(sample.matrix.keep)
summary(sample.matrix.keep)

sample.df.keep <- data.frame(sample.matrix.keep) 
dim(sample.df.keep)
names(sample.df.keep)

# Species dataframe for keep.dat
head(fe.gbif); dim(fe.gbif); str(fe.gbif) 
fe.gbif.df <- fe.gbif[, (colnames(fe.gbif) %in% keep.dat)] # retrieve species records
head(fe.gbif.df); dim(fe.gbif.df) 
fe.gbif.df$mask <- 1 # Add mask column
names(fe.gbif.df)

### Create directory

mainDirMaxent <- "D:/Github/BGE-SDM/Output/"

### 6. MAXENT bioclim ####

### CHECK FOLDER NAMES !!!
### Logistic
#error?
swd <- rbind(as.data.frame(fe.gbif.df), sample.df.keep); dim(swd) 


pa <- c(rep(1, nrow(fe.gbif.df)), rep(0, nrow(sample.df.keep))); length(pa) # presence/absence vector
me <- maxent(swd, pa, args = c("noproduct", "nothreshold", "nohinge", "noextrapolate", "outputformat=logistic", "jackknife", "applyThresholdRule=10 percentile training presence", "projectionlayers=D:/Github/BGE-SDM/RDATA/wc5/bio.asc", "redoifexists"), path=file.path(mainDirMaxent))

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

eval <- evaluate(me, p=fe.gbif.df, a=sample.df.keep)
eval
str(eval)
AUC <- eval@auc
AUC # 0.8198168
threshold(eval)
plot(eval, 'ROC')

# single response bio02 #BIO2 = Mean Diurnal Range (Mean of monthly (max temp - min temp))
head(swd); str(swd)
swd.bio02 <- as.data.frame(swd[,c('bio02', 'mask')])
head(swd.bio02)
me.bio02 <- maxent(swd.bio02, pa, args = c("noproduct", "nothreshold", "nohinge", "noextrapolate", "outputformat=logistic", "jackknife", "applyThresholdRule=10 percentile training presence", "projectionlayers=D:/GIS/Worldclim/Present/5arcmin/bio.asc", "redoifexists"), path=paste(file.path(mainDirMaxent), '/single.response', sep=""))

me.bio02
response(me.bio02, 'bio02')

# single response bio16 #BIO16 = Precipitation of Wettest Quarter

### Null-model ###
## Run null-model from source

### 7. SDM on ISRIC Soil within climate niche ####

### 8. Get ISRIC soil data for fe gbif collections within bioclim model prediction ####
##??



### 9. Restrict ISRIC background to 500km buffer original dataset ####

### 10. MAXENT ISRIC ####

### 11. Project to future bioclim ####

### 11. Project to future ISRIC ####
