rm(list = ls(all=T))
setwd("D:/Github/BGE-SDM")
#save(list=ls(all=TRUE), file="D:/Github/BGE-SDM/SDM.RData") # save RDATA for later use
load("D:/Github/BGE-SDM/SDM.RData")
#### Install required libraries ####
library(rgbif)
library(raster) 
library(sp)
library(mapr)

#### 1. GBIF DATA ####
#get the key of species by its name
key = name_backbone(name="Oeneis jutta")$speciesKey

#search data according to key. **occ_data vs occ_Search
dat <- occ_data(taxonKey = key,limit=1000)
head(dat)

names(dat)
names(dat$data)

#plot occurences
map_plot(dat, lon = "decimalLongitude", lat = "decimalLatitude", size = 1, pch = 3)

#assign data to a variable
dat.data <- dat$data
dat.data

# get the columns that matter for mapping and cleaning the occurrence data:
dat.filter <- dat.data[,c("key", "scientificName", "decimalLatitude", "decimalLongitude", "basisOfRecord", "speciesKey", "species", "year","coordinateUncertaintyInMeters")] #alternative columns
dat.filter

### Extent of the study area
extent <- extent(-43, 109, 25, 81) ##NEED TO BE REARRANGE


# countries boundaries
countries <- rgdal::readOGR("D:/Github/BGE-SDM/GISDATA/Study Area SHP")
dev.off()
proj4string(countries) <- "+proj=longlat +datum=WGS84"
countries <- crop(countries, extent)
plot(countries)
str(countries)
countries@bbox

#plot occurences data on the study area
plot(countries); points(dat.data$decimalLongitude, dat.data$decimalLatitude, pch=19, col="red")

plot(countries)

#prepare data for maxent
fe.gbif <- dat.data[, c('species', 'decimalLongitude', 'decimalLatitude')]
head(fe.gbif); dim(fe.gbif)
duplicates <- duplicated(fe.gbif)
fe.gbif <- fe.gbif[!duplicates,]
head(fe.gbif); dim(fe.gbif)
write.csv(fe.gbif, 'D:/Github/BGE-SDM/Output/Oeneis_jutta.csv', row.names=F)
coordinates(fe.gbif) <- ~decimalLongitude+decimalLatitude
str(fe.gbif)
P4S.latlon <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
fe.gbif@proj4string <- P4S.latlon
plot(countries); plot(fe.gbif, col='purple', add=T)


### 2. Compile climate data ####
### WORLDCLIM DATA ####

# Download precipitation/Bio/+++ data from WorldClim

global.clim <- getData("worldclim", var="bio", res=5, download=T, path="RDATA")

files.present.bil <- list.files('D:/Github/BGE-SDM/RDATA/wc5/', pattern="[.]bil$", full.names=T) # alternatives for pattern (c|C)(e|E)(l|L)$

files.present.bil

### Loop for cropping with extent and write as ascii ####
for(i in files.present.bil)  {
  raster <- raster(i)
  raster <- crop(raster, extent)
  writeRaster(raster,
              filename  = paste("D:/Github/BGE-SDM/RDATA/clipped/", "fe_buffer_", basename(i), sep = ""),
              format    = 'ascii',
              NAflag    = -9999,
              overwrite = T)
}

# read worldclim ####
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
image(present.df, 'fe_buffer_bio01')

str(present.df)
class(present.df)
present.df@coords
present.df@proj4string <- P4S.latlon


### 3. Get abiotic bioclim data ####

str(fe.gbif); str(present.df)
fe.gbif.abiotic <- over(fe.gbif, present.df) # Get climate variables + grid.index
str(fe.gbif.abiotic) # 173  20
fe.gbif.abiotic
class(present.df)
head(fe.gbif.abiotic)
dim(fe.gbif.abiotic); dim(fe.gbif) # 173  20  1
head(fe.gbif); str(fe.gbif)
fe.gbif <- cbind(fe.gbif, fe.gbif.abiotic) # Link species col and climate data
head(fe.gbif)

#error?? = Error in duplicated.default(fe.gbif[, c("species", "grid.index")]) :duplicated() applies only to vectors
duplicates <- duplicated(fe.gbif[,c("species", "grid.index")]) # Duplicates on grid.index


table(duplicates) # 114 F 59 T
plot(raster(present.df, 'fe_buffer_bio01')); points(fe.gbif$decimalLongitude, fe.gbif$decimalLatitude)

boxplot(fe.gbif$fe_buffer_bio01, main = "Mean annual Temperature", ylab="Temperature x 10 (in °C)")
boxplot(fe.gbif$fe_buffer_bio12, main = "Annual precipitation", ylab="Precipitation (in mm)")


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
writeRaster(mask, filename  = "D:/Github/BGE-SDM/Output/mask.asc", format = 'ascii', NAflag = -9999, overwrite = T)
mask <- raster('D:/Github/BGE-SDM/Output/mask.asc')
plot(mask, col='red')

# add mask to present.df
present.df@data$mask <- as.data.frame(stack('D:/Github/BGE-SDM/Output/mask.asc')) # Add mask layer to spdf
head(present.df); dim(present.df) # 7776000      21

head(fe.gbif)


fe.gbif.maxent <- fe.gbif$data[, c("species", "decimalLongitude", "decimalLatitude")]
#ERROR?
names(fe.gbif.maxent) <- c("species", "lon", "lat")
head(fe.gbif.maxent); dim(fe.gbif.maxent) # 103   3
plot(countries); points(fe.gbif.maxent$lon, fe.gbif.maxent$lat, pch=19, col='red')

write.csv(fe.gbif.maxent, 'D:/Github/BGE-SDM/Output/fe.gbif.maxent.csv', row.names=F) # write species points file


### 5. Select uncorrelated variables using VIF ####
### 6. MAXENT bioclim ####
### 7. SDM on ISRIC Soil within climate niche ####

### 8. Get ISRIC soil data for fe gbif collections within bioclim model prediction ####
##??
fe.gbif.isric <- subset(fe.gbif, fe.gbif$decimalLatitude >= native.extent@ymin & fe.gbif$decimalLatitude <= native.extent@ymax & fe.gbif$decimalLongitude >= native.extent@xmin & fe.gbif$decimalLongitude <= native.extent@xmax)
head(fe.gbif.isric); dim(fe.gbif.isric) # 173 3
plot(countries); points(fe.gbif.isric$decimalLongitude, fe.gbif.isric$decimalLatitude, pch=19, col="red")

coordinates(fe.gbif.isric) <- ~decimalLongitude+decimalLatitude
str(fe.gbif.isric)
fe.gbif.isric@proj4string <- P4S.latlon
str(fe.gbif.isric); str(isric.soil.df)
dim(fe.gbif.isric) # 173 1
plot(countries); plot(fe.gbif.isric, col='red', add=T)




### 9. Restrict ISRIC background to 500km buffer original dataset ####
### 10. MAXENT ISRIC ####
### 11. Project to future bioclim ####
### 11. Project to future ISRIC ####
