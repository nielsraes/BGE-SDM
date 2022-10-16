rm(list = ls(all=T))

#### Install required libraries ####

library(rgbif)
library(raster) 
library(sp)

#### 1. GBIF DATA ####

#get the key of species by its name
key = name_backbone(name="Boloria alaskensis")$speciesKey

#search data according to key. occ_data vs occ_Search
dat <- occ_data(taxonKey = key,limit=10000)

#assign data to a variable
dat.data <- dat$data

# get the columns that matter for mapping and cleaning the occurrence data:
dat.filter <- dat.data[,c("key", "scientificName", "decimalLatitude", "decimalLongitude", "basisOfRecord", "speciesKey", "species", "year","coordinateUncertaintyInMeters")] #alternative columns
dat.filter


### 2. WORLDCLIM DATA ####

# Download precipitation/Bio/+++ data from WorldClim
global.precip <- getData("worldclim", var="prec", res=5, download=T, path="RDATA")
global.clim <- getData("worldclim", var="bio", res=5, download=T, path="RDATA")


# read worldclim ####
files.present <- list.files('D:/Github/BGE-SDM/RDATA/wc5', pattern="[.]bil$", full.names=T) # alternatives for pattern (c|C)(e|E)(l|L)$
# files <- list.files('Z:/World/Climate/Worldclim/05arcmin/Present/bio/', pattern="[.]bil$", full.names=T) # alternatives for pattern (c|C)(e|E)(l|L)$
files.present <- files.present[2:20] # remove alt
files.present

### Extent of the study area
extent <- extent(-43, 109, 25, 81)

# countries boundaries
countries <- rgdal::readOGR("D:/Github/BGE-SDM/GISDATA/Study Area SHP")
proj4string(countries) <- "+proj=longlat +datum=WGS84"
plot(countries)
str(countries)

### Loop for cropping with extent and write as ascii ####
for(i in files.present)  {
  raster <- raster(i)
  raster <- crop(raster, extent)
  writeRaster(raster,
              filename  = paste("D:/Github/BGE-SDM/RDATA/clipped/", "sunda.", basename(i), sep = ""),
              format    = 'ascii',
              NAflag    = -9999,
              overwrite = T)
}
