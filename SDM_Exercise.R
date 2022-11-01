#### Install required libraries ####

library(rgbif)
library(raster) 
library(sp)

#### 1. GBIF DATA ####

#Importing occurrence data

#get the key of species by its name
key = name_backbone(name="Oeneis jutta")$speciesKey

#search data according to key. occ_data vs occ_Search
dat <- occ_data(taxonKey = key,limit=5000)
dat.data<- dat$data
dim(dat$data)

colnames(da)

##select the records that have longitude and latitude data
acgeo <- subset(dat.data, !is.na(decimalLongitude) & !is.na(decimalLatitude))
dim(acgeo)


#take all the points in the extent of study area
dat.extent <- dat.data[dat.data$decimalLongitude >= -43 & dat.data$decimalLatitude <= 109, ]

library(geodata)
wrld <- world(path=".")
plot(wrld, xlim=c(-10,60), ylim=c(10,82), col="light yellow", border="light gray")
# add the points
points(acgeo$decimalLongitude, acgeo$decimalLatitude, col='red', pch=20)

#Data cleaning

# differentiating by (sub) species
# dups2 <- duplicated(acgeo[, c('species', 'lon', 'lat')])
# ignoring (sub) species and other naming variation
dups2 <- duplicated(acgeo[, c('decimalLongitude', 'decimalLatitude')])
# number of duplicates
sum(dups2)
## [1] 483
# keep the records that are _not_ duplicated
acg <- acgeo[!dups2, ]

#Cross-checking
library(terra)


acv <- vect(acg, geom=c("decimalLongitude", "decimalLatitude"), crs="+proj=longlat +datum=WGS84")
class(acv)

ovr <- extract(acv, wrld)
head(ovr)

cntr <- ovr$NAME_0

i <- which(is.na(cntr))
i
## integer(0)
j <- which(cntr != acv$country)
# for the mismatches, bind the country names of the polygons and points
m <- cbind(cntr[j], acg$country[j])
colnames(m) <- c("polygons", "acaule")
m
## polygons acaule

plot(acv)
lines(wrld, col='blue', lwd=2)
points(acv[j, ], col='red', pch=20, cex=2)

georef <- subset(dat.data, (is.na(decimalLongitude) | is.na(decimalLatitude)) & ! is.na(locality) )
dim(georef)


#bias
# create a SpatRaster with the extent of acgeo
r <- rast(acv)
# set the resolution of the cells to (for example) 1 degree
res(r) <- 1
# extend (expand) the extent of the SpatRaster a little
r <- extend(r, ext(r)+1)

# sample:
set.seed(13)
acsel <- spatSample(acv, size=1, "random", strata=r)
# to illustrate the method and show the result
p <- as.polygons(r)
plot(p, border='gray')
points(acv)
# selected points in red
points(acsel, cex=1, col='red', pch='x')

file <- paste(system.file(package="predicts"), '/ex/acaule.rds', sep='')
acsel <- vect(readRDS(file))

library(predicts)
# get the predictors filename
f1 <- system.file("ex/bio.tif", package="predicts")
f2 <- system.file("ex/biome.tif", package="predicts")
r <- rast(c(f1, f2))
plot
# have the same random sample.
set.seed(1963)

bg <- spatSample(r, 1000, "random", na.rm=TRUE, as.points=TRUE)

plot(r, 1)
points(bg, cex=0.5)

e <- ext(-80, -53, -39, -22)
bg2 <- spatSample(r, 500, "random", na.rm=FALSE, as.points=TRUE, ext=e)
plot(r, 1)
lines(e, col="red", lwd=2)
points(bg2, cex=0.5)

acfile <- file.path(system.file(package="predicts"), "ex/acaule.rds")
ac <- readRDS(acfile)

ac <- vect(ac)

# circles with a radius of 50 km
x <- buffer(ac, 50000)
pol <- aggregate(x)
# sample randomly from all circles
set.seed(999)
samp1 <- spatSample(pol, 500, "random")

pcells <- cells(r, samp1)
# remote duplicates
pcells <- unique(pcells[,2])
# back to coordinates
xy <- xyFromCell(r, pcells)

plot(pol, axes=TRUE)
points(samp1, pch="+", cex=.5)
points(xy, cex=0.75, pch=20, col='blue')

spxy <- vect(xy, crs="+proj=longlat +datum=WGS84")
xyInside <- intersect(spxy, x)

# make a new, empty, smaller raster
m <- crop(rast(r[[1]]), ext(x)+1)
# extract cell numbers for the circles
v <- cells(m, x)
# get unique cell numbers from which you could sample
v <- unique(v[,"cell"])
head(v)
## [1] 1505 1541 918 919 954 955
# to display the results
m[v] <- 1
plot(m)
lines(x)



library(predicts)
# get the predictors filename
f1 <- system.file("ex/bio.tif", package="predicts")
f2 <- system.file("ex/biome.tif", package="predicts")
r <- rast(c(f1, f2))
# select 5000 random points
# set seed to assure that the examples will always
# have the same random sample.
set.seed(1963)
bg <- spatSample(r, 5000, "random", na.rm=FALSE, as.points=TRUE)

plot(r, 1)
points(bg, cex=0.5)


# read worldclim
files.present <- list.files('D:/Github/BGE-SDM/RDATA/clipped/', pattern="[.]asc$", full.names=T) # alternatives for pattern (c|C)(e|E)(l|L)$
files.present
present.stack <- stack(files.present)
head(present.stack)

predictors <- rast(present.stack)
predictors

names(predictors)

library(geodata)
wrld <- world(path=".")

file <- paste("D:/Github/BGE-SDM/Output/fe.gbif.maxent.csv", sep="")
file

bradypus <- read.table(file, header=TRUE, sep=',')
# we do not need the first column
bradypus <- bradypus[,-1]

plot(predictors, 1)
lines(wrld)
points(bradypus, col='blue')

#Extracting values from rasters
presvals <- extract(predictors, bradypus)
# remove the ID variable
presvals <- presvals[,-1]
# setting random seed to always create the same
# random set of points for this example
set.seed(0)
backgr <- spatSample(predictors, 500, "random", as.points=TRUE)
absvals <- values(backgr)

pb <- c(rep(1, nrow(presvals)), rep(0, nrow(absvals)))
sdmdata <- data.frame(cbind(pb, rbind(presvals, absvals)))
head(sdmdata)

pairs(sdmdata[,2:5], cex=0.1)

saveRDS(sdmdata, "sdm.Rds")
saveRDS(presvals, "pvals.Rds")

#Model fitting
sdmdata <- readRDS("sdm.Rds")
presvals <- readRDS("pvals.Rds")

m1 <- glm(pb ~ fe_buffer_bio01 + fe_buffer_bio05 + fe_buffer_bio12, data=sdmdata)
class(m1)
summary(m1)

m2 = glm(pb ~ ., data=sdmdata)
m2

bc <- envelope(presvals[,c("fe_buffer_bio01", "fe_buffer_bio05", "fe_buffer_bio12")])
bc

#Model prediction

names(predictors)

p <- predict(predictors, m1)
plot(p)


samp <- sample(nrow(sdmdata), round(0.75 * nrow(sdmdata)))
traindata <- sdmdata[samp,]
traindata <- traindata[traindata[,1] == 1, 2:9]
testdata <- sdmdata[-samp,]
bc <- envelope(traindata)
e <- pa_evaluate(testdata[testdata==1,], testdata[testdata==0,], bc)
e

plot(e, "ROC")

pres <- sdmdata[sdmdata[,1] == 1, 2:9]
back <- sdmdata[sdmdata[,1] == 0, 2:9]

k <- 5
group <- make_folds(pres, k)
group[1:10]

e <- list()
for (i in 1:k) {
  train <- pres[group != i,]
  test <- pres[group == i,]
  bc <- envelope(train)
  p <- predict(bc, test)
  a <- predict(bc, back)
  e[[i]] <- pa_evaluate(p, a)
}

stats <- lapply(e, function(x){x@stats})
stats <- do.call(rbind, stats)
colMeans(stats)

bradypus <- bradypus[,-1]
presvals <- extract(predictors, bradypus)
set.seed(0)
backgr <- spatSample(predictors, 500, xy=TRUE, values=FALSE)

nr <- nrow(bradypus)
s <- sample(nr, 0.25 * nr)
pres_train <- bradypus[-s, ]
pres_test <- bradypus[s, ]

nr <- nrow(backgr)
set.seed(9)
s <- sample(nr, 0.25 * nr)
back_train <- backgr[-s, ]
back_test <- backgr[s, ]

sb <- dismo::ssb(pres_test, back_test, pres_train)
sb[,1] / sb[,2]


i <- pwd_sample(pres_test, back_test, pres_train, n=1, tr=0.1)
pres_test_pwd <- pres_test[!is.na(i[,1]), ]
back_test_pwd <- back_test[na.omit(as.vector(i)), ]
sb2 <- dismo::ssb(pres_test_pwd, back_test_pwd, pres_train)
sb2[1]/ sb2[2]
