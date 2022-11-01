install.packages(c("terra", "remotes"))

library(predicts)
filename <- file.path(system.file(package="predicts"), "ex/bradypus.csv")
# this is the file we will use:
basename(filename)


bradypus <- read.csv(filename)
# first rows
head(bradypus)

bradypus <- bradypus[,2:3]
head(bradypus)

acaule2 <- geodata::sp_occurrence("solanum", "acaule", geo=FALSE)

# load the saved S. acaule data
acfile <- file.path(system.file(package="predicts"), "ex/acaule.csv")
acaule <- read.csv(acfile)
dim(acaule)

colnames(acaule)
acgeo <- subset(acaule, !is.na(lon) & !is.na(lat))
dim(acgeo)

library(geodata)
wrld <- world(path=".")
plot(wrld, xlim=c(-110,60), ylim=c(-80,40), col="light yellow", border="light gray")
# add the points
points(acgeo$lon, acgeo$lat, col='red', pch=20)

library(predicts)
# get the predictors filename
f1 <- system.file("ex/bio.tif", package="predicts")
f2 <- system.file("ex/biome.tif", package="predicts")
r <- rast(c(f1, f2))
# select 5000 random points
# set seed to assure that the examples will always
# have the same random sample.
set.seed(1963)
bg <- spatSample(r, 5000, "random", as.points=TRUE)

plot(r, 1)
points(bg, cex=0.5)

library(predicts)
f <- system.file("ex/bio.tif", package="predicts")
predictors <- rast(f)
predictors
names(predictors)
plot(predictors)

library(geodata)
wrld <- world(path=".")
file <- paste(system.file(package="predicts"), "/ex/bradypus.csv", sep="")
bradypus <- read.table(file,  header=TRUE,  sep=',')
# we do not need the first column
bradypus  <- bradypus[,-1]

# first layer of the SpatRaster
plot(predictors, 1)
lines(wrld)
points(bradypus, col='blue')

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

sdmdata <- readRDS("sdm.Rds")
presvals <- readRDS("pvals.Rds")

m1 <- glm(pb ~ bio1 + bio5 + bio12, data=sdmdata)
class(m1)

summary(m1)
m2 = glm(pb ~ ., data=sdmdata)


bc <- envelope(presvals[,c("bio1", "bio5", "bio12")])
bc

bio1 <- c(40, 150, 200)
bio5 <- c(60, 115, 290)
bio12 <- c(600, 1600, 1700)
pd <- data.frame(cbind(bio1, bio5, bio12))
pd

pr <- partialResponse(bc, presvals, "bio1")
plot(pr, type="l")


f <- system.file("ex/bio.tif", package="predicts")
predictors <- rast(f)
names(predictors)
## [1] "bio1"  "bio5"  "bio6"  "bio7"  "bio8"  "bio9"  "bio12" "bio16" "bio17"
p <- predict(predictors, m1)
plot(p)

p <- rnorm(50, mean=0.7, sd=0.3)
a <- rnorm(50, mean=0.4, sd=0.4)
par(mfrow=c(1, 2))
plot(sort(p), col='red', pch=21)
points(sort(a), col='blue', pch=24)
legend(1, 0.95 * max(a,p), c('presence', 'absence'),
       pch=c(21,24), col=c('red', 'blue'))
comb <- c(p,a)
group <- c(rep('presence', length(p)), rep('absence', length(a)))
boxplot(comb~group, col=c('blue', 'red'))

group = c(rep(1, length(p)), rep(0, length(a)))
cor.test(comb, group)$estimate
##       cor
## 0.3418309
mv <- wilcox.test(p,a)
auc <- as.numeric(mv$statistic) / (length(p) * length(a))
auc

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

fsp <- system.file("/ex/bradypus.csv", package="predicts")
bradypus <- read.csv(fsp)
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


i <- pwd_sample(pres_test, back_test, pres_train, n=1, tr=0.1)
pres_test_pwd <- pres_test[!is.na(i[,1]), ]
back_test_pwd <- back_test[na.omit(as.vector(i)), ]
sb2 <- dismo::ssb(pres_test_pwd, back_test_pwd, pres_train)

bc <- envelope(predictors, pres_train)
ptst <- predict(bc, extract(predictors, pres_test))
btst <- predict(bc, extract(predictors, back_test))
pa_evaluate(ptst, btst)@stats

pwdptst <- predict(bc, extract(predictors, pres_test_pwd))
pwdbtst <- predict(bc, extract(predictors, back_test_pwd))
pa_evaluate(pwdptst, pwdbtst)@stats
