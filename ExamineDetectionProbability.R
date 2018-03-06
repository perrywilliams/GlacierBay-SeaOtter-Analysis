rm(list=ls())

############################################################################
############################################################################
############################################################################
###
### Real Data Analysis
###
############################################################################
############################################################################
############################################################################

###
### Packages and dependencies
###

required.packages=c("coda",
                    "fBasics",
                    "fields",
                    "ggmap",
                    "ggplot2",
                    "gridExtra",
                    "gstat",
                    "inline",
                    "maptools",
                    "raster",
                    "rasterVis",
                    "RCurl",
                    "RColorBrewer",
                    "RcppArmadillo",
                    "rgdal",
                    "rgeos")
lapply(required.packages,library,character.only=TRUE)

############################################################################
### Load Bathymetry Data
############################################################################

## This file is available on GitHub
load("~/Dropbox/Post-Doc/RFiles/lutris/data/bath.raster.RData")

###
### Study area extent
###

xmin=404000
xmax=460000
ymin=6460000
ymax=6528000

###
### Initial condition
###

d=c(442000,6465000)

###
### Discretization values
###

dt=1/200
time.steps=1/dt*length(1993:2012)
keep=seq(1,time.steps,time.steps/length(1993:2012))

###
### Homogenization value
###

us.fact=10

###
### Bathymetry data
###

bath.tmp=crop(bath.raster,extent(xmin,xmax,ymin,ymax))
Bath.r=aggregate(bath.tmp,c(16,16),na.rm=TRUE)
Bath.r[94,92]=NA  # add boulder island
x=dim(Bath.r)[2]
y=dim(Bath.r)[1]
q=x*y

###
### Aggregate values
###

bath=crop(Bath.r,extent(xmin,xmax,ymin,ymax))

###
### Convert Data to Raster
###

cell=raster(,nrows=y,ncols=x,xmn=xmin,xmx=xmax,
            ymn=ymin,ymx=ymax,crs=NA)
CISU=N.r=
    ISUind.r=
        Y.r=Counts=
            Boundary=
                BoundaryDist=
                    DistCov=
                        DepthCov=
                            gamma=delta=lambda0=SurvR=SimulatedDataR=cell


cell[]=1:q

rand.2017=read.csv("~/Dropbox/Post-Doc/OriginalDataFiles/SO_D_2017/SO_D_20170719_R.csv")
rand.2017=rand.2017[!is.na(rand.2017$LATITUDE_WGS84),]
xyc=cbind(rand.2017$LONGITUDE_WGS84,
          rand.2017$LATITUDE_WGS84)
DDcoor=SpatialPoints(xyc,CRS(
                           "+proj=longlat  +datum=WGS84")
                       )

utmcoor=spTransform(DDcoor,
                       CRS(
                           "+proj=utm +zone=8  +datum=NAD27 +units=m")
                       )
head(utmcoor@coords)
count.2017=rand.2017$COUNT_ADULT+rand.2017$COUNT_PUP
cell.ID=extract(cell,utmcoor)
mat=cbind(utmcoor@coords,cell.ID,count.2017)
head(mat)
mat=mat[order(mat[,3]),]
uniq.ID=unique(mat[,3])
max.photos=rep(NA,length(uniq.ID))
for(i in uniq.ID){
    data.tmp=subset(mat,cell.ID==i)
    max.photos[i]=dim(data.tmp)[1]
}
max.dim.Y2=max(max.photos,na.rm=TRUE)

Y.tmp=matrix(,length(uniq.ID),max.dim.Y2)
for(i in 1:length(uniq.ID)){
    data.tmp=subset(mat,mat[,3]==uniq.ID[i])
    Y.tmp[i,]=c(data.tmp[,4],rep(NA,dim(Y.tmp)[2]-length(data.tmp[,4])))
}
Y=cbind(uniq.ID,Y.tmp)

Y.2017=matrix(NA,q,max.dim.Y2)

Y.2017[uniq.ID,]=Y.tmp


###
###
###

## Priors
q.p=1
r.p=1
alpha=0.001
beta=0.001
n.iter=100000
checkpoint=5000
thin=1

###
### Run algorithm
###

name="~/Nmix.RData"
source("~/Dropbox/MCMCAlgorithms/NMixtureModel/SpatialVaryingNMixtureModel/MCMCAlgorithm.R")
Nmixmcmc(Y=Y.tmp,q.p,r.p,alpha,beta,n.iter,checkpoint,name,thin)

###
### Load output and plot results
###

load(name)
(status=sum(!is.na(out[[1]])))
p=out[[1]][1:status,]
plot(p,type='l')

N=out[[2]][1:status,]
head(N)
Y.tmp[452,]
plot(N[,452],type='l')
