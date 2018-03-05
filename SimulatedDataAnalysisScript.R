rm(list=ls())

##################################################################
###
### Glacier Bay sea otter analysis script
###
##################################################################


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
                    "RColorBrewer",
                    "RcppArmadillo",
                    "rgdal",
                    "rgeos")
## install.packages(required.packages)
lapply(required.packages,library,character.only=TRUE)

## C++ sampler for projection
    code <- '
    arma::mat Hmat = Rcpp::as<arma::mat>(H);
    arma::mat c0mat = Rcpp::as<arma::mat>(c0);
    int steps = Rcpp::as<int>(timesteps);
    int n = steps;
    int k = Hmat.n_rows;
    arma::mat call(k,n);
    call.col(0)=Hmat*c0mat;
    for(int i = 1; i < n; ++i){
      call.col(i)=Hmat*call.col(i-1);
    }
    return Rcpp::wrap(call);
  '
    calcc = cxxfunction(signature(H="numeric",c0="numeric",timesteps="numeric"),
                        body=code, plugin="RcppArmadillo")

## Neighborhood matrix calculator
neighborhood=function(raster, boundary){
    nn=matrix(,length(raster[]),4)
    for(i in 1:dim(nn)[1]){
        if(extract(boundary,i)==1){
            loc=adjacent(raster,i)[,2]
            ln=loc[which((loc+1)==i)]
            rn=loc[which((loc-1)==i)]
            bn=loc[which((loc-dim(raster)[2])==i)]
            tn=loc[which((loc+dim(raster)[2])==i)]
            nn[i,1]=if(length(ln)>0 & extract(boundary,ln)>0){ln}else{0}
            nn[i,2]=if(length(rn)>0 & extract(boundary,rn)>0){rn}else{0}
            nn[i,3]=if(length(bn)>0 & extract(boundary,bn)>0){bn}else{0}
            nn[i,4]=if(length(tn)>0 & extract(boundary,tn)>0){tn}else{0}
        }else{nn[i,]=0}
    }
    nn
}

## Propagator matrix for plain diffusion PDE
propagator.plainZF=function(NN,delta,gamma,dx,dy,dt){
    H <- matrix(0,dim(NN)[1],dim(NN)[1])
    for(i in 1:dim(H)[1]){
        if(length(which(NN[i,]>0))>0){
            ind.tmp=ifelse(NN[i,]>0,1,0)
            H[i,i]=1-2*delta[i]*
                (dt/dx^2+dt/dy^2)+dt*gamma[i]+
                dt/dx^2*delta[i]*(1-ind.tmp[1])+
                          dt/dx^2*delta[i]*(1-ind.tmp[2])+
                                    dt/dy^2*delta[i]*(1-ind.tmp[3])+
                                              dt/dy^2*delta[i]*
                                                        (1-ind.tmp[4])
            H[i,NN[i,1]]=dt/dx^2*delta[i]*ind.tmp[1]
            H[i,NN[i,2]]=dt/dx^2*delta[i]*ind.tmp[2]
            H[i,NN[i,3]]=dt/dy^2*delta[i]*ind.tmp[3]
            H[i,NN[i,4]]=dt/dy^2*delta[i]*ind.tmp[4]
        }else{H[i,i]=0}
    }
    H
}

###
### Load Glacier Bay Sea Otter Data
###

setwd("~/Dropbox/Post-Doc/OriginalDataFiles")
OrigData.tmp=read.csv("GLBA_noDots.csv",header=TRUE)
OrigData=subset(OrigData.tmp,OrigData.tmp$group!=0)

###
### Get google map
###

map=get_map(
    location = c(-136.5,58.75),
    source = "google",
    zoom = 9,
    maptype = "terrain",
    color="bw"
)
background.map=ggmap(map)
## background.map

###
### Subset data
###

year=OrigData$year[!is.na(OrigData$GROUP_X)]
counts=OrigData$adults[!is.na(OrigData$GROUP_X)]+
    OrigData$pups[!is.na(OrigData$GROUP_X)]
x=OrigData$GROUP_X[!is.na(OrigData$GROUP_X)]
y=OrigData$GROUP_Y[!is.na(OrigData$GROUP_X)]

###
### Describe projection
###

utmcoor<-SpatialPoints(cbind(x,y),
                       proj4string=CRS(
                           "+proj=utm +zone=8  +datum=NAD27 +units=m")
                       )
lon.UTM=utmcoor$x
lat.UTM=utmcoor$y
loc.data.UTM=data.frame(lat.UTM,lon.UTM)

###
### Convert coordinates to Decimal Degrees (lat,long)
###

longlatcoor<-spTransform(utmcoor,CRS("+proj=longlat +datum=NAD83"))
lon.DD=longlatcoor$x
lat.DD=longlatcoor$y
loc.data.DD=data.frame(lat.DD,lon.DD)


###
### Plot map and all locations
###

location.map=background.map +
    geom_point(data=loc.data.DD,
               aes(y=lat.DD,
                   x=lon.DD,
                   color=log(counts)
                   ),
               size=0.1
               ) +
    scale_colour_gradientn(colours = terrain.colors(10)) +
    theme(legend.position="none")
## location.map


###
### 1999 locations
###

p.1999=background.map +
    geom_point(data=subset(loc.data.DD,year==1999),
               aes(y=lat.DD,
                   x=lon.DD,
                   color=log(subset(counts,year==1999))
                   ),
               size=0.1
               ) +
    scale_colour_gradientn(colours = terrain.colors(10)) +
    theme(legend.position="none")
## p.1999

###
### 2000 locations
###

p.2000=background.map +
    geom_point(data=subset(loc.data.DD,year==2000),
               aes(y=lat.DD,
                   x=lon.DD,
                   color=log(subset(counts,year==2000))
                   ),
               size=0.1
               ) +
    scale_colour_gradientn(colours = terrain.colors(10)) +
    theme(legend.position="none")
## p.2000

###
### 2001 locations
###

p.2001=background.map +
    geom_point(data=subset(loc.data.DD,year==2001),
               aes(y=lat.DD,
                   x=lon.DD,
                   color=log(subset(counts,year==2001))
                   ),
               size=0.1
               ) +
  scale_colour_gradientn(colours = terrain.colors(10)) +
    theme(legend.position="none")
## p.2001

###
### 2002 locations
###

p.2002=background.map +
    geom_point(data=subset(loc.data.DD,year==2002),
               aes(y=lat.DD,
                   x=lon.DD,
                   color=log(subset(counts,year==2002))
                   ),
               size=0.1
               ) +
  scale_colour_gradientn(colours = terrain.colors(10)) +
    theme(legend.position="none")
## p.2002

###
### 2003 locations
###

p.2003=background.map +
    geom_point(data=subset(loc.data.DD,year==2003),
               aes(y=lat.DD,
                   x=lon.DD,
                   color=log(subset(counts,year==2003))
                   ),
               size=0.1
               ) +
  scale_colour_gradientn(colours = terrain.colors(10)) +
    theme(legend.position="none")
## p.2003

###
### 2004 locations
###

p.2004=background.map +
    geom_point(data=subset(loc.data.DD,year==2004),
               aes(y=lat.DD,
                   x=lon.DD,
                   color=log(subset(counts,year==2004))
                   ),
               size=0.1
               ) +
  scale_colour_gradientn(colours = terrain.colors(10)) +
    theme(legend.position="none")
## p.2004

###
### 2006 locations
###

p.2006=background.map +
    geom_point(data=subset(loc.data.DD,year==2006),
               aes(y=lat.DD,
                   x=lon.DD,
                   color=log(subset(counts,year==2006))
                   ),
               size=0.1
               ) +
  scale_colour_gradientn(colours = terrain.colors(10)) +
    theme(legend.position="none")
## p.2006

###
### 2012 locations
###

p.2012=background.map +
    geom_point(data=subset(loc.data.DD,year==2012),
               aes(y=lat.DD,
                   x=lon.DD,
                   color=log(subset(counts,year==2012))
                   ),
               size=0.1
               ) +
  scale_colour_gradientn(colours = terrain.colors(10)) +
    theme(legend.position="none")
## p.2012

###
### Plot all years with google maps
###

## grid.arrange(p.1999,p.2000,p.2001,p.2002,p.2003,p.2004,
##              p.2006,p.2012, ncol=2)

###
### Load ISU data
###

all.ISU.data=read.csv(paste("~/Dropbox/Post-Doc/",
                        "OriginalDataFiles/",
                        "All ISUs 1999_2012_02032016.csv",
                        sep=""))
ISU.N=all.ISU.data$circle.adt+
    all.ISU.data$circle.pup
ISU.Y=all.ISU.data$strip.adt+
    all.ISU.data$strip.pup
ISU.N[ISU.Y>ISU.N]=ISU.Y[ISU.Y>ISU.N]
sum(ISU.N<ISU.Y)
ISU.year=all.ISU.data$year
ISU.data=data.frame(ISU.N,ISU.Y,ISU.year)
Y.1999=ISU.data$ISU.Y[ISU.data$ISU.year==1999]
N.1999=ISU.data$ISU.N[ISU.data$ISU.year==1999]
Y.2000=ISU.data$ISU.Y[ISU.data$ISU.year==2000]
N.2000=ISU.data$ISU.N[ISU.data$ISU.year==2000]
Y.2001=ISU.data$ISU.Y[ISU.data$ISU.year==2001]
N.2001=ISU.data$ISU.N[ISU.data$ISU.year==2001]
Y.2002=ISU.data$ISU.Y[ISU.data$ISU.year==2002]
N.2002=ISU.data$ISU.N[ISU.data$ISU.year==2002]
Y.2003=ISU.data$ISU.Y[ISU.data$ISU.year==2003]
N.2003=ISU.data$ISU.N[ISU.data$ISU.year==2003]
Y.2004=ISU.data$ISU.Y[ISU.data$ISU.year==2004]
N.2004=ISU.data$ISU.N[ISU.data$ISU.year==2004]
Y.2006=ISU.data$ISU.Y[ISU.data$ISU.year==2006]
N.2006=ISU.data$ISU.N[ISU.data$ISU.year==2006]
Y.2012=ISU.data$ISU.Y[ISU.data$ISU.year==2012]
N.2012=ISU.data$ISU.N[ISU.data$ISU.year==2012]

ISU=list(Y.1999=Y.1999,
         N.1999=N.1999,
         Y.2000=Y.2000,
         N.2000=N.2000,
         Y.2001=Y.2001,
         N.2001=N.2001,
         Y.2002=Y.2002,
         N.2002=N.2002,
         Y.2003=Y.2003,
         N.2003=N.2003,
         Y.2004=Y.2004,
         N.2004=N.2004,
         Y.2006=Y.2006,
         N.2006=N.2006,
         Y.2012=Y.2012,
         N.2012=N.2012)

#############################################
### Load Bathymetry Data
#############################################

load("~/Dropbox/Post-Doc/RFiles/lutris/data/bath.raster.RData")

###
### Study area extent
###

xmin=404000
xmax=460000
ymin=6460000
ymax=6528000

###
### Bathymetry data
###

bath.tmp=crop(bath.raster,extent(xmin,xmax,ymin,ymax))
Bath.r=aggregate(bath.tmp,c(16,16),na.rm=TRUE)
Bath.r[94,92]=NA  # add boulder island

###
### Plot Bath.r and sea otter locations
###

## plot(Bath.r)
## points(loc.data.UTM$lon,loc.data.UTM$lat,cex=0.2,pch=16)


#################################################
#################################################
###
### Simulated Data for Analysis
###
#################################################
#################################################

###
### Parameter values
###

theta=500
kappa=5.95
beta0=17.72
beta1=-1.48
beta2=0.77
beta3=-0.35
beta4=1.0
beta=c(beta0,beta1,beta2,beta3,beta4)
gamma=0.20
p.1999=0.7953388
p.2000=0.75
p.2001=0.8583887
p.2002=0.8503951
p.2003=0.7606767
p.2004=0.7694584
p.2006=0.7470478
p.2012=0.5833859
odp=0.5

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
### Aggregate values
###

bath=crop(Bath.r,extent(xmin,xmax,ymin,ymax))
## plot(bath)
## points(loc.data.UTM$lon,loc.data.UTM$lat,cex=0.2,pch=16)

###
### Dimensions
###

x=dim(bath)[2]
y=dim(bath)[1]
q=x*y

###
### Create `Cell' raster
###

cell=raster(,nrows=y,ncols=x,xmn=xmin,xmx=xmax,
               ymn=ymin,ymx=ymax,crs=NA)
cell[]=1:q
## plot(cell)

###
### Create Bounday raster
###

Boundary=bath
Boundary[is.na(Bath.r)]=0
Boundary[Bath.r<0]=1
Boundary[Bath.r>=0]=0
## plot(Boundary)
## points(loc.data.UTM$lon,loc.data.UTM$lat,cex=0.2,pch=16)

BoundaryNA=Boundary
BoundaryNA[Boundary==0]=NA
## plot(BoundaryNA)
## points(loc.data.UTM$lon,loc.data.UTM$lat,cex=0.2,pch=16)

BoundaryInf=BoundaryNA
BoundaryInf[120:170,0:80]=NA
BoundaryInf[150:170,0:140]=NA
BoundaryInf[134:170,116:140]=NA
BoundaryInf.v=rep(BoundaryInf[],length(keep))
## plot(BoundaryInf)

###
### Up-scaled boundary layer surrounded by zeros
###

Boundary.us=aggregate(Boundary,fact=us.fact,na.rm=TRUE,fun=max)
Boundary.us[Boundary.us[]>1]=1
Boundary.us[1,]=0
Boundary.us[,1]=0
Boundary.us[dim(Boundary.us)[1],]=0
Boundary.us[,dim(Boundary.us)[2]]=0
 ## plot(Boundary.us)

###
### Depth covariate
###

depth=40
DepthCov=bath
DepthCov[bath[]< -depth]=0
DepthCov[bath[]>= -depth]=1
DepthCov[is.na(bath[])] = 0
DepthCov=DepthCov*Boundary
## plot(DepthCov)
## points(loc.data.UTM$lon,loc.data.UTM$lat,cex=0.2,pch=16)

###
### Distance to shore covariate
###

DistCov.tmp1=Boundary
DistCov.tmp1[Boundary==1]=NA
DistCov.tmp2=distance(DistCov.tmp1)
DistCov.tmp2[is.na(DistCov.tmp2)]=0
DistCov=scale(DistCov.tmp2,center=TRUE,scale=TRUE)
DistCov=DistCov*Boundary
## plot(DistCov)
## points(loc.data.UTM$lon,loc.data.UTM$lat,cex=0.2,pch=16)

###
### Bottom slope
###

SlopeCov.tmp=terrain(bath,opt='slope',unit='degrees',neighbors=8)
SlopeCov.tmp[is.na(SlopeCov.tmp)]=0
SlopeCov=scale(SlopeCov.tmp, center=TRUE, scale=TRUE)
SlopeCov=SlopeCov*Boundary*DepthCov
## plot(SlopeCov)
## points(loc.data.UTM$lon,loc.data.UTM$lat,cex=0.2,pch=16)

###
### Shoreline complexity
###

ShoreCov.tmp=boundaries(DistCov.tmp1,type='inner')
ShoreCov.tmp[is.na(ShoreCov.tmp)]=0
ShoreCov=focal(ShoreCov.tmp,w=matrix(1,nr=11,nc=11))*Boundary
ShoreCov[is.na(ShoreCov)]=0
ShoreCov=scale(ShoreCov, center=TRUE, scale=TRUE)
ShoreCov=ShoreCov*Boundary
## plot(ShoreCov)
## points(loc.data.UTM$lon,loc.data.UTM$lat,cex=0.2,pch=16)

###
### Create co-variate matrices
###

X=cbind(1,DepthCov[],DistCov[],SlopeCov[]*DepthCov[],ShoreCov[])
W=matrix(1,nr=q,nc=1)

###
### Diffusion rate for pde
###

delta=cell
delta[]=exp(X%*%beta)
## plot(delta*BoundaryNA)

###
### Growth rate for pde
###

gamma[]=W%*%gamma

###
### Diffusion rate for homogenized pde
###

delta.bar=aggregate(delta,fact=us.fact,na.rm=TRUE,fun=function(x,na.rm){1/mean(1/x)})
## plot(delta.bar)

###
### Growth rate for homogenized pde
###

gamma.bar=delta.bar*aggregate(gamma/delta,fact=us.fact,fun=mean,na.rm=TRUE)
## plot(gamma.bar)

###
### First-order neighborhood matrix
###

NN=neighborhood(delta.bar,Boundary.us)

###
### Propagator matrix
###

H=propagator.plainZF(NN=NN,delta=delta.bar[],gamma=gamma.bar[],
                     dx=res(delta.bar)[1],dy=res(delta.bar)[2],dt=dt)

###
### Initial condition (D is in km)
###

c0=raster(,nrows=y/us.fact,ncols=x/us.fact,
            xmn=xmin,xmx=xmax,ymn=ymin,ymx=ymax,crs=NA)
lambda0=cell
D=rdist(data.frame(SpatialPoints(lambda0)),matrix(d,1,2))/1000
lambda0[]=(exp(-D^2/kappa^2)/sum(exp(-D^2/kappa^2))*theta)
c0[]=extract(delta*lambda0,SpatialPoints(delta.bar))
c0.star=c0
## plot(lambda0*BoundaryInf)
## plot(c0*Boundary.us)

###
### Stack raster layers
###

c.all=brick(nrows=dim(c0)[1], ncols=dim(c0)[2], xmn=xmin, xmx=xmax,
                ymn=ymin, ymx=ymax)
c.all=setValues(c.all, calcc(H, vec(c0[]),time.steps)[,keep])
lambda.all=disaggregate(c.all,us.fact)/delta
## plot(c.all)
## plot(lambda.all*BoundaryNA)
## plot(lambda.all[[20]]*BoundaryInf)
## points(loc.data.UTM$lon,loc.data.UTM$lat,cex=0.2,pch=16)

###
### Expected abundance
###

EA.v=numeric(20)
for(i in 1:20){
    EA.v[i]=sum(BoundaryInf[]*lambda.all[[i]][],na.rm=TRUE)
}
EA.v

###
### Generate abundance from expected abundance (lambda)
###

N=lambda0
N[]=rnbinom(n=length(lambda.all[[1]][]),size=odp,,mu=lambda.all[[1]][])
N=stack(mget(rep("N",length(keep))))
for(i in 2:length(keep)){
    N[[i]][]=rnbinom(n=length(lambda.all[[i]][]),size=odp,,mu=lambda.all[[i]][])
}

###
### Generate data from abundance
###

p=0.74
p.vs=c(p,
       p,
       p,
       p,
       p,
       p,
       p.1999,
       p.2000,
       p.2001,
       p.2002,
       p.2003,
       p.2004,
       p,
       p.2006,
       p,
       p,
       p,
       p,
       p,
       p.2012)
Y=N[[1]]
Y[]=rbinom(q, N[[1]][],p.vs[1])
Y=stack(mget(rep("Y", length(keep))))
for (i in 2:length(keep)) {
    Y[[i]][] = rbinom(q, N[[i]][],p.vs[i])
}

###
### Create spatial lines object
###

## 1999
TP.surv.1999=subset(OrigData.tmp,OrigData.tmp$year==1999)
l.1999=list()
Sl.1999=list()
S.1999=list()

for(i in 1:dim(TP.surv.1999)[1]){
    l.1999[[i]]=cbind(c(TP.surv.1999$west_long[i],
                        TP.surv.1999$east_long[i]),
                      c(TP.surv.1999$west_lat[i],
                        TP.surv.1999$east_lat[i]))
    Sl.1999[[i]]=Line(l.1999[[i]])
    S.1999[[i]]=Lines(list(Sl.1999[[i]]), ID=i)
}
Sb.1999=SpatialLines(S.1999)
## plot(Sb.1999)

## 2000
TP.surv.2000=subset(OrigData.tmp,OrigData.tmp$year==2000)
l.2000=list()
Sl.2000=list()
S.2000=list()

for(i in 1:dim(TP.surv.2000)[1]){
    l.2000[[i]]=cbind(c(TP.surv.2000$west_long[i],
                        TP.surv.2000$east_long[i]),
                      c(TP.surv.2000$west_lat[i],
                        TP.surv.2000$east_lat[i]))
    Sl.2000[[i]]=Line(l.2000[[i]])
    S.2000[[i]]=Lines(list(Sl.2000[[i]]), ID=i)
}
Sb.2000=SpatialLines(S.2000)
## plot(Sb.2000)

## 2001
TP.surv.2001=subset(OrigData.tmp,OrigData.tmp$year==2001)
l.2001=list()
Sl.2001=list()
S.2001=list()

for(i in 1:dim(TP.surv.2001)[1]){
    l.2001[[i]]=cbind(c(TP.surv.2001$west_long[i],
                        TP.surv.2001$east_long[i]),
                      c(TP.surv.2001$west_lat[i],
                        TP.surv.2001$east_lat[i]))
    Sl.2001[[i]]=Line(l.2001[[i]])
    S.2001[[i]]=Lines(list(Sl.2001[[i]]), ID=i)
}
Sb.2001=SpatialLines(S.2001)
## plot(Sb.2001)

## 2002
TP.surv.2002=subset(OrigData.tmp,OrigData.tmp$year==2002)
l.2002=list()
Sl.2002=list()
S.2002=list()

for(i in 1:dim(TP.surv.2002)[1]){
    l.2002[[i]]=cbind(c(TP.surv.2002$west_long[i],
                        TP.surv.2002$east_long[i]),
                      c(TP.surv.2002$west_lat[i],
                        TP.surv.2002$east_lat[i]))
    Sl.2002[[i]]=Line(l.2002[[i]])
    S.2002[[i]]=Lines(list(Sl.2002[[i]]), ID=i)
}
Sb.2002=SpatialLines(S.2002)
## plot(Sb.2002)

## 2003
TP.surv.2003=subset(OrigData.tmp,OrigData.tmp$year==2003)
l.2003=list()
Sl.2003=list()
S.2003=list()

for(i in 1:dim(TP.surv.2003)[1]){
    l.2003[[i]]=cbind(c(TP.surv.2003$west_long[i],
                        TP.surv.2003$east_long[i]),
                      c(TP.surv.2003$west_lat[i],
                        TP.surv.2003$east_lat[i]))
    Sl.2003[[i]]=Line(l.2003[[i]])
    S.2003[[i]]=Lines(list(Sl.2003[[i]]), ID=i)
}
Sb.2003=SpatialLines(S.2003)
## plot(Sb.2003)

## 2004
TP.surv.2004=subset(OrigData.tmp,OrigData.tmp$year==2004)
l.2004=list()
Sl.2004=list()
S.2004=list()

for(i in 1:dim(TP.surv.2004)[1]){
    l.2004[[i]]=cbind(c(TP.surv.2004$west_long[i],
                        TP.surv.2004$east_long[i]),
                      c(TP.surv.2004$west_lat[i],
                        TP.surv.2004$east_lat[i]))
    Sl.2004[[i]]=Line(l.2004[[i]])
    S.2004[[i]]=Lines(list(Sl.2004[[i]]), ID=i)
}
Sb.2004=SpatialLines(S.2004)
## plot(Sb.2004)

## 2006
TP.surv.2006=subset(OrigData.tmp,OrigData.tmp$year==2006)
l.2006=list()
Sl.2006=list()
S.2006=list()

for(i in 1:dim(TP.surv.2006)[1]){
    l.2006[[i]]=cbind(c(TP.surv.2006$west_long[i],
                        TP.surv.2006$east_long[i]),
                      c(TP.surv.2006$west_lat[i],
                        TP.surv.2006$east_lat[i]))
    Sl.2006[[i]]=Line(l.2006[[i]])
    S.2006[[i]]=Lines(list(Sl.2006[[i]]), ID=i)
}
Sb.2006=SpatialLines(S.2006)
## plot(Sb.2006)

## 2012
TP.surv.2012=subset(OrigData.tmp,OrigData.tmp$year==2012)
l.2012=list()
Sl.2012=list()
S.2012=list()

for(i in 1:dim(TP.surv.2012)[1]){
    l.2012[[i]]=cbind(c(TP.surv.2012$west_long[i],
                        TP.surv.2012$east_long[i]),
                      c(TP.surv.2012$west_lat[i],
                        TP.surv.2012$east_lat[i]))
    Sl.2012[[i]]=Line(l.2012[[i]])
    S.2012[[i]]=Lines(list(Sl.2012[[i]]), ID=i)
}
Sb.2012=SpatialLines(S.2012)
## plot(Sb.2012)

## plot(Boundary)
## plot(Sb.2012,add=TRUE)
loc.data.UTM.2012=subset(loc.data.UTM,year==2012)
## points(loc.data.UTM.2012$lon,loc.data.UTM.2012$lat,cex=0.2,pch=16)

###
### Convert the spatial lines object Sb into a raster
###

transects.r=raster(,nrows=y,ncols=x,xmn=xmin,
                   xmx=xmax,ymn=ymin,ymx=ymax)

## Transects.r.1993=rasterize(Sb.1993,transects.r)
Transects.r.1999=rasterize(Sb.1999,transects.r)
Transects.r.2000=rasterize(Sb.2000,transects.r)
Transects.r.2001=rasterize(Sb.2001,transects.r)
Transects.r.2002=rasterize(Sb.2002,transects.r)
Transects.r.2003=rasterize(Sb.2003,transects.r)
Transects.r.2004=rasterize(Sb.2004,transects.r)
Transects.r.2006=rasterize(Sb.2006,transects.r)
Transects.r.2012=rasterize(Sb.2012,transects.r)

## Transects.r.1993[!is.na(Transects.r.1993[])]=1
Transects.r.1999[!is.na(Transects.r.1999[])]=1
Transects.r.2000[!is.na(Transects.r.2000[])]=1
Transects.r.2001[!is.na(Transects.r.2001[])]=1
Transects.r.2002[!is.na(Transects.r.2002[])]=1
Transects.r.2003[!is.na(Transects.r.2003[])]=1
Transects.r.2004[!is.na(Transects.r.2004[])]=1
Transects.r.2006[!is.na(Transects.r.2006[])]=1
Transects.r.2012[!is.na(Transects.r.2012[])]=1

###
### Capture simulated Y at transects
###

## Data.1993.r=Transects.r.1993
Data.1999.r=Transects.r.1999
Data.2000.r=Transects.r.2000
Data.2001.r=Transects.r.2001
Data.2002.r=Transects.r.2002
Data.2003.r=Transects.r.2003
Data.2004.r=Transects.r.2004
Data.2006.r=Transects.r.2006
Data.2012.r=Transects.r.2012

## Data.1993.r[]=Transects.r.1993[]*Y[[1]][]
Data.1999.r[]=Transects.r.1999[]*Y[[7]][]
Data.2000.r[]=Transects.r.2000[]*Y[[8]][]
Data.2001.r[]=Transects.r.2001[]*Y[[9]][]
Data.2002.r[]=Transects.r.2002[]*Y[[10]][]
Data.2003.r[]=Transects.r.2003[]*Y[[11]][]
Data.2004.r[]=Transects.r.2004[]*Y[[12]][]
Data.2006.r[]=Transects.r.2006[]*Y[[14]][]
Data.2012.r[]=Transects.r.2012[]*Y[[20]][]
## plot(Data.1999.r)
## plot(Data.2012.r)

###
### Data for analysis
###

Y=c(## Data.1993.r[],
    rep(NA,length(Data.1999.r[])*(1999-1993)),
    Data.1999.r[],
    Data.2000.r[],
    Data.2001.r[],
    Data.2002.r[],
    Data.2003.r[],
    Data.2004.r[],
    rep(NA,length(Data.1999.r[])),
    Data.2006.r[],
    rep(NA,length(Data.1999.r[])*(2012-2007)),
    Data.2012.r[])

############################################################################
############################################################################
###
### Fit the model to the simulated data using MCMC
###
############################################################################
############################################################################

###
### Bundle data
###

data=list(Y=Y,
          ISU=ISU,
          X=X)

###
### Spatio-temporal settings
###

dt=1/200
time.frame=1993:2012
us.fact=10
res=400
d=c(442000,6465000)
xmin=404000
xmax=460000
ymin=6460000
ymax=6528000
extent=c(xmin,xmax,ymin,ymax)
st.info=list(dt=dt,
             time.frame=time.frame,
             us.fact=us.fact,
             res=res,
             d=d,
             extent=extent
             )

###
### MCMC Settings
###

n.iter=50000
checkpoint=10

###
### Priors
###

## beta
mu.beta=0
var.beta=10^2

## gamma
q.gamma=-0.5
r.gamma=0.5

## theta
mu.theta=500
var.theta=500

## kappa
mu.kappa=5.95
var.kappa=5

## p
q.p=1
r.p=1

## tau
q.tau=0
r.tau=1

priors=list(
    beta.prior=c(mu.beta,var.beta),
    gamma.prior=c(q.gamma,r.gamma),
    theta.prior=c(mu.theta,var.theta),
    kappa.prior=c(mu.kappa,var.kappa),
    p.prior=c(q.p,r.p),
    tau.prior=c(q.tau,r.tau)
)

###
### Starting values
###

gamma=0.20
beta=c(17.72,-1.48,0.77,-0.35,1.0)
theta=500
kappa=5.95
p.1999=0.80
p.2000=0.75
p.2001=0.86
p.2002=0.85
p.2003=0.76
p.2004=0.77
p.2006=0.75
p.2012=0.58
p=c(rep(NA,6*q),
    rep(p.1999,q),
    rep(p.2000,q),
    rep(p.2001,q),
    rep(p.2002,q),
    rep(p.2003,q),
    rep(p.2004,q),
    rep(NA,q),
    rep(p.2006,q),
    rep(NA,q*5),
    rep(p.2012,q))
odp=0.5

inits=list(gamma=gamma,
           beta=beta,
           theta=theta,
           kappa=kappa,
           odp=odp,
           p=p
           )

###
### Parameters to monitor
###

parameters=c("gamma",
             "beta",
             "kappa",
             "theta",
             "odp",
             "p",
             "n.tot")

###
### Output location
###

output.location="~/GB_Output.RData"

Bathymetry=bath
rm(list=ls()[!ls() %in% c("data",
                          "priors",
                          "inits",
                          "parameters",
                          "st.info",
                          "Bathymetry",
                          "n.iter",
                          "checkpoint",
                          "output.location")])

###
### Run MCMC
###

source(paste("~/Dropbox/GitHub/GlacierBay-SeaOtter-Analysis/",
             "GlacierBay-SeaOtter-Analysis/MCMC.R",sep=""))
MCMC(data,
     priors,
     inits,
     parameters,
     st.info,
     Bathymetry,
     n.iter=n.iter,
     checkpoint=checkpoint,
     output.location)



############################################################################
### Load MCMC results
############################################################################

rm(list=ls())
output.location="~/GB_Output.RData"
load(output.location)
(status=sum(!is.na(MCMC.Chains[[3]])))

###
### True values
###

beta0=17.72
beta1=-1.48
beta2=0.77
beta3=-0.35
beta4=1.0
p.truth=c(0.7953388,
          0.75,
          0.8583887,
          0.8503951,
          0.7606767,
          0.7694584,
          0.7470478,
          0.5833859)
odp=0.5

###
### Extract parameters
###

gamma.est=MCMC.Chains$gamma[1:status]
beta.est=MCMC.Chains$beta[1:status,]
theta.est=MCMC.Chains$theta[1:status]
kappa.est=MCMC.Chains$kappa[1:status]
odp.est=MCMC.Chains$odp[1:status]

###
### Plot parameters
###

par(mfrow=c(5,2),mar=c(4,4,1,1))
plot(gamma.est,type='l',ylim=c(0.15,0.25))
abline(h=0.2,col=2)
plot(beta.est[,1],type='l')
abline(h=beta0,col=2)
plot(beta.est[,2],type='l')
abline(h=beta1,col=2)
plot(beta.est[,3],type='l')
abline(h=beta2,col=2)
plot(beta.est[,4],type='l')
abline(h=beta3,col=2)
plot(beta.est[,5],type='l')
abline(h=beta4,col=2)
plot(theta.est,type='l')
abline(h=500,col=2)
plot(kappa.est,type='l')
abline(h=5.95,col=2)
plot(odp.est,type='l')
abline(h=0.5,col=2)

###
### Calculate abundance
###

N.est=MCMC.Chains$n.tot[1:status,]
mean.n=apply(N.est,2,mean,na.rm=TRUE)
lb.n=apply(N.est,2,quantile,0.025)
ub.n=apply(N.est,2,quantile,0.975)
cbind(lb.n,mean.n,ub.n)

###
### Plot abundance
###

par(mfrow=c(5,4),mar=c(2,2,0,1))
for(i in 1:20){
    plot(N.est[,i],type='l')
}

###
### Plot p
###

p.est=MCMC.Chains$p[1:status,]
for(i in 1:8){
    plot(p.est[,i],type='l');abline(h=p.truth[i],col=2)
}
