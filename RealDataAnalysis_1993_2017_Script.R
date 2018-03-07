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

############################################################################
### Load Glacier Bay Sea Otter Data
############################################################################

data=read.csv("~/Dropbox/Post-Doc/OriginalDataFiles/GLBA_noDots.csv")
dataDistr=read.csv(paste("~/Dropbox/Post-Doc/OriginalDataFiles/",
                         "DistribSurveysC.csv",
                         sep=""))
ind.pred.tmp=1:dim(dataDistr)[1]
ind.pred=ind.pred.tmp[dataDistr$year==2004|dataDistr$year==2006]
dataDistr=dataDistr[-ind.pred,]
ind=1993:2012

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
c0=raster(,nrows=y/us.fact,ncols=x/us.fact,
          xmn=xmin,xmx=xmax,ymn=ymin,ymx=ymax,
          crs=NA)
data.r=raster(,nrows=y,ncols=x,xmn=xmin,
              xmx=xmax,ymn=ymin,ymx=ymax,crs=0)
data1.sub=subset(data,year=="1993")
data1.sub=data1.sub[!is.na(data1.sub$GROUP_X),]
data2.sub=subset(dataDistr,year=="1993")
x_loc=c(data1.sub$Group_X,data2.sub$POINT_X)
y_loc=c(data1.sub$Group_Y,data2.sub$POINT_Y)
Counts.tmp=c(data1.sub$animals,data2.sub$animals)
Counts.sub=rasterize(x=cbind(x_loc,y_loc),y=data.r,
                     field=Counts.tmp,fun="max",na.rm=TRUE,background=0)
Counts[]=Counts.sub[]
Counts=stack(mget(rep("Counts",25)))

years=1994:2012
ind=2
for(t in years){
    data1.sub=subset(data,year==t)
    data1.sub=data1.sub[!is.na(data1.sub$GROUP_X),]
    data2.sub=subset(dataDistr,year==t)

    x_loc=c(data1.sub$GROUP_X,data2.sub$POINT_X)
    y_loc=c(data1.sub$GROUP_Y,data2.sub$POINT_Y)
    Counts.tmp=c(data1.sub$animals,data2.sub$animals)
    if(length(x_loc)>0){
        Counts.sub=rasterize(x=cbind(x_loc,y_loc),
                             y=data.r,field=Counts.tmp,fun="max",na.rm=TRUE,
                             background=0)
        Counts[[ind]][]=Counts.sub[]
    }else{
        Counts[[ind]][]=rep(NA,q)
    }
    ind=ind+1
}

############################################################################
### Add new (2017) aerial photograph data
############################################################################

rand.2017=read.csv("~/Dropbox/Post-Doc/OriginalDataFiles/SO_D_2017/SO_D_20170719_R.csv")
rand.2017=rand.2017[!is.na(rand.2017$LATITUDE_WGS84),]
opt.2017=read.csv("~/Dropbox/Post-Doc/OriginalDataFiles/SO_D_2017/SO_D_20170721_O.csv")
opt.2017=opt.2017[!is.na(opt.2017$LATITUDE_WGS84),]
abund.2017=read.csv("~/Dropbox/Post-Doc/OriginalDataFiles/SO_D_2017/SO_D_20170728_A.csv")
abund.2017=abund.2017[!is.na(abund.2017$LATITUDE_WGS84),]

all.2017=rbind(rand.2017,opt.2017,abund.2017)
xyc=cbind(all.2017$LONGITUDE_WGS84,
          all.2017$LATITUDE_WGS84)

###
### Assign projection
###

DDcoor=SpatialPoints(xyc,CRS(
                           "+proj=longlat  +datum=WGS84")
                       )

###
### Change projection
###

utmcoor=spTransform(DDcoor,
                       CRS(
                           "+proj=utm +zone=8  +datum=NAD27 +units=m")
                       )
head(utmcoor@coords)

###
### Counts of sea otters
###

count.2017=all.2017$COUNT_ADULT+all.2017$COUNT_PUP

###
### Identify cell id
###

cell[]=1:q
cell.ID=extract(cell,utmcoor@coords)

###
### build matrix with coords, cell id, and count values
###

mat=cbind(utmcoor@coords,cell.ID,count.2017)
head(mat)
mat=mat[order(mat[,3]),]

###
### Create Y matrix
###

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
### Locations of observations
###

## obs.r=read.csv(paste("~/Dropbox/Post-Doc/",
##                      "OriginalDataFiles/SO_I_2017/SO_I_20170719_R.csv",
##                       sep=""))

## obs.o=read.csv(paste("~/Dropbox/Post-Doc/",
##                      "OriginalDataFiles/SO_I_2017/SO_I_20170721_O.csv",
##                       sep=""))

## obs.a=read.csv(paste("~/Dropbox/Post-Doc/",
##                      "OriginalDataFiles/SO_I_2017/SO_I_20170728_A.csv",
##                       sep=""))

## x_loc=c(obs.r$LONGITUDE_WGS84,
##         obs.o$LONGITUDE_WGS84,
##         obs.a$LONGITUDE_WGS84)
## y_loc=c(obs.r$LATITUDE_WGS84,
##         obs.o$LATITUDE_WGS84,
##         obs.a$LATITUDE_WGS84)
## Counts.tmp=c(obs.r$COUNT_ADULT+obs.r$COUNT_PUP,
##              obs.o$COUNT_ADULT+obs.o$COUNT_PUP,
##              obs.a$COUNT_ADULT+obs.a$COUNT_PUP
##              )

## xyc=cbind(x_loc,y_loc)
## xyc=xyc[!is.na(x_loc),]
## DDcoor=SpatialPoints(xyc,CRS(
##                            "+proj=longlat  +datum=WGS84")
##                        )

## utmcoor=spTransform(DDcoor,
##                        CRS(
##                            "+proj=utm +zone=8  +datum=NAD27 +units=m")
##                        )

## Counts.sub=rasterize(x=utmcoor@coords,
##                      y=data.r,
##                      field=Counts.tmp,fun="max",na.rm=TRUE,
##                      background=0)

## plot(Counts.sub)
## Counts[[21]][]=rep(NA,q)
## Counts[[22]][]=rep(NA,q)
## Counts[[23]][]=rep(NA,q)
## Counts[[24]][]=rep(NA,q)
## Counts[[25]][]=Counts.sub[]



###################################################################
### Transects
###################################################################

OrigData=read.csv("~/Dropbox/Post-Doc/OriginalDataFiles/GLBA_noDots.csv",
                  header=TRUE)
year=OrigData$year

###
### 1999
###

d.99.tmp=subset(OrigData,year==1999)

###
### Unique transects flown in 1999
###

ID=cumsum(!duplicated(d.99.tmp[,20:23]))
d.99=cbind(d.99.tmp,ID)

###
### Subset observation data for each transect
###

l.99=list()
Sl.99=list()
S.99=list()
## 324, 342
for(i in unique(ID)){
    t.sub=subset(d.99,ID==i)

    ##
    ## Store start, end, group location coordinates
    ##

    coords.tmp1=matrix(,dim(t.sub)*3,2)
    colnames(coords.tmp1)=c("x","y")

    ## West end point
    coords.tmp1[1:dim(t.sub)[1],1]=t.sub$west_long
    coords.tmp1[1:dim(t.sub)[1],2]=t.sub$west_lat

    ## Group locations
    coords.tmp1[(dim(t.sub)[1]+1):(2*dim(t.sub)[1]),1]=t.sub$GROUP_X
    coords.tmp1[(dim(t.sub)[1]+1):(2*dim(t.sub)[1]),2]=t.sub$GROUP_Y

    ## East end point
    coords.tmp1[(2*dim(t.sub)[1]+1):(3*dim(t.sub)[1]),1]=t.sub$east_long
    coords.tmp1[(2*dim(t.sub)[1]+1):(3*dim(t.sub)[1]),2]=t.sub$east_lat

    ## Remove duplicates
    coords.tmp2=unique(coords.tmp1)

    ## Remove NA rows
    coords.tmp3=coords.tmp2[!is.na(coords.tmp2[,1]),]

    ## Order from left to right
    coords=coords.tmp3[order(coords.tmp3[,1]),]

    l.99[[i]]=coords
    Sl.99[[i]]=Line(l.99[[i]])
    S.99[[i]]=Lines(list(Sl.99[[i]]), ID=i)
    Sb.99=SpatialLines(S.99)
    ## plot(Sb.99)
    ## readline()

}
## plot(Sb.99)


###
### 2000
###

d.00.tmp=subset(OrigData,year==2000)

###
### Unique transects flown in 2000
###

ID=cumsum(!duplicated(d.00.tmp[,20:23]))
d.00=cbind(d.00.tmp,ID)

###
### Subset observation data for each transect
###

l.00=list()
Sl.00=list()
S.00=list()

for(i in unique(ID)){ #493
    t.sub=subset(d.00,ID==i)

    ##
    ## Store start, end, group location coordinates
    ##

    coords.tmp1=matrix(,dim(t.sub)*3,2)
    colnames(coords.tmp1)=c("x","y")

    ## West end point
    coords.tmp1[1:dim(t.sub)[1],1]=t.sub$west_long
    coords.tmp1[1:dim(t.sub)[1],2]=t.sub$west_lat

    ## Group locations
    coords.tmp1[(dim(t.sub)[1]+1):(2*dim(t.sub)[1]),1]=t.sub$GROUP_X
    coords.tmp1[(dim(t.sub)[1]+1):(2*dim(t.sub)[1]),2]=t.sub$GROUP_Y

    ## East end point
    coords.tmp1[(2*dim(t.sub)[1]+1):(3*dim(t.sub)[1]),1]=t.sub$east_long
    coords.tmp1[(2*dim(t.sub)[1]+1):(3*dim(t.sub)[1]),2]=t.sub$east_lat

    ## Remove duplicates
    coords.tmp2=unique(coords.tmp1)

    ## Remove NA rows
    coords.tmp3=coords.tmp2[!is.na(coords.tmp2[,1]),]

    ## Order from left to right
    coords=coords.tmp3[order(coords.tmp3[,1]),]

    l.00[[i]]=coords
    Sl.00[[i]]=Line(l.00[[i]])
    S.00[[i]]=Lines(list(Sl.00[[i]]), ID=i)
    Sb.00=SpatialLines(S.00)
    ## plot(Sb.00)
    ## readline()

}

## plot(Sb.00)

###
### 2001
###

d.01.tmp=subset(OrigData,year==2001)

###
### Unique transects flown in 2001
###

ID=cumsum(!duplicated(d.01.tmp[,20:23]))
d.01=cbind(d.01.tmp,ID)

###
### Subset observation data for each transect
###

l.01=list()
Sl.01=list()
S.01=list()

for(i in unique(ID)[1:343]){  ## 343
    t.sub=subset(d.01,ID==i)

    ##
    ## Store start, end, group location coordinates
    ##

    coords.tmp1=matrix(,dim(t.sub)*3,2)
    colnames(coords.tmp1)=c("x","y")

    ## West end point
    coords.tmp1[1:dim(t.sub)[1],1]=t.sub$west_long
    coords.tmp1[1:dim(t.sub)[1],2]=t.sub$west_lat

    ## Group locations
    coords.tmp1[(dim(t.sub)[1]+1):(2*dim(t.sub)[1]),1]=t.sub$GROUP_X
    coords.tmp1[(dim(t.sub)[1]+1):(2*dim(t.sub)[1]),2]=t.sub$GROUP_Y

    ## East end point
    coords.tmp1[(2*dim(t.sub)[1]+1):(3*dim(t.sub)[1]),1]=t.sub$east_long
    coords.tmp1[(2*dim(t.sub)[1]+1):(3*dim(t.sub)[1]),2]=t.sub$east_lat

    ## Remove duplicates
    coords.tmp2=unique(coords.tmp1)

    ## Remove NA rows
    coords.tmp3=coords.tmp2[!is.na(coords.tmp2[,1]),]

    ## Order from left to right
    coords=coords.tmp3[order(coords.tmp3[,1]),]

    l.01[[i]]=coords
    Sl.01[[i]]=Line(l.01[[i]])
    S.01[[i]]=Lines(list(Sl.01[[i]]), ID=i)
    Sb.01=SpatialLines(S.01)
    ## plot(Sb.01)
    ## readline()

}
## plot(Sb.01)

###
### 2002
###

d.02.tmp=subset(OrigData,year==2002)

###
### Unique transects flown in 2002
###

ID=cumsum(!duplicated(d.02.tmp[,20:23]))
d.02=cbind(d.02.tmp,ID)

###
### Subset observation data for each transect
###

l.02=list()
Sl.02=list()
S.02=list()

for(i in unique(ID)){  ## 514
    t.sub=subset(d.02,ID==i)

    ##
    ## Store start, end, group location coordinates
    ##

    coords.tmp1=matrix(,dim(t.sub)*3,2)
    colnames(coords.tmp1)=c("x","y")

    ## West end point
    coords.tmp1[1:dim(t.sub)[1],1]=t.sub$west_long
    coords.tmp1[1:dim(t.sub)[1],2]=t.sub$west_lat

    ## Group locations
    coords.tmp1[(dim(t.sub)[1]+1):(2*dim(t.sub)[1]),1]=t.sub$GROUP_X
    coords.tmp1[(dim(t.sub)[1]+1):(2*dim(t.sub)[1]),2]=t.sub$GROUP_Y

    ## East end point
    coords.tmp1[(2*dim(t.sub)[1]+1):(3*dim(t.sub)[1]),1]=t.sub$east_long
    coords.tmp1[(2*dim(t.sub)[1]+1):(3*dim(t.sub)[1]),2]=t.sub$east_lat

    ## Remove duplicates
    coords.tmp2=unique(coords.tmp1)

    ## Remove NA rows
    coords.tmp3=coords.tmp2[!is.na(coords.tmp2[,1]),]

    ## Order from left to right
    coords=coords.tmp3[order(coords.tmp3[,1]),]

    l.02[[i]]=coords
    Sl.02[[i]]=Line(l.02[[i]])
    S.02[[i]]=Lines(list(Sl.02[[i]]), ID=i)
    Sb.02=SpatialLines(S.02)
    ## plot(Sb.02)
    ## readline()

}
## plot(Sb.02)

###
### 2003
###

d.03.tmp=subset(OrigData,year==2003)

###
### Unique transects flown in 2003
###

ID=cumsum(!duplicated(d.03.tmp[,20:23]))
d.03=cbind(d.03.tmp,ID)

###
### Subset observation data for each transect
###

l.03=list()
Sl.03=list()
S.03=list()

for(i in unique(ID)){  ## 514
    t.sub=subset(d.03,ID==i)

    ##
    ## Store start, end, group location coordinates
    ##

    coords.tmp1=matrix(,dim(t.sub)*3,2)
    colnames(coords.tmp1)=c("x","y")

    ## West end point
    coords.tmp1[1:dim(t.sub)[1],1]=t.sub$west_long
    coords.tmp1[1:dim(t.sub)[1],2]=t.sub$west_lat

    ## Group locations
    coords.tmp1[(dim(t.sub)[1]+1):(2*dim(t.sub)[1]),1]=t.sub$GROUP_X
    coords.tmp1[(dim(t.sub)[1]+1):(2*dim(t.sub)[1]),2]=t.sub$GROUP_Y

    ## East end point
    coords.tmp1[(2*dim(t.sub)[1]+1):(3*dim(t.sub)[1]),1]=t.sub$east_long
    coords.tmp1[(2*dim(t.sub)[1]+1):(3*dim(t.sub)[1]),2]=t.sub$east_lat

    ## Remove duplicates
    coords.tmp2=unique(coords.tmp1)

    ## Remove NA rows
    coords.tmp3=coords.tmp2[!is.na(coords.tmp2[,1]),]

    ## Order from left to right
    coords=coords.tmp3[order(coords.tmp3[,1]),]

    l.03[[i]]=coords
    Sl.03[[i]]=Line(l.03[[i]])
    S.03[[i]]=Lines(list(Sl.03[[i]]), ID=i)
    Sb.03=SpatialLines(S.03)
    ## plot(Sb.03)
    ## readline()

}
## plot(Sb.03)

###
### 2004
###

d.04.tmp=subset(OrigData,year==2004)

###
### Unique transects flown in 2004
###

ID=cumsum(!duplicated(d.04.tmp[,20:23]))
d.04=cbind(d.04.tmp,ID)

###
### Subset observation data for each transect
###

l.04=list()
Sl.04=list()
S.04=list()

for(i in unique(ID)){  ## 417
    t.sub=subset(d.04,ID==i)

    ##
    ## Store start, end, group location coordinates
    ##

    coords.tmp1=matrix(,dim(t.sub)*3,2)
    colnames(coords.tmp1)=c("x","y")

    ## West end point
    coords.tmp1[1:dim(t.sub)[1],1]=t.sub$west_long
    coords.tmp1[1:dim(t.sub)[1],2]=t.sub$west_lat

    ## Group locations
    coords.tmp1[(dim(t.sub)[1]+1):(2*dim(t.sub)[1]),1]=t.sub$GROUP_X
    coords.tmp1[(dim(t.sub)[1]+1):(2*dim(t.sub)[1]),2]=t.sub$GROUP_Y

    ## East end point
    coords.tmp1[(2*dim(t.sub)[1]+1):(3*dim(t.sub)[1]),1]=t.sub$east_long
    coords.tmp1[(2*dim(t.sub)[1]+1):(3*dim(t.sub)[1]),2]=t.sub$east_lat

    ## Remove duplicates
    coords.tmp2=unique(coords.tmp1)

    ## Remove NA rows
    coords.tmp3=coords.tmp2[!is.na(coords.tmp2[,1]),]

    ## Order from left to right
    coords=coords.tmp3[order(coords.tmp3[,1]),]

    l.04[[i]]=coords
    Sl.04[[i]]=Line(l.04[[i]])
    S.04[[i]]=Lines(list(Sl.04[[i]]), ID=i)
    Sb.04=SpatialLines(S.04)

}
## plot(Sb.04)

###
### 2006
###

d.06.tmp=subset(OrigData,year==2006)

###
### Unique transects flown in 2006
###

ID=cumsum(!duplicated(d.06.tmp[,20:23]))
d.06=cbind(d.06.tmp,ID)

###
### Subset observation data for each transect
###

l.06=list()
Sl.06=list()
S.06=list()

for(i in unique(ID)){  ## 365
    t.sub=subset(d.06,ID==i)

    ##
    ## Store start, end, group location coordinates
    ##

    coords.tmp1=matrix(,dim(t.sub)*3,2)
    colnames(coords.tmp1)=c("x","y")

    ## West end point
    coords.tmp1[1:dim(t.sub)[1],1]=t.sub$west_long
    coords.tmp1[1:dim(t.sub)[1],2]=t.sub$west_lat

    ## Group locations
    coords.tmp1[(dim(t.sub)[1]+1):(2*dim(t.sub)[1]),1]=t.sub$GROUP_X
    coords.tmp1[(dim(t.sub)[1]+1):(2*dim(t.sub)[1]),2]=t.sub$GROUP_Y

    ## East end point
    coords.tmp1[(2*dim(t.sub)[1]+1):(3*dim(t.sub)[1]),1]=t.sub$east_long
    coords.tmp1[(2*dim(t.sub)[1]+1):(3*dim(t.sub)[1]),2]=t.sub$east_lat

    ## Remove duplicates
    coords.tmp2=unique(coords.tmp1)

    ## Remove NA rows
    coords.tmp3=coords.tmp2[!is.na(coords.tmp2[,1]),]

    ## Order from left to right
    coords=coords.tmp3[order(coords.tmp3[,1]),]

    l.06[[i]]=coords
    Sl.06[[i]]=Line(l.06[[i]])
    S.06[[i]]=Lines(list(Sl.06[[i]]), ID=i)
    Sb.06=SpatialLines(S.06)

}
## plot(Sb.06)

###
### 2012
###

d.12.tmp=subset(OrigData,year==2012)

###
### Unique transects flown in 2012
###

ID=cumsum(!duplicated(d.12.tmp[,20:23]))
d.12=cbind(d.12.tmp,ID)

###
### Subset observation data for each transect
###

l.12=list()
Sl.12=list()
S.12=list()

for(i in unique(ID)[1:324]){  ## 347
    t.sub=subset(d.12,ID==i)

    ##
    ## Store start, end, group location coordinates
    ##

    coords.tmp1=matrix(,dim(t.sub)*3,2)
    colnames(coords.tmp1)=c("x","y")

    ## West end point
    coords.tmp1[1:dim(t.sub)[1],1]=t.sub$west_long
    coords.tmp1[1:dim(t.sub)[1],2]=t.sub$west_lat

    ## Group locations
    coords.tmp1[(dim(t.sub)[1]+1):(2*dim(t.sub)[1]),1]=t.sub$GROUP_X
    coords.tmp1[(dim(t.sub)[1]+1):(2*dim(t.sub)[1]),2]=t.sub$GROUP_Y

    ## East end point
    coords.tmp1[(2*dim(t.sub)[1]+1):(3*dim(t.sub)[1]),1]=t.sub$east_long
    coords.tmp1[(2*dim(t.sub)[1]+1):(3*dim(t.sub)[1]),2]=t.sub$east_lat

    ## Remove duplicates
    coords.tmp2=unique(coords.tmp1)

    ## Remove NA rows
    coords.tmp3=coords.tmp2[!is.na(coords.tmp2[,1]),]

    ## Order from left to right
    coords=coords.tmp3[order(coords.tmp3[,1]),]

    l.12[[i]]=coords
    Sl.12[[i]]=Line(l.12[[i]])
    S.12[[i]]=Lines(list(Sl.12[[i]]), ID=i)
    Sb.12=SpatialLines(S.12)

}
## plot(Sb.12)

## ###
## ### 2017
## ###

## d.17.tmp=read.csv(paste("~/Dropbox/Post-Doc/",
##                         "SurveyRoute_2017/2017AllTransects.csv",
##                         sep=""))


## DDcoor=SpatialPoints(d.17.tmp[,7:6],CRS(
##                            "+proj=longlat  +datum=WGS84")
##                        )

## utmcoor=spTransform(DDcoor,
##                        CRS(
##                            "+proj=utm +zone=8  +datum=NAD27 +units=m")
##                        )

## d.17.tmp[,7:6]=utmcoor@coords

## ###
## ### Unique transects flown in 2017
## ###

## ID=cumsum(!duplicated(d.17.tmp[,1:2]))
## d.17=cbind(d.17.tmp,ID)


## ###
## ### Subset observation data for each transect
## ###

## l.17=list()
## Sl.17=list()
## S.17=list()

## for(i in unique(ID)){  ## 347
##     t.sub=subset(d.17,ID==i)

##     ##
##     ## Store start, end, group location coordinates
##     ##

##     coords.tmp1=matrix(,dim(t.sub)[1]*3,2)
##     colnames(coords.tmp1)=c("x","y")

##     ## West end point
##     coords.tmp1[1:dim(t.sub)[1],1]=t.sub[t.sub$Side=="W",]$Longitude
##     coords.tmp1[1:dim(t.sub)[1],2]=t.sub[t.sub$Side=="W",]$Latitude

##     ## ## Group locations
##     ## coords.tmp1[(dim(t.sub)[1]+1):(2*dim(t.sub)[1]),1]=t.sub$GROUP_X
##     ## coords.tmp1[(dim(t.sub)[1]+1):(2*dim(t.sub)[1]),2]=t.sub$GROUP_Y

##     ## East end point
##     coords.tmp1[(2*dim(t.sub)[1]+1):(3*dim(t.sub)[1]),1]=t.sub[t.sub$Side=="E",]$Longitude
##     coords.tmp1[(2*dim(t.sub)[1]+1):(3*dim(t.sub)[1]),2]=t.sub[t.sub$Side=="E",]$Latitude

##     ## Remove duplicates
##     coords.tmp2=unique(coords.tmp1)

##     ## Remove NA rows
##     coords.tmp3=coords.tmp2[!is.na(coords.tmp2[,1]),]

##     ## Order from left to right
##     coords=coords.tmp3[order(coords.tmp3[,1]),]

##     l.17[[i]]=coords
##     Sl.17[[i]]=Line(l.17[[i]])
##     S.17[[i]]=Lines(list(Sl.17[[i]]), ID=i)
##     Sb.17=SpatialLines(S.17)

## }
## plot(Sb.17)

##################################################
### Rasterize transects
##################################################

transectAll=raster(,nrows=y,ncols=x,xmn=xmin,
                  xmx=xmax,ymn=ymin,ymx=ymax,crs=NA)
transectAll[]=rep(NA,q)
transectAll=stack(mget(rep("transectAll",25)))
transectAll[[1]][]=rep(NA,q)
transectAll[[2]][]=rep(NA,q)
transectAll[[3]][]=rep(NA,q)
transectAll[[4]][]=rep(NA,q)
transectAll[[5]][]=rep(NA,q)
transectAll[[6]][]=rep(NA,q)
transectAll[[7]]=rasterize(x=Sb.99,y=transectAll[[7]],
                         field=1,fun="last",background=NA)
transectAll[[8]]=rasterize(x=Sb.00,y=transectAll[[8]],
                         field=1,fun="last",background=NA)
transectAll[[9]]=rasterize(x=Sb.01,y=transectAll[[9]],
                         field=1,fun="last",background=NA)
transectAll[[10]]=rasterize(x=Sb.02,y=transectAll[[10]],
                         field=1,fun="last",background=NA)
transectAll[[11]]=rasterize(x=Sb.03,y=transectAll[[11]],
                         field=1,fun="last",background=NA)
transectAll[[12]]=rasterize(x=Sb.04,y=transectAll[[12]],
                         field=1,fun="last",background=NA)
transectAll[[13]][]=rep(NA,q)
transectAll[[14]]=rasterize(x=Sb.06,y=transectAll[[14]],
                         field=1,fun="last",background=NA)
transectAll[[15]][]=rep(NA,q)
transectAll[[16]][]=rep(NA,q)
transectAll[[17]][]=rep(NA,q)
transectAll[[18]][]=rep(NA,q)
transectAll[[19]][]=rep(NA,q)
transectAll[[20]]=rasterize(x=Sb.12,y=transectAll[[20]],
                         field=1,fun="last",background=NA)
## transectAll[[21]][]=rep(NA,q)
## transectAll[[22]][]=rep(NA,q)
## transectAll[[23]][]=rep(NA,q)
## transectAll[[24]][]=rep(NA,q)
## transectAll[[25]]=rasterize(x=Sb.17,y=transectAll[[25]],
##                          field=1,fun="last",background=NA)

## plot(Counts[[25]])
## plot(transectAll[[25]],add=TRUE)

#####################################################
### Combine transect data and count data
#####################################################

Y.r[]=Counts[[7]][]*transectAll[[7]][]
Y.r=stack(mget(rep("Y.r",20)))
for (i in 1:20) {
  Y.r[[i]][]=ifelse(Counts[[i]][]>0,Counts[[i]][],
    transectAll[[i]][]-1)
}
## plot(Y.r[[20]])
## plot(Y.r[[25]])

######################################################
### Load ISU data
######################################################

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

###
### Create `Cell' raster
###

cell[]=1:q
## plot(cell)

###
### Create Bounday raster
###

Boundary[is.na(bath)]=0
Boundary[bath<0]=1
Boundary[bath>=0]=0
## plot(Boundary)

BoundaryNA=Boundary
BoundaryNA[Boundary==0]=NA
## plot(BoundaryNA)

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
## plot(Boundary.us)
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
DepthCov[is.na(bath[])]=0
DepthCov=DepthCov*Boundary
## plot(DepthCov)

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

###
### Bottom slope
###

SlopeCov.tmp=terrain(bath,opt='slope',unit='degrees',neighbors=8)
SlopeCov.tmp[is.na(SlopeCov.tmp)]=0
SlopeCov=scale(SlopeCov.tmp, center=TRUE, scale=TRUE)
SlopeCov=SlopeCov*Boundary*DepthCov
## plot(SlopeCov)

###
### Shoreline complexity
###

ShoreCov.tmp=boundaries(DistCov.tmp1,type='inner')
ShoreCov.tmp[is.na(ShoreCov.tmp)]=0
ShoreCov=focal(ShoreCov.tmp,w=matrix(1,nr=11,nc=11))*Boundary
ShoreCov[is.na(ShoreCov)]=0
ShoreCov=scale(ShoreCov,center=TRUE, scale=TRUE)
ShoreCov=ShoreCov*Boundary
## plot(ShoreCov)
## points(loc.data.UTM$lon,loc.data.UTM$lat,cex=0.2,pch=16)

###
### Create co-variate matrices
###

X=cbind(1,DepthCov[],DistCov[],SlopeCov[]*DepthCov[],ShoreCov[])
## cor(X, use="pairwise.complete.obs")
W=matrix(1,nr=q,nc=1)


###
### Data for analysis
###

Y=c(Y.r[[1]][],
    Y.r[[2]][],
    Y.r[[3]][],
    Y.r[[4]][],
    Y.r[[5]][],
    Y.r[[6]][],
    Y.r[[7]][],
    Y.r[[8]][],
    Y.r[[9]][],
    Y.r[[10]][],
    Y.r[[11]][],
    Y.r[[12]][],
    Y.r[[13]][],
    Y.r[[14]][],
    Y.r[[15]][],
    Y.r[[16]][],
    Y.r[[17]][],
    Y.r[[18]][],
    Y.r[[19]][],
    Y.r[[20]][]## ,
    ## Y.r[[21]][],
    ## Y.r[[22]][],
    ## Y.r[[23]][],
    ## Y.r[[24]][],
    ## Y.r[[25]][]
    )

#################################################
#################################################
###
### Fit the model to the real data using MCMC
###
#################################################
#################################################


###
### Bundle data
###

data=list(Y=Y,
          Y.2017=Y.2017,
          ISU=ISU,
          X=X)

###
### Spatio-temporal settings
###

dt=1/200
time.frame=1993:2018
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
checkpoint=1000

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
p.2017=0.05
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
    rep(p.2012,q),
    rep(NA,q*4),
    rep(p.2017,q),
    rep(NA,q)
    )
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

output.location="~/GB_2018_DataAnalysis.RData"

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

## script=getURL(
##     paste("https://raw.githubusercontent.com/perrywilliams/",
##           "GlacierBay-SeaOtter-Analysis/master/MCMC.R",sep=""),
##     ssl.verifypeer=FALSE)
## eval(parse(text = script))

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
output.location="~/GB_2018_DataAnalysis.RData"
load(output.location)
(status=sum(!is.na(MCMC.Chains[[3]])))

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

par(mfrow=c(3,2),mar=c(2,2,0,1))
for(i in 21:25){
    plot(N.est[,i],type='l')
}

###
### Plot p
###

par(mfrow=c(3,3))
p.est=MCMC.Chains$p[1:status,]
for(i in 1:9){
    plot(p.est[,i],type='l')
}
