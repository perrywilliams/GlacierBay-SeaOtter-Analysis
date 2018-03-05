############################################################################
############################################################################
############################################################################
###
### Real Data Analysis
###
############################################################################
############################################################################
############################################################################

rm(list=ls())

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
## plot(bath)
## points(loc.data.UTM$lon,loc.data.UTM$lat,cex=0.2,pch=16)

###
### Load Glacier Bay Sea Otter Data
###

data=read.csv("~/Dropbox/Post-Doc/OriginalDataFiles/GLBA_noDots.csv")
dataDistr=read.csv("~/Dropbox/Post-Doc/OriginalDataFiles/DistribSurveysC.csv")
ind.pred.tmp=1:dim(dataDistr)[1]
ind.pred=ind.pred.tmp[dataDistr$year==2004|dataDistr$year==2006]
dataDistr=dataDistr[-ind.pred,]
ind=1993:2012
counts=numeric(20)
for(i in 1:length(ind)){
    dat=subset(dataDistr,dataDistr$year==ind[i])
    counts[i]=sum(dat$animals,na.rm=TRUE)
}
counts

###
### Convert Data to Raster
###

cell=raster(,nrows=y,ncols=x,xmn=xmin,xmx=xmax,
            ymn=ymin,ymx=ymax,crs=NA)
CISU=N.r=ISUind.r=Y.r=Counts=Boundary=BoundaryDist=DistCov=DepthCov=gamma=delta=lambda0=SurvR=SimulatedDataR=cell
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
Counts=stack(mget(rep("Counts",20)))
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
## plot(Counts)

###################################################################
### Transects
###################################################################

OrigData=read.csv("~/Dropbox/Post-Doc/OriginalDataFiles/GLBA_noDots.csv",header=TRUE)
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

##################################################
### Rasterize transects
##################################################

transectAll=raster(,nrows=y,ncols=x,xmn=xmin,
                  xmx=xmax,ymn=ymin,ymx=ymax,crs=NA)
transectAll[]=rep(NA,q)
transectAll=stack(mget(rep("transectAll",20)))
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

## plot(transectAll)

#####################################################
### Combine transect data and count data
#####################################################

Y.r[]=Counts[[7]][]*transectAll[[7]][]
Y.r=stack(mget(rep("Y.r",20)))
for (i in 1:20) {
  Y.r[[i]][]=ifelse(Counts[[i]][]>0,Counts[[i]][],
    transectAll[[i]][]-1)
}
## plot(Y.r[[7]])

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
    Y.r[[20]][])

#################################################
#################################################
###
### Fit the model to the real data using MCMC
###
#################################################
#################################################

###
### MCMC settings
###

n.mcmc=200000
checkpoint=100

###
### Priors
###

## Gamma
q.gamma=0.01
r.gamma=0.01

## Beta
mu.beta=0
sigma2.beta=10^2

## Theta
mu.theta=500
sigma2.theta=500

## kappa
mu.kappa=5
sigma2.kappa=5

## detection probability
q.p=1
r.p=1

## Overdispersion parameter
l.odp.mn=log(2)
l.odp.sd=log(3)

###
### Temporal settings
###

dt=1/200
time.frame=1993:2012
time.steps=1/dt*length(time.frame)
keep=seq(1,time.steps,time.steps/length(1993:2012))

###
### Spatial settings
###

us.fact=10
dx=400*us.fact
dy=400*us.fact

###
### Tuning parameters
###

gamma.tune=0.002
beta.tune=c(0.1392345, 0.04616871, 0.1127799, 0.05642842, 0.03912024)
theta.tune=40.81335
kappa.tune=0.9653723
odp.tune=0.02357948

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
odp=2
delta[]=exp(X%*%beta)
delta.star=delta
delta.bar=aggregate(delta,fact=us.fact,
                        fun=function(x,na.rm){
                            (1/mean(1/x,na.rm=TRUE))
                        })
gamma.bar=delta.bar*aggregate(gamma/delta,
                              fact=us.fact,
                              fun=mean,
                              na.rm=TRUE
                              )

NN=neighborhood(delta.bar,Boundary.us)

H=propagator.plainZF(NN=NN,
                     delta=delta.bar[],
                     gamma=gamma.bar[],
                     dx=dx,dy=dy,dt=dt)

accept.gamma=0
accept.beta=rep(0,length(beta))
accept.theta=0
accept.kappa=0
accept.N=rep(0,length(Y))
accept.odp=0

###
### Initial condition
###

D=rdist(data.frame(SpatialPoints(lambda0)),matrix(d,1,2))/1000
us.cells=aggregate(cell,fact=us.fact,fun=mean)
us.cells[]=1:length(us.cells[])
lambda0[]=(exp(-D^2/kappa^2)/sum(exp(-D^2/kappa^2))*theta)
lambda0.star=lambda0
c0[]=extract(delta*lambda0,SpatialPoints(us.cells))
c0.star=c0

###
### c.all and lambda.all and N.all
###

c.all=brick(nrows=dim(c0)[1], ncols=dim(c0)[2], xmn=xmin, xmx=xmax,
            ymn=ymin, ymx=ymax)
c.all=setValues(c.all, calcc(H, vec(c0[]),time.steps)[,keep])
lambda=vec(disaggregate(c.all,us.fact)/delta)
N=rnbinom(length(lambda),size=odp,,mu=lambda)
N=ifelse(N<Y,Y,N)

###
### Results objects
###

gamma.save=rep(NA,n.mcmc)
gamma.tune.save=rep(NA,n.mcmc)
beta.save=matrix(NA,n.mcmc,length(beta))
beta.tune.save=matrix(NA,n.mcmc,length(beta))
theta.save=rep(NA,n.mcmc)
theta.tune.save=rep(NA,n.mcmc)
kappa.save=rep(NA,n.mcmc)
kappa.tune.save=rep(NA,n.mcmc)
p.save=matrix(NA,n.mcmc,8)
odp.save=rep(NA,n.mcmc)
odp.tune.save=rep(NA,n.mcmc)
n.tot.save=matrix(NA,n.mcmc,20)

###
### Begin MCMC loop
###

for(k in 1:n.mcmc){

    ##
    ## Sample gamma
    ##

    gamma.star=rnorm(n=1,mean=gamma,sd=gamma.tune)
    if(gamma.star>0){
        gamma.bar.star=delta.bar*aggregate(gamma.star/delta,
                                           fact=us.fact,fun=function(x,na.rm){
                                               mean(x,na.rm=TRUE)
                                           })
        H.star=propagator.plainZF(NN,delta.bar[],gamma.bar.star[],
                                  dx=dx,dy=dy,dt=dt)
        if(min(range(H.star,na.rm=TRUE))>=0){
            c.all.star=setValues(c.all,calcc(H.star,vec(c0[]),time.steps)[,keep])
            lambda.star=vec(disaggregate(c.all.star,us.fact)/delta)
            mh1=sum(dnbinom(x=N,size=odp,,mu=lambda.star,log=TRUE),na.rm=TRUE)+
                sum(dbeta(gamma.star,q.gamma,r.gamma,log=TRUE))
            mh2=sum(dnbinom(x=N,size=odp,,mu=lambda,log=TRUE),na.rm=TRUE)+
                sum(dbeta(gamma,q.gamma,r.gamma,log=TRUE))
            mh=exp(mh1-mh2)
            if(mh>runif(1)){
                gamma=gamma.star
                gamma.bar=gamma.bar.star
                H=H.star
                c.all=c.all.star
                lambda=lambda.star
                accept.gamma=accept.gamma+1
            }
        }
    }

    ##
    ## Sample beta0
    ##

    beta0.star=rnorm(1,mean=beta[1],sd=beta.tune[1])
    beta.star=c(beta0.star,beta[2:5])
    delta.star[]=exp(X%*%beta.star)
    delta.bar.star=aggregate(delta.star,
                             fact=us.fact,
                             fun=function(x,na.rm){
                                 (1/mean(1/x,na.rm=TRUE))
                             }
                             )
    gamma.bar.star=delta.bar.star*aggregate(gamma/delta.star,
                                            fact=us.fact,fun=function(x,na.rm){
                                                mean(x,na.rm=TRUE)
                                            })
    H.star=propagator.plainZF(NN,delta.bar.star[],gamma.bar.star[],
                                  dx=dx,dy=dy,dt=dt)
    if(min(H.star,na.rm=TRUE)>=0){
        c0.star[]=extract(delta.star*lambda0,
                        SpatialPoints(us.cells))
        c.all.star=setValues(c.all,calcc(H.star,vec(c0.star[]),time.steps)[,keep])
        lambda.star=vec(disaggregate(c.all.star,us.fact)/delta.star)
        mh1=sum(dnbinom(x=N,size=odp,,mu=lambda.star,log=TRUE),na.rm=TRUE)+
            sum(dnorm(beta.star,mu.beta,sigma2.beta^0.5,log=TRUE))
        mh2=sum(dnbinom(x=N,size=odp,,mu=lambda,log=TRUE),na.rm=TRUE)+
            sum(dnorm(beta,mu.beta,sigma2.beta^0.5,log=TRUE))
        mh=exp(mh1-mh2)
        if(mh>runif(1)){
            beta=beta.star
            delta=delta.star
            delta.bar=delta.bar.star
            gamma.bar=gamma.bar.star
            H=H.star
            c0=c0.star
            c.all=c.all.star
            lambda=lambda.star
            accept.beta[1]=accept.beta[1]+1
        }
    }

    ##
    ## Sample beta1
    ##

    beta1.star=rnorm(1,mean=beta[2],sd=beta.tune[2])
    beta.star=c(beta[1],beta1.star,beta[3:5])
    delta.star[]=exp(X%*%beta.star)
    delta.bar.star=aggregate(delta.star,
                             fact=us.fact,
                             fun=function(x,na.rm){
                                 (1/mean(1/x,na.rm=TRUE))
                             }
                             )
    gamma.bar.star=delta.bar.star*aggregate(gamma/delta.star,
                                            fact=us.fact,fun=function(x,na.rm){
                                                mean(x,na.rm=TRUE)
                                            })
    H.star=propagator.plainZF(NN,delta.bar.star[],gamma.bar.star[],
                                  dx=dx,dy=dy,dt=dt)
    if(min(H.star,na.rm=TRUE)>=0){
        c0.star[]=extract(delta.star*lambda0,
                        SpatialPoints(us.cells))
        c.all.star=setValues(c.all,calcc(H.star,vec(c0.star[]),time.steps)[,keep])
        lambda.star=vec(disaggregate(c.all.star,us.fact)/delta.star)
        mh1=sum(dnbinom(x=N,size=odp,,mu=lambda.star,log=TRUE),na.rm=TRUE)+
            sum(dnorm(beta.star,mu.beta,sigma2.beta^0.5,log=TRUE))
        mh2=sum(dnbinom(x=N,size=odp,,mu=lambda,log=TRUE),na.rm=TRUE)+
            sum(dnorm(beta,mu.beta,sigma2.beta^0.5,log=TRUE))
        mh=exp(mh1-mh2)
        if(mh>runif(1)){
            beta=beta.star
            delta=delta.star
            delta.bar=delta.bar.star
            gamma.bar=gamma.bar.star
            H=H.star
            c0=c0.star
            c.all=c.all.star
            lambda=lambda.star
            accept.beta[2]=accept.beta[2]+1
        }
    }

    ##
    ## Sample beta2
    ##

    beta2.star=rnorm(1,mean=beta[3],sd=beta.tune[3])
    beta.star=c(beta[1:2],beta2.star,beta[4:5])
    delta.star[]=exp(X%*%beta.star)
    delta.bar.star=aggregate(delta.star,
                             fact=us.fact,
                             fun=function(x,na.rm){
                                 (1/mean(1/x,na.rm=TRUE))
                             }
                             )
    gamma.bar.star=delta.bar.star*aggregate(gamma/delta.star,
                                            fact=us.fact,fun=function(x,na.rm){
                                                mean(x,na.rm=TRUE)
                                            })
    H.star=propagator.plainZF(NN,delta.bar.star[],gamma.bar.star[],
                                  dx=dx,dy=dy,dt=dt)
    if(min(H.star,na.rm=TRUE)>=0){
        c0.star[]=extract(delta.star*lambda0,
                        SpatialPoints(us.cells))
        c.all.star=setValues(c.all,calcc(H.star,vec(c0.star[]),time.steps)[,keep])
        lambda.star=vec(disaggregate(c.all.star,us.fact)/delta.star)
        mh1=sum(dnbinom(x=N,size=odp,,mu=lambda.star,log=TRUE),na.rm=TRUE)+
            sum(dnorm(beta.star,mu.beta,sigma2.beta^0.5,log=TRUE))
        mh2=sum(dnbinom(x=N,size=odp,,mu=lambda,log=TRUE),na.rm=TRUE)+
            sum(dnorm(beta,mu.beta,sigma2.beta^0.5,log=TRUE))
        mh=exp(mh1-mh2)
        if(mh>runif(1)){
            beta=beta.star
            delta=delta.star
            delta.bar=delta.bar.star
            gamma.bar=gamma.bar.star
            H=H.star
            c0=c0.star
            c.all=c.all.star
            lambda=lambda.star
            accept.beta[3]=accept.beta[3]+1
        }
    }

    ##
    ## Sample beta3
    ##

    beta3.star=rnorm(1,mean=beta[4],sd=beta.tune[4])
    beta.star=c(beta[1:3],beta3.star,beta[5])
    delta.star[]=exp(X%*%beta.star)
    delta.bar.star=aggregate(delta.star,
                             fact=us.fact,
                             fun=function(x,na.rm){
                                 (1/mean(1/x,na.rm=TRUE))
                             }
                             )
    gamma.bar.star=delta.bar.star*aggregate(gamma/delta.star,
                                            fact=us.fact,fun=function(x,na.rm){
                                                mean(x,na.rm=TRUE)
                                            })
    H.star=propagator.plainZF(NN,delta.bar.star[],gamma.bar.star[],
                                  dx=dx,dy=dy,dt=dt)
    if(min(H.star,na.rm=TRUE)>=0){
        c0.star[]=extract(delta.star*lambda0,
                        SpatialPoints(us.cells))
        c.all.star=setValues(c.all,calcc(H.star,vec(c0.star[]),time.steps)[,keep])
        lambda.star=vec(disaggregate(c.all.star,us.fact)/delta.star)
        mh1=sum(dnbinom(x=N,size=odp,,mu=lambda.star,log=TRUE),na.rm=TRUE)+
            sum(dnorm(beta.star,mu.beta,sigma2.beta^0.5,log=TRUE))
        mh2=sum(dnbinom(x=N,size=odp,,mu=lambda,log=TRUE),na.rm=TRUE)+
            sum(dnorm(beta,mu.beta,sigma2.beta^0.5,log=TRUE))
        mh=exp(mh1-mh2)
        if(mh>runif(1)){
            beta=beta.star
            delta=delta.star
            delta.bar=delta.bar.star
            gamma.bar=gamma.bar.star
            H=H.star
            c0=c0.star
            c.all=c.all.star
            lambda=lambda.star
            accept.beta[4]=accept.beta[4]+1
        }
    }

    ##
    ## Sample beta4
    ##

    beta4.star=rnorm(1,mean=beta[5],sd=beta.tune[5])
    beta.star=c(beta[1:4],beta4.star)
    delta.star[]=exp(X%*%beta.star)
    delta.bar.star=aggregate(delta.star,
                             fact=us.fact,
                             fun=function(x,na.rm){
                                 (1/mean(1/x,na.rm=TRUE))
                             }
                             )
    gamma.bar.star=delta.bar.star*aggregate(gamma/delta.star,
                                            fact=us.fact,fun=function(x,na.rm){
                                                mean(x,na.rm=TRUE)
                                            })
    H.star=propagator.plainZF(NN,delta.bar.star[],gamma.bar.star[],
                                  dx=dx,dy=dy,dt=dt)
    if(min(H.star,na.rm=TRUE)>=0){
        c0.star[]=extract(delta.star*lambda0,
                        SpatialPoints(us.cells))
        c.all.star=setValues(c.all,calcc(H.star,vec(c0.star[]),time.steps)[,keep])
        lambda.star=vec(disaggregate(c.all.star,us.fact)/delta.star)
        mh1=sum(dnbinom(x=N,size=odp,,mu=lambda.star,log=TRUE),na.rm=TRUE)+
            sum(dnorm(beta.star,mu.beta,sigma2.beta^0.5,log=TRUE))
        mh2=sum(dnbinom(x=N,size=odp,,mu=lambda,log=TRUE),na.rm=TRUE)+
            sum(dnorm(beta,mu.beta,sigma2.beta^0.5,log=TRUE))
        mh=exp(mh1-mh2)
        if(mh>runif(1)){
            beta=beta.star
            delta=delta.star
            delta.bar=delta.bar.star
            gamma.bar=gamma.bar.star
            H=H.star
            c0=c0.star
            c.all=c.all.star
            lambda=lambda.star
            accept.beta[5]=accept.beta[5]+1
        }
    }

    ##
    ## Sample theta
    ##

    theta.star=rnorm(1,theta,theta.tune)
    if(theta.star>0){
        lambda0.star[]=exp(-D^2/kappa^2)/sum(exp(-D^2/kappa^2))*theta.star
        c0.star=extract(delta*lambda0.star,
                        SpatialPoints(us.cells))
        c.all.star=setValues(c.all,calcc(H,vec(c0.star[]),time.steps)[,keep])
        lambda.star=vec(disaggregate(c.all.star,us.fact)/delta)
        mh1=sum(dnbinom(x=N,size=odp,,mu=lambda.star,log=TRUE),na.rm=TRUE)+
            dnorm(theta.star,mu.theta,sigma2.theta^0.5,log=TRUE)
        mh2=sum(dnbinom(x=N,size=odp,,mu=lambda,log=TRUE),na.rm=TRUE)+
            dnorm(theta,mu.theta,sigma2.theta^0.5,log=TRUE)
        mh=exp(mh1-mh2)
        if(mh>runif(1)){
            theta=theta.star
            lambda0=lambda0.star
            c0=c0.star
            c.all=c.all.star
            lambda=lambda.star
            accept.theta=accept.theta+1
        }
    }

    ##
    ## Sample kappa
    ##

    kappa.star=rnorm(1,kappa,kappa.tune)
    if(kappa.star>0){
        lambda0.star[]=exp(-D^2/kappa.star^2)/sum(exp(-D^2/kappa.star^2))*theta
        c0.star=extract(delta*lambda0.star,
                        SpatialPoints(us.cells))
        c.all.star=setValues(c.all,calcc(H,vec(c0.star[]),time.steps)[,keep])
        lambda.star=vec(disaggregate(c.all.star,us.fact)/delta)
        mh1=sum(dnbinom(x=N,size=odp,,mu=lambda.star,log=TRUE),na.rm=TRUE)+
            dnorm(kappa.star,mu.kappa,sigma2.kappa^0.5,log=TRUE)
        mh2=sum(dnbinom(x=N,size=odp,,mu=lambda,log=TRUE),na.rm=TRUE)+
            dnorm(kappa,mu.kappa,sigma2.kappa^0.5,log=TRUE)
        mh=exp(mh1-mh2)
        if(mh>runif(1)){
            kappa=kappa.star
            lambda0=lambda0.star
            c0=c0.star
            c.all=c.all.star
            lambda=lambda.star
            accept.kappa=accept.kappa+1
        }
    }

    ##
    ## Sample N
    ##

    N.star=sample(-1:1,length(N),replace=TRUE)+N
    N.star=ifelse(N.star<0,0,N.star)
    mh1=dbinom(x=Y,size=N.star,prob=p,log=TRUE)+
        dnbinom(x=N.star,size=odp,,mu=lambda,log=TRUE)
    mh2=dbinom(x=Y,size=N,prob=p,log=TRUE)+
        dnbinom(x=N,size=odp,,mu=lambda,log=TRUE)
    mh=exp(mh1-mh2)
    ru.tmp=runif(length(mh))
    N=ifelse(mh>ru.tmp,
             N.star,
             N)
    accept.N=ifelse(mh>ru.tmp,
                    accept.N+1,
                    accept.N)

    ##
    ## Sample odp
    ##

    odp.star=exp(rnorm(1,log(odp),odp.tune))
    mh1=sum(dnbinom(x=N,size=odp.star,,mu=lambda,log=TRUE),na.rm=TRUE)+
        dnorm(log(odp.star),l.odp.mn,l.odp.sd,log=TRUE)
    mh2=sum(dnbinom(x=N,size=odp,,mu=lambda,log=TRUE),na.rm=TRUE)+
        dnorm(log(odp),l.odp.mn,l.odp.sd,log=TRUE)
    mh=exp(mh1-mh2)
    if(mh>runif(1)){
        odp=odp.star
        accept.odp=accept.odp+1
    }

    ##
    ## Sample p
    ##

    p.1999=rbeta(1,sum(Y.1999)+
                   q.p,sum(N.1999-Y.1999)+r.p)
    p.2000=rbeta(1,sum(Y.2000)+
                   q.p,sum(N.2000-Y.2000)+r.p)
    p.2001=rbeta(1,sum(Y.2001)+
                   q.p,sum(N.2001-Y.2001)+r.p)
    p.2002=rbeta(1,sum(Y.2002)+
                   q.p,sum(N.2002-Y.2002)+r.p)
    p.2003=rbeta(1,sum(Y.2003)+
                   q.p,sum(N.2003-Y.2003)+r.p)
    p.2004=rbeta(1,sum(Y.2004)+
                   q.p,sum(N.2004-Y.2004)+r.p)
    p.2006=rbeta(1,sum(Y.2006)+
                   q.p,sum(N.2006-Y.2006)+r.p)
    p.2012=rbeta(1,sum(Y.2012)+
                   q.p,sum(N.2012-Y.2012)+r.p)
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

    ##
    ## Derived parameters
    ##

    n.tot.v=rnbinom(n=length(lambda),size=odp,,mu=lambda*BoundaryInf.v)
    n.tot=unname(tapply(n.tot.v,(seq_along(n.tot.v)-1)%/%q,sum,na.rm=TRUE))

    ##
    ## Store results
    ##

    ## parameters
    gamma.save[k]=gamma
    beta.save[k,]=beta
    theta.save[k]=theta
    kappa.save[k]=kappa
    odp.save[k]=odp
    p.save[k,]=c(p.1999,p.2000,p.2001,p.2002,p.2003,p.2004,p.2006,p.2012)
    n.tot.save[k,]=n.tot

    ## tuners
    gamma.tune.save[k]=gamma.tune
    beta.tune.save[k,]=beta.tune
    theta.tune.save[k]=theta.tune
    kappa.tune.save[k]=kappa.tune
    odp.tune.save[k]=odp.tune

    ##
    ## Checkpoint
    ##

    if(k%%checkpoint==0){

        ##
        ## Update tuning parameters
        ##

        if(accept.gamma/k<0.3){
            gamma.tune=gamma.tune*0.9
        }
        if(accept.gamma/k>0.5){
            gamma.tune=gamma.tune*1.1
        }

        beta.tune=ifelse(accept.beta/k<0.3,
                         beta.tune*0.9,
                  ifelse(accept.beta/k>0.5,
                         beta.tune*1.1,
                         beta.tune))

        if(accept.theta/k<0.3){
            theta.tune=theta.tune*0.9
        }
        if(accept.theta/k>0.5){
            theta.tune=theta.tune*1.1
        }

        if(accept.kappa/k<0.3){
            kappa.tune=kappa.tune*0.9
        }
        if(accept.kappa/k>0.5){
            kappa.tune=kappa.tune*1.1
        }

        if(accept.odp/k<0.3){
            odp.tune=odp.tune*0.9
        }
        if(accept.odp/k>0.5){
            odp.tune=odp.tune*1.1
        }

        ##
        ## Save to list
        ##

        out=list(gamma.save,
                 beta.save,
                 gamma.tune.save,
                 beta.tune.save,
                 theta.save,
                 theta.tune.save,
                 kappa.save,
                 kappa.tune.save,
                 accept.N,
                 odp.save,
                 odp.tune.save,
                 p.save,
                 n.tot.save
                 )
        save(out,file=paste("~/Dropbox/Post-Doc/manuscripts/",
                            "SeaOtterEcology/DandD/Revision/",
                            "MCMCNBRealDataResults1.RData",
                            sep=""))
        cat(k,"")
    }

###
### End MCMC
###

}

############################################################################
###
### Results Files:
###        MCMCNBRealDataResults1.RData (12 Feb 2018): All parameters (worked)
###
############################################################################

###
### Load MCMC resullts
###

load(paste("~/Dropbox/Post-Doc/manuscripts/",
           "SeaOtterEcology/DandD/Revision/",
           "MCMCNBRealDataResults1.RData",
           sep=""))

(status=sum(!is.na(out[[1]])))

###
### Plot growth rate parameters and tuning values
###

par(mfrow=c(1,2))
gamma.est=out[[1]][1:status]
mean(gamma.est)
plot(gamma.est,type='l',ylim=c(0.15,0.25))
abline(h=0.2,col=2)
gamma.tune.save=out[[3]][1:status]
plot(gamma.tune.save,type='l')
## tail(gamma.tune.save)

###
### Plot motility parameters and tuning values
###

par(mfrow=c(5,2),mar=c(2,4,1,1))
beta0=17.72
beta1=-1.48
beta2=0.77
beta3=-0.35
beta4=1.0

beta.est=out[[2]][1:status,]
mean(beta.est[,1])
plot(beta.est[,1],type='l')
abline(h=beta0,col=2)
beta.tune.save=out[[4]][1:status,]
plot(beta.tune.save[,1],type='l')

mean(beta.est[,2])
plot(beta.est[,2],type='l')
abline(h=beta1,col=2)
plot(beta.tune.save[,2],type='l')

mean(beta.est[,3])
plot(beta.est[,3],type='l')
abline(h=beta2,col=2)
plot(beta.tune.save[,3],type='l')

mean(beta.est[,4])
plot(beta.est[,4],type='l')
abline(h=beta3,col=2)
plot(beta.tune.save[,4],type='l')

mean(beta.est[,5])
plot(beta.est[,5],type='l')
abline(h=beta4,col=2)
plot(beta.tune.save[,5],type='l')
## tail(beta.tune.save)

###
### Plot initial conditions
###

## Theta
par(mfrow=c(2,2),mar=c(4,3,1,1))
theta.est=out[[5]][1:status]
mean(theta.est)
plot(theta.est,type='l')
abline(h=500,col=2)
theta.tune.save=out[[6]][1:status]
plot(theta.tune.save,type='l')
## tail(theta.tune.save)

## Kappa
kappa.est=out[[7]][1:status]
mean(kappa.est)
plot(kappa.est,type='l')
abline(h=5.95,col=2)
kappa.tune.save=out[[8]][1:status]
plot(kappa.tune.save,type='l')
## tail(kappa.tune.save)

###
### Plot overdispersion parameter
###

odp.est=out[[10]][1:status]
mean(odp.est)
plot(odp.est,type='l',ylim=c(0,0.1))
abline(h=2,col=2)
odp.tune.save=out[[11]][1:status]
plot(odp.tune.save,type='l')
## tail(odp.tune.save)

###
### Calculate abundance
###

N.est=out[[13]][1:status,]
mean.n=apply(N.est,2,mean,na.rm=TRUE)
lb.n=apply(N.est,2,quantile,0.025)
ub.n=apply(N.est,2,quantile,0.975)
cbind(lb.n,mean.n,ub.n)

par(mfrow=c(1,1))
plot(1993:2012,mean.n,pch=16,ylim=c(0,15000))
segments(x0=1993:2012,
         y0=mean.n,
         x1=1993:2012,
         y1=lb.n)
segments(x0=1993:2012,
         y0=mean.n,
         x1=1993:2012,
         y1=ub.n)

counts[counts==0]=NA
points(1993:2012,counts,col=4,pch=16)
dbe=c(NA,NA,NA,NA,NA,NA,384,554,1238,1266,1866,2381,NA,2785,NA,NA,NA,NA,NA,8508)
dbeSE=c(NA,NA,NA,NA,NA,NA,111,97,143,196,458,594,NA,361,NA,NA,NA,NA,NA,2243)
dbelq=dbe+2*dbeSE
dbeuq=dbe-2*dbeSE
offs=0.1
points(1993:2012+offs,dbe,pch=16,col=2)
segments(x0=1993:2012+offs,
         y0=dbe,
         x1=1993:2012+offs,
         y1=dbelq,col=2)
segments(x0=1993:2012+offs,
         y0=dbe,
         x1=1993:2012+offs,
         y1=dbeuq,col=2)

###
### Calculate p
###

p.1999=0.7953388
p.2000=0.75
p.2001=0.8583887
p.2002=0.8503951
p.2003=0.7606767
p.2004=0.7694584
p.2006=0.7470478
p.2012=0.5833859
p.truth=c(p.1999,p.2000,p.2001,p.2002,p.2003,p.2004,p.2006,p.2012)

p.est=out[[12]][1:status,]
for(i in 1:8){
    plot(p.est[,i],type='l');abline(h=p.truth[i],col=2)
    readline()
}
