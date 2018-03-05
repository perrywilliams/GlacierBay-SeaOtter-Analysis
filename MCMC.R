#############################################################################
#############################################################################
###
### Function to fit the Glacier Bay model to data using MCMC
###
#############################################################################
#############################################################################

MCMC=function(data,
              priors,
              inits,
              parameters,
              st.info,
              Bathymetry,
              n.iter=5000,
              checkpoint=1000,
              output.location){

    ##
    ##  Subroutines and Packages
    ##

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
    calcc = cxxfunction(signature(
        H="numeric",c0="numeric",timesteps="numeric"),
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

    ##
    ## Priors
    ##

    q.gamma=priors$gamma.prior[1]
    r.gamma=priors$gamma.prior[2]
    beta.mean=priors$beta.prior[1]
    beta.var=priors$beta.prior[2]
    theta.mean=priors$theta.prior[1]
    theta.var=priors$theta.prior[2]
    kappa.mean=priors$kappa.prior[1]
    kappa.var=priors$kappa.prior[2]
    q.p=priors$p.prior[1]
    r.p=priors$p.prior[2]
    q.tau=priors$tau.prior[1]
    r.tau=priors$tau.prior[2]

    ##
    ## Dimensions
    ##

    dt=st.info$dt
    time.frame=st.info$time.frame
    us.fact=st.info$us.fact
    res=st.info$res
    d=st.info$d
    time.steps=1/dt*length(time.frame)
    keep=seq(1,time.steps,time.steps/length(time.frame))
    dx=res*us.fact
    dy=res*us.fact
    xmin=st.info$extent[1]
    xmax=st.info$extent[2]
    ymin=st.info$extent[3]
    ymax=st.info$extent[4]
    y=170
    x=140
    q=x*y

    ##
    ## Boundary layers
    ##

    Boundary=Bathymetry
    Boundary[is.na(Bathymetry)]=0
    Boundary[Bathymetry<0]=1
    Boundary[Bathymetry>=0]=0
    BoundaryNA=Boundary
    BoundaryNA[Boundary==0]=NA
    BoundaryInf=BoundaryNA
    BoundaryInf[120:170,0:80]=NA
    BoundaryInf[150:170,0:140]=NA
    BoundaryInf[134:170,116:140]=NA
    BoundaryInf.v=rep(BoundaryInf[],length(keep))
    Boundary.us=aggregate(Boundary,fact=us.fact,na.rm=TRUE,fun=max)
    Boundary.us[Boundary.us[]>1]=1
    Boundary.us[1,]=0
    Boundary.us[,1]=0
    Boundary.us[dim(Boundary.us)[1],]=0
    Boundary.us[,dim(Boundary.us)[2]]=0

    ##
    ## Data
    ##

    Y=data$Y
    X=data$X

    Y.1999=data$ISU$Y.1999
    N.1999=data$ISU$N.1999
    Y.2000=data$ISU$Y.2000
    N.2000=data$ISU$N.2000
    Y.2001=data$ISU$Y.2001
    N.2001=data$ISU$N.2001
    Y.2002=data$ISU$Y.2002
    N.2002=data$ISU$N.2002
    Y.2003=data$ISU$Y.2003
    N.2003=data$ISU$N.2003
    Y.2004=data$ISU$Y.2004
    N.2004=data$ISU$N.2004
    Y.2006=data$ISU$Y.2006
    N.2006=data$ISU$N.2006
    Y.2012=data$ISU$Y.2012
    N.2012=data$ISU$N.2012

    ##
    ## Starting values
    ##

    gamma=inits$gamma
    beta=inits$beta
    theta=inits$theta
    kappa=inits$kappa
    p=inits$p
    odp=inits$odp

    accept.gamma=0
    accept.beta=rep(0,length(beta))
    accept.theta=0
    accept.kappa=0
    accept.N=rep(0,length(Y))
    accept.odp=0

    delta=raster(,nrows=y,ncols=x,xmn=xmin,xmx=xmax,
                 ymn=ymin,ymx=ymax,crs=NA)
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

    cell=raster(,nrows=y,ncols=x,xmn=xmin,xmx=xmax,
                ymn=ymin,ymx=ymax,crs=NA)
    cell[]=1:q
    lambda0=raster(,nrows=y,ncols=x,xmn=xmin,xmx=xmax,
                   ymn=ymin,ymx=ymax,crs=NA)
    c0=raster(,nrows=y/us.fact,ncols=x/us.fact,
              xmn=xmin,xmx=xmax,ymn=ymin,ymx=ymax,crs=NA)
    D=rdist(data.frame(SpatialPoints(lambda0)),matrix(d,1,2))/1000
    us.cells=aggregate(cell,fact=us.fact,fun=mean)
    us.cells[]=1:length(us.cells[])
    lambda0[]=(exp(-D^2/kappa^2)/sum(exp(-D^2/kappa^2))*theta)
    lambda0.star=lambda0
    c0[]=extract(delta*lambda0,SpatialPoints(us.cells))
    c0.star=c0

    ## c.all and lambda.all and N.all
    c.all=brick(nrows=dim(c0)[1],ncols=dim(c0)[2],xmn=xmin,xmx=xmax,
                ymn=ymin,ymx=ymax)
    c.all=setValues(c.all, calcc(H, vec(c0[]),time.steps)[,keep])
    lambda=vec(disaggregate(c.all,us.fact)/delta)
    N=rnbinom(length(lambda),size=odp,,mu=lambda)
    N=ifelse(N<Y,Y,N)

    ##
    ## Tuning parameters
    ##

    gamma.tune=0.002
    beta.tune=c(0.1392345, 0.04616871, 0.1127799, 0.05642842, 0.03912024)
    theta.tune=40
    kappa.tune=1
    odp.tune=0.024

    ##
    ## Containers
    ##

    MCMC.Chains=vector('list',length(parameters))
    names(MCMC.Chains)=parameters
    if('gamma'%in%parameters){
        MCMC.Chains$gamma=matrix(,n.iter,1)
    }
    if('beta'%in%parameters){
        MCMC.Chains$beta=matrix(,n.iter,length(beta))
    }
    if('theta'%in%parameters){
        MCMC.Chains$theta=matrix(,n.iter,1)
    }
    if('kappa'%in%parameters){
        MCMC.Chains$kappa=matrix(,n.iter,1)
    }
    if('p'%in%parameters){
        MCMC.Chains$p=matrix(,n.iter,8)
    }
    if('odp'%in%parameters){
        MCMC.Chains$odp=matrix(,n.iter,1)
    }
    if('n.tot'%in%parameters){
        MCMC.Chains$n.tot=matrix(,n.iter,length(keep))
    }

    ##
    ## Gibbs loop
    ##

    for(k in 1:n.iter){

        ##
        ## Sample gamma
        ##

        gamma.star=rnorm(1,gamma,gamma.tune)
        if(gamma.star>q.gamma&gamma.star<r.gamma){
            gamma.bar.star=delta.bar*aggregate(gamma.star/delta,
                                               fact=us.fact,
                                               fun=function(x,na.rm){
                                                   mean(x,na.rm=TRUE)
                                               })
            H.star=propagator.plainZF(NN,delta.bar[],gamma.bar.star[],
                                      dx=dx,dy=dy,dt=dt)
            if(min(range(H.star,na.rm=TRUE))>=0){
                c.all.star=setValues(
                    c.all,calcc(H.star,vec(c0[]),time.steps)[,keep])
                lambda.star=vec(disaggregate(c.all.star,us.fact)/delta)
                mh1=sum(dnbinom(x=N,size=odp,,mu=lambda.star,log=TRUE),
                        na.rm=TRUE)
                mh2=sum(dnbinom(x=N,size=odp,,mu=lambda,log=TRUE),
                        na.rm=TRUE)
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
        ## Sample beta
        ##

        for(i in 1:length(beta)){
            beta.tmp=rnorm(1,mean=beta[i],sd=beta.tune[i])
            beta.star=beta
            beta.star[i]=beta.tmp
            delta.star[]=exp(X%*%beta.star)
            delta.bar.star=aggregate(delta.star,
                                     fact=us.fact,
                                     fun=function(x,na.rm){
                                         (1/mean(1/x,na.rm=TRUE))
                                     }
                                     )
            gamma.bar.star=delta.bar.star*aggregate(gamma/delta.star,
                                                    fact=us.fact,
                                                    fun=function(x,na.rm){
                                                        mean(x,na.rm=TRUE)
                                                    })
            H.star=propagator.plainZF(NN,delta.bar.star[],gamma.bar.star[],
                                      dx=dx,dy=dy,dt=dt)
            if(min(H.star,na.rm=TRUE)>=0){
                c0.star[]=extract(delta.star*lambda0,
                                  SpatialPoints(us.cells))
                c.all.star=setValues(
                    c.all,calcc(H.star,vec(c0.star[]),time.steps)[,keep])
                lambda.star=vec(disaggregate(c.all.star,us.fact)/delta.star)
                mh1=sum(dnbinom(x=N,size=odp,,mu=lambda.star,log=TRUE),
                        na.rm=TRUE)+
                    sum(dnorm(beta.star,beta.mean,beta.var^0.5,log=TRUE))
                mh2=sum(dnbinom(x=N,size=odp,,mu=lambda,log=TRUE),
                        na.rm=TRUE)+
                    sum(dnorm(beta,beta.mean,beta.var^0.5,log=TRUE))
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
                    accept.beta[i]=accept.beta[i]+1
                }
            }
        }

        ##
        ## Sample theta
        ##

        theta.star=rnorm(1,theta,theta.tune)
        if(theta.star>0){
            lambda0.star[]=exp(-D^2/kappa^2)/
                sum(exp(-D^2/kappa^2))*theta.star
            c0.star=extract(delta*lambda0.star,
                            SpatialPoints(us.cells))
            c.all.star=setValues(
                c.all,calcc(H,vec(c0.star[]),time.steps)[,keep])
            lambda.star=vec(disaggregate(c.all.star,us.fact)/delta)
            mh1=sum(dnbinom(x=N,size=odp,,mu=lambda.star,log=TRUE),
                    na.rm=TRUE)+
                dnorm(theta.star,theta.mean,theta.var^0.5,log=TRUE)
            mh2=sum(dnbinom(x=N,size=odp,,mu=lambda,log=TRUE),na.rm=TRUE)+
                dnorm(theta,theta.mean,theta.var^0.5,log=TRUE)
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
            lambda0.star[]=exp(-D^2/kappa.star^2)/
                sum(exp(-D^2/kappa.star^2))*theta
            c0.star=extract(delta*lambda0.star,
                            SpatialPoints(us.cells))
            c.all.star=setValues(
                c.all,calcc(H,vec(c0.star[]),time.steps)[,keep])
            lambda.star=vec(disaggregate(c.all.star,us.fact)/delta)
            mh1=sum(dnbinom(x=N,size=odp,,mu=lambda.star,log=TRUE),
                    na.rm=TRUE)+
                dnorm(kappa.star,kappa.mean,kappa.var^0.5,log=TRUE)
            mh2=sum(dnbinom(x=N,size=odp,,mu=lambda,log=TRUE),na.rm=TRUE)+
                dnorm(kappa,kappa.mean,kappa.var^0.5,log=TRUE)
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

        odp.star=rnorm(1,odp,odp.tune)
        if(odp.star>0&odp.star<1){
            mh1=sum(dnbinom(x=N,size=odp.star,,mu=lambda,log=TRUE),
                    na.rm=TRUE)
            mh2=sum(dnbinom(x=N,size=odp,,mu=lambda,log=TRUE),na.rm=TRUE)
            mh=exp(mh1-mh2)
            if(mh>runif(1)){
                odp=odp.star
                accept.odp=accept.odp+1
            }
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
        n.tot=unname(
            tapply(n.tot.v,(seq_along(n.tot.v)-1)%/%q,sum,na.rm=TRUE)
        )

        ##
        ## Store results
        ##

        if('gamma'%in%parameters){
            MCMC.Chains$gamma[k,]=gamma
        }
        if('beta'%in%parameters){
            MCMC.Chains$beta[k,]=beta
        }
        if('theta'%in%parameters){
            MCMC.Chains$theta[k,]=theta
        }
        if('kappa'%in%parameters){
            MCMC.Chains$kappa[k,]=kappa
        }
        if('p'%in%parameters){
            MCMC.Chains$p[k,]=c(p.1999,p.2000,p.2001,p.2002,p.2003,
                                p.2004,p.2006,p.2012)
        }
        if('odp'%in%parameters){
            MCMC.Chains$odp[k,]=odp
        }
        if('n.tot'%in%parameters){
            MCMC.Chains$n.tot[k,]=n.tot
        }

        ## ## tuners
        ## gamma.tune.save[k]=gamma.tune
        ## beta.tune.save[k,]=beta.tune
        ## theta.tune.save[k]=theta.tune
        ## kappa.tune.save[k]=kappa.tune
        ## odp.tune.save[k]=odp.tune

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

            save(MCMC.Chains,file=output.location)

            cat(k,"")
        }
    }
}

