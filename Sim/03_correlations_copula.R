rm(list=ls(all=TRUE))
load('data_u.Rdata')

nn      = dim(data$ret)[1]
dm      = dim(data$ret)[2]
oos     = 5*12 # last 5 years


## Look at rolling window yearly correlations vs RCor

C       = combinat::combn(1:dm,2)
RCroll1 = RCroll2 = matrix(NA,ncol=ncol(C),nrow=nn)

for(i in 1:ncol(C)){
    for(t in 2:nn){
        RCroll1[t,i] = cor(data$u_arma_garch[max(c(1,t-11)):t,C[1,i]],
                           data$u_arma_garch[max(c(1,t-11)):t,C[2,i]])
        RCroll2[t,i] = cor(data$u_arma_rv[max(c(1,t-11)):t,C[1,i]],
                           data$u_arma_rv[max(c(1,t-11)):t,C[2,i]])
    }
}

par(mfrow=c(3,ceiling(ncol(C)/3)))
for(i in 1:ncol(C)){
    plot(RCroll1[,i],type='l',ylim=c(-1,1),main=cor(RCroll1[-1,i],data$RCor[C[1,i],C[2,i],2:nn]))
    lines(RCroll2[,i],col=2)
    lines(data$RCor[C[1,i],C[2,i],],col=4)
    abline(h=c(-1,0,1))
}


##----
## DCC-Gaussian, DCC-t and AIW copulas
## estimate for in-sample period only
## use fixed parameters to roll for out-of-sample
## produce t+1 forecasts (N*dm matrices of u)
## that will be used in portfolio
##----
# inputs:
# u
# M
# oos
# post_sample
u=data$u_arma_garch
oos=5*12
M=100
post_sample=20 # cannot be larger than M

###


if (post_sample>M)stop('posterior sample larger than n.iter. M')
nn     <- dim(u)[1]
nnis   <- nn-oos
dm     <- dim(u)[2]
x_all  <- qnorm(u)
x_is   <- x_all[1:nnis,]
bi     <- M
Qold   <- array(NA,c(dm, dm, nnis))
R      <- array(NA,c(dm, dm, nnis))
Qold[,,1] <- cor(x_is)
R[,,1] <- diag(diag(Qold[,,1])^{-1/2})%*%Qold[,,1]%*%diag(diag(Qold[,,1])^{-1/2})
aold   <- 0.01
bold   <- 0.95
llold  <- rep(0,nnis)
resdcc <- matrix(NA,ncol=2,nrow=(bi+M))
accdcc <- rep(0,bi+M)
Vpred  = Qpred = list()
S      <- cov(x_is)

for(t in 2:nnis){
    Qold[,,t]   <- S*(1-aold-bold)+aold*(x_is[t-1,]%*%t(x_is[t-1,]))+bold*Qold[,,(t-1)]
    R[,,t]   <- diag(diag(Qold[,,t])^{-1/2})%*%Qold[,,t]%*%diag(diag(Qold[,,t])^{-1/2})
    llold[t] <- mvtnorm::dmvnorm(x_is[t,], rep(0,dm), R[,,t], log=T)
}

for(m in 1:(M+bi)){

    parnew = rnorm(2,c(aold,bold),0.005) # changed from 0.0005

    anew  <- parnew[1]
    bnew  <- parnew[2]
    llnew <- rep(0,nnis)
    Qnew  = Qold

    for(t in 2:nnis){
        Qnew[,,t]   <- S*(1-anew-bnew)+anew*(x_is[t-1,]%*%t(x_is[t-1,]))+bnew*Qnew[,,(t-1)]
        R[,,t]   <- diag(diag(Qnew[,,t])^{-1/2})%*%Qnew[,,t]%*%diag(diag(Qnew[,,t])^{-1/2})
        llnew[t] <- mvtnorm::dmvnorm(x_is[t,], rep(0,dm), R[,,t], log=T)
    }

    if((sum(llnew)-sum(llold)+
        dbeta(anew,3,10,log=T)-dbeta(aold,3,10,log=T)+
        dbeta(bnew,10,3,log=T)-dbeta(bold,10,3,log=T))>log(runif(1))
       &&
       (sum(parnew)<1))
    {
        llold  = llnew
        aold   = anew
        bold   = bnew
        Qold   = Qnew
        accdcc[m] = 1
    }

    resdcc[m,] <- c(aold,bold)
    Qpred[[m]] <- S*(1-aold-bold)+aold*(x_is[nnis,]%*%t(x_is[nnis,]))+bold*Qold[,,nnis]
    Vpred[[m]] <- diag(diag(Qpred[[m]])^{-1/2})%*%Qpred[[m]]%*%diag(diag(Qpred[[m]])^{-1/2})

    # if(m>bi){save(Qold,file=paste('temp/dcc_',m-bi,'.Rdata',sep=''))}
}

ind_m=round(seq(1,M,length=post_sample))

res = list(resdcc[(ind_m+bi),],accdcc[ind_m+bi],oos,x_all)
names(res) = c('resdcc','accdcc','oos','X')
save(res,file=paste('temp/results_scalar_dcc.Rdata',sep=''))

## Roll and predict the correlations

rm(list=ls(all=TRUE))

load('temp/results_scalar_dcc.Rdata')

nn      = dim(res$X)[1]
dm      = dim(res$X)[2]
oos     = res$oos
x_is    = res$X[1:(nn-oos),]
miter   = length(res$accdcc)

Q_all=R_all=array(NA,c(dm, dm, nn,miter))
Q_all[,,1,] <- cor(x_is)
R_all[,,1,] <- diag(diag(Q_all[,,1,1])^{-1/2})%*%Q_all[,,1,1]%*%diag(diag(Q_all[,,1,1])^{-1/2})
S          <- cov(x_is)


a = res$resdcc[,1]
b = res$resdcc[,2]

for(m in 1:miter){
    for(t in 2:nn){
        Q_all[,,t,m]   <- S*(1-a[m]-b[m])+a[m]*(res$X[t-1,]%*%t(res$X[t-1,]))+b[m]*Q_all[,,(t-1),m]
        R_all[,,t,m]   <- diag(diag(Q_all[,,t,m])^{-1/2})%*%Q_all[,,t,m]%*%diag(diag(Q_all[,,t,m])^{-1/2})
    }
}



