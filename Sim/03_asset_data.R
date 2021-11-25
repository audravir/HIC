rm(list=ls(all=TRUE))
library(readxl)
library(moments)
library(rugarch)
library(tseries)
library(sfsmisc)
library(zoo)
library(plyr)
library(edfun)

assets  = c( "SP500","FTSE100","IBOVESPA","NIKKEI")
dm      = length(assets)
df      = NULL
df$Date = seq(as.Date("1975/1/1"), as.Date(Sys.Date()), by = "day")


for (i in 1:dm){
    data  = read_excel(paste('asset_data/',assets[i],'.xlsx',sep=''))
    tmp   = 100*(log(data$PX_LAST[-1])-log(data$PX_LAST[-length(data$PX_LAST[-1])]))
    Date  = as.Date(data$Date[-1])
    tmpdf = data.frame(tmp,Date)
    colnames(tmpdf) = c(assets[i],"Date")
    df    = merge.data.frame(df,tmpdf,by='Date',all=TRUE)
}

df        = na.omit(df)

freq      = "month"

Week      = as.Date(cut(df$Date, freq))
indw      = unique(Week)
nnw       = length(indw)
ret_w     = RV_w = matrix(NA,ncol=dm,nrow=nnw)

for(i in 1:dm){
    for(t in 1:nnw){
        ret_w[t,i] = sum(df[which(Week==indw[t]),i+1])
        RV_w[t,i]  = sum((df[which(Week==indw[t]),i+1])^2)
    }
}

RV_w[which(RV_w==0)] = 0.0001


par(mfrow=c(2,ceiling(dm/2)))
for(i in 1:dm){
    plot(indw,ret_w[,i],type='l',main=assets[i])
}

par(mfrow=c(2,dm))
for(i in 1:dm){
    acf(ret_w[,i],lag.max = 50,ylim=c(-0.15,0.15))
}
for(i in 1:dm){
    pacf(ret_w[,i],lag.max = 50,ylim=c(-0.15,0.15))
}

head(indw)
tail(indw)

# Remove ARMA(1,1) to obtain de-meaned returns
# dmret_w and remove ARCH effects via GARCH


dmret_w = matrix(0,ncol=dm,nrow=nnw)
coef    = ts = NULL

for(i in 1:dm){
    model = arma(ret_w[,i],order = c(1, 1),include.intercept = TRUE)
    coef = rbind(coef,model$coef)
    ts   = rbind(ts,model$coef/sqrt(diag(model$vcov)))
    for(t in 2:nnw){
        dmret_w[t,i] = ret_w[t,i] - coef[i,3]-
        coef[i,1]*ret_w[t-1,i]-coef[i,2]*dmret_w[t-1,i]}

}

par(mfrow=c(2,dm))
for(i in 1:dm){
    acf(dmret_w[,i]^2,lag.max = 50,ylim=c(-0.15,0.3))
}
for(i in 1:dm){
    pacf(dmret_w[,i]^2,lag.max = 50,ylim=c(-0.15,0.3))
}

spec <- ugarchspec(variance.model = list(model = "sGARCH",
                                         garchOrder = c(1, 1),
                                         submodel = NULL,
                                         external.regressors = NULL,
                                         variance.targeting = FALSE),

                   mean.model     = list(armaOrder = c(0, 0),
                                         external.regressors = NULL,
                                         include.mean=FALSE),
                   distribution.model = "norm")


resid_garch = hs_garch = matrix(NA,ncol=dm,nrow=nnw)
resid_manual = hs_manual = matrix(NA,ncol=dm,nrow=nnw)
coef    = ts = NULL

for(i in 1:dm){
    garch <- ugarchfit(spec = spec, data = dmret_w[,i], solver.control = list(trace=0))
    coef      = rbind(coef,garch@fit$coef)
    ts       = rbind(ts,garch@fit$tval)
    hs_garch[,i]   = garch@fit$sigma
    resid_garch[,i] = garch@fit$z
}


for(i in 1:dm){
    hs_manual[1,i] = var(dmret_w[,i])*((nnw-1)/nnw)
    for(t in 2:nnw)
        hs_manual[t,i]=coef[i,1]+coef[i,2]*dmret_w[t-1,i]^2+coef[i,3]*hs_manual[t-1,i]
}

par(mfrow=c(2,ceiling(dm/2)))
for(i in 1:dm){
    plot(hs_manual[,i],type='l')
    lines((hs_garch[,i])^2,col=2)
}


par(mfrow=c(2,ceiling(dm/2)))
for(i in 1:dm){
    plot(resid_garch[,i],type='l')
    lines(dmret_w[,i]/sqrt(hs_manual[,i]),col=2)
}

garch_res = list(coef=coef,ts=ts,resid=resid_garch)

apply(garch_res$resid ,2,mean)
apply(garch_res$resid ,2,sd)






##--------------------
## Using RV
##--------------------

RV1_resid = RV_resid = ret_w/sqrt(RV_w)

par(mfrow=c(2,dm))
for(i in 1:dm){
    acf(RV1_resid[,i]^2,lag.max = 50,ylim=c(-0.15,0.15))
}
for(i in 1:dm){
    pacf(RV1_resid[,i]^2,lag.max = 50,ylim=c(-0.15,0.15))
}

apply(RV1_resid, 2,mean)
apply(RV1_resid, 2,sd)

coef=ts=NULL
for(i in 1:dm){
    model = arma(RV1_resid[,i],order = c(1, 1),include.intercept = TRUE)
    coef = rbind(coef,model$coef)
    ts   = rbind(ts,model$coef/sqrt(diag(model$vcov)))
    for(t in 2:nnw){
        dmret_w[t,i] = RV1_resid[t,i] - coef[i,3]-
            coef[i,1]*RV1_resid[t-1,i]-coef[i,2]*dmret_w[t-1,i]}

    RV_resid[,i] =  dmret_w[,i]/sd(dmret_w[,i])
}

RV_res = list(coef=coef,ts=ts,resid=RV_resid)

apply(RV_res$resid ,2,mean)
apply(RV_res$resid ,2,sd)

apply(garch_res$resid ,2,mean)
apply(garch_res$resid ,2,sd)


par(mfrow=c(2,dm))
for(i in 1:dm){
    acf(RV_res$resid[,i]^2)
}
for(i in 1:dm){
    pacf(RV_res$resid[,i]^2)
}



resid = garch_res$resid


desc = matrix(NA,ncol=14,nrow=dm)
for(i in 1:dm){
    desc[i,1] = mean(resid[,i])
    desc[i,2] = median(resid[,i])
    desc[i,3] = sd(resid[,i])
    desc[i,4] = skewness(resid[,i])
    desc[i,5] = kurtosis(resid[,i])
    t = shapiro.test(resid[,i])
    desc[i,6] = t$p.value
    t = ks.test(resid[,i],y='pnorm')
    desc[i,7] = t$p.value
    t = jarque.test(resid[,i])
    desc[i,8] = t$p.value

    t = Box.test (resid[,i], lag = 1, type = "Ljung")
    desc[i,9] = t$p.value
    t = Box.test (resid[,i], lag = 5, type = "Ljung")
    desc[i,10] = t$p.value
    t = Box.test (resid[,i], lag = 10, type = "Ljung")
    desc[i,11] = t$p.value

    t = Box.test (resid[,i]^2, lag = 1, type = "Ljung")
    desc[i,12] = t$p.value
    t = Box.test (resid[,i]^2, lag = 5, type = "Ljung")
    desc[i,13] = t$p.value
    t = Box.test (resid[,i]^2, lag = 10, type = "Ljung")
    desc[i,14] = t$p.value
}

rownames(desc) = assets
desc           = cbind(desc)
colnames(desc) = c('mean','median','sd','skew.','kurt.',
                   'SW','KS','JB','LB(1)','LB(5)','LB(10)',
                   'LB(1)','LB(5)','LB(10)')

round(desc,2)

grid=seq(-5,5,length=1000)

par(mfrow=c(2,ceiling(dm/2)))
for(i in 1:dm) {hist((resid[,i]),main=assets[i],
                     ylab = '',xlab = '',freq=FALSE)
    lines(grid,dnorm(grid),lwd=2)}


par(mfrow=c(2,ceiling(dm/2)))
for(i in 1:dm) {hist(pnorm(resid[,i]),main=assets[i],
                     ylab = '',xlab = '',freq=FALSE)
    abline(h=1,lwd=2)}


par(mfrow=c(2,ceiling(dm/2)))
for(i in 1:dm) {plot((resid[,i]),main=assets[i],
                     ylab = '',xlab = '',type='l')}




###


lrets = as.matrix(df[,-1])

RCov = RCor = array(NA,c(dm,dm,nnw))

for(t in 1:nnw){
    RCov[,,t] = t(lrets[which(Week==indw[t]),])%*%lrets[which(Week==indw[t]),]
    RCov[which(diag(RCov[,,t])==0),which(diag(RCov[,,t])==0),t] = 0.0001
    RCor[,,t] = nearcor(round(cov2cor(RCov[,,t]),8))$cor
}

Sig    = alply(RCor,3)
Sbar   = Reduce('+',Sig)/length(Sig)

C = combinat::combn(1:dm,2)

RCroll1 = RCroll2 = matrix(NA,ncol=ncol(C),nrow=nnw)

for(i in 1:ncol(C)){
    for(t in 2:nnw){
        RCroll1[t,i] = cor(ret_w[max(c(1,t-20)):t,C[1,i]],
                           ret_w[max(c(1,t-20)):t,C[2,i]])
        RCroll2[t,i] = cor(u_garch_N[max(c(1,t-20)):t,C[1,i]],
                           u_garch_N[max(c(1,t-20)):t,C[2,i]])
    }
}



#pdf(paste('RCor_',freq,'.pdf',sep=''),width=20,height=10)
par(mfrow=c(3,ceiling(ncol(C)/3)))
for(i in 1:ncol(C)){
    plot(indw, RCor[C[1,i],C[2,i],],type='l',ylim=c(-1.2,1.2),
         main=c(assets[C[1,i]],assets[C[2,i]]),ylab='',
         xlab='')
    lines(indw,rollapply(RCor[C[1,i],C[2,i],],
                         width=13,FUN=mean,fill=c(NA,NA,NA),
                         align = c("center")),col=2,lwd=2)
    abline(h=0,lwd=2,lty=2,col='gray80')
    abline(h=Sbar[C[1,i],C[2,i]],lwd=2,col=3)
    lines(indw,RCroll1[,i],col=4)
    lines(indw,RCroll2[,i],col=6)
}
#dev.off()

u_garch_N = pnorm(garch_res$resid)
u_RV_N    = pnorm(RV_res$resid)

par(mfrow=c(2,ceiling(dm/2)))
for(i in 1:dm){
    hist(u_garch_N[,i],main=assets[i])
}

par(mfrow=c(2,ceiling(dm/2)))
for(i in 1:dm){
    hist(u_RV_N[,i],main=assets[i])
}

dist_RV=dist_garch=list()
u_garch_e=u_RV_e=matrix(NA,ncol=dm,nrow=nnw)

for(i in 1:dm){
    res     = list((edfun(RV_res$resid[,i])))
    dist_RV = c(dist_RV,res)
    u_RV_e[,i] = res[[1]]$pfun(RV_res$resid[,i])
    res     = list((edfun(garch_res$resid[,i])))
    dist_garch = c(dist_garch,res)
    u_garch_e[,i] = res[[1]]$pfun(garch_res$resid[,i])
}

par(mfrow=c(2,ceiling(dm/2)))
for(i in 1:dm){
    hist(u_RV_e[,i],main=assets[i])
}

par(mfrow=c(2,ceiling(dm/2)))
for(i in 1:dm){
    hist(u_garch_e[,i],main=assets[i])
}

data = list(Date=indw,RCor=RCor,lret=ret_w,
            u_garch_e=u_garch_e,u_garch_N=u_garch_N,
            u_RV_e=u_RV_e,u_RV_N=u_RV_N)
names(data)

#save(data,file=paste('data_',freq,'.Rdata',sep=''))

