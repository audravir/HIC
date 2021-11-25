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

rm(tmpdf,data)

df    = na.omit(df)

freq  = "month"

FQ   = as.Date(cut(df$Date, freq))
indw = unique(FQ)
nn   = length(indw)
ret  = RV = matrix(NA,ncol=dm,nrow=nn)

for(i in 1:dm){
    for(t in 1:nn){
        ret[t,i] = sum(df[which(FQ==indw[t]),i+1])
        RV[t,i]  = sum((df[which(FQ==indw[t]),i+1])^2)
    }
}

RV[which(RV==0)] = 0.0001


par(mfrow=c(2,ceiling(dm/2)))
for(i in 1:dm){
    plot(indw,ret[,i],type='l',main=assets[i])
}

par(mfrow=c(2,dm))
for(i in 1:dm){
    acf(ret[,i],lag.max = 50,ylim=c(-0.15,0.15),
        main=assets[i])
}
for(i in 1:dm){
    pacf(ret[,i],lag.max = 50,ylim=c(-0.15,0.15),
         main=assets[i])
}

head(indw)
tail(indw)


###----
## RCov and RCor
###----

lrets = as.matrix(df[,-1])

RCov = RCor = array(NA,c(dm,dm,nn))

for(t in 1:nn){
    RCov[,,t] = t(lrets[which(FQ==indw[t]),])%*%lrets[which(FQ==indw[t]),]
    RCov[which(diag(RCov[,,t])==0),which(diag(RCov[,,t])==0),t] = 0.0001
    RCor[,,t] = nearcor(round(cov2cor(RCov[,,t]),8))$cor
}

Sig    = alply(RCor,3)
Sbar   = Reduce('+',Sig)/length(Sig)

C = combinat::combn(1:dm,2)

#pdf(paste('RCor_',freq,'.pdf',sep=''),width=20,height=10)
par(mfrow=c(2,ceiling(ncol(C)/2)))
for(i in 1:ncol(C)){
    plot(indw, RCor[C[1,i],C[2,i],],type='l',ylim=c(-1.2,1.2),
         main=c(assets[C[1,i]],assets[C[2,i]]),ylab='',
         xlab='')
    lines(indw,rollapply(RCor[C[1,i],C[2,i],],
                         width=13,FUN=mean,fill=c(NA,NA,NA),
                         align = c("center")),col=2,lwd=2)
    abline(h=0,lwd=2,lty=2,col='gray80')
    abline(h=Sbar[C[1,i],C[2,i]],lwd=2,col=3)
    abline(h=c(-1,1),lty=2,col='gray80')
}
#dev.off()

# Check if diagonal of RCov is the same as RV
# Plots should coincide
par(mfrow=c(2,ceiling(dm/2)))
for(i in 1:dm){
    plot(indw,sqrt(RV[,i]),type='l',main=assets[i],lwd=3)
    lines(indw,sqrt(RCov[i,i,]),col=2,lwd=2)
}

data = list(date=indw,ret=ret,RV=RV,RCor=RCor,assets=assets)
save(data,file='data.RData')
