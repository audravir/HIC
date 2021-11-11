rm(list=ls(all=TRUE))

assets  = c( "BAC", "JPM",  "IBM",  "MSFT",
             "XOM",  "AA",   "AXP",  "DD",
             "GE",   "KO" )
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

freq      = "week"

Week      = as.Date(cut(df$Date, freq))
indw      = unique(Week)
nnw       = length(indw)
ret_w     = RV_w = matrix(NA,ncol=dm,nrow=nnw)

for(i in 1:dm){
    for(t in 1:nnw){
        ret_w[t,i] = sum(df[which(Week==indw[t]),i+1])
        RV_w[t,i] = sum((df[which(Week==indw[t]),i+1])^2)
    }
}

par(mfrow=c(2,ceiling(dm/2)))
for (i in 1:dm){
    plot(indw,ret_w[,i],type='l',main=assets[i])
    lines(indw,sqrt(RV_w[,i]),col=2,lwd=2)
}


stand_w = ret_w/sqrt(RV_w)

ar_order  = rep(NA,dm)
for(i in 1:dm){
    model = ar(stand_w[,i],aic = TRUE)
    stand_w[,i] = model$resid
    ar_order[i] = model$order
}


stand_w = stand_w/matrix(rep(apply(stand_w, 2, sd,na.rm=TRUE),nnw),nrow=nnw,byrow = TRUE)
stand_w = stand_w[complete.cases(stand_w),]


desc = matrix(NA,ncol=14,nrow=dm)
for(i in 1:dm){
    desc[i,1] = mean(stand_w[,i])
    desc[i,2] = median(stand_w[,i])
    desc[i,3] = sd(stand_w[,i])
    desc[i,4] = skewness(stand_w[,i])
    desc[i,5] = kurtosis(stand_w[,i])
    t = shapiro.test(stand_w[,i])
    desc[i,6] = t$p.value
    t = ks.test(stand_w[,i],y='pnorm')
    desc[i,7] = t$p.value
    t = jarque.test(stand_w[,i])
    desc[i,8] = t$p.value

    t = Box.test (stand_w[,i], lag = 1, type = "Ljung")
    desc[i,9] = t$p.value
    t = Box.test (stand_w[,i], lag = 5, type = "Ljung")
    desc[i,10] = t$p.value
    t = Box.test (stand_w[,i], lag = 10, type = "Ljung")
    desc[i,11] = t$p.value

    t = Box.test (stand_w[,i]^2, lag = 1, type = "Ljung")
    desc[i,12] = t$p.value
    t = Box.test (stand_w[,i]^2, lag = 5, type = "Ljung")
    desc[i,13] = t$p.value
    t = Box.test (stand_w[,i]^2, lag = 10, type = "Ljung")
    desc[i,14] = t$p.value
}

rownames(desc) = assets
desc           = cbind(desc,ar_order)
colnames(desc) = c('mean','median','sd','skew.','kurt.',
                   'SW','KS','JB','LB(1)','LB(5)','LB(10)',
                   'LB(1)','LB(5)','LB(10)','AR order')

round(desc,2)


par(mfrow=c(2,ceiling(dm/2)))
for(i in 1:dm) {hist(pnorm(stand_w[,i]),main=assets[i],
                     ylab = '',xlab = '')}


par(mfrow=c(2,ceiling(dm/2)))
for(i in 1:dm) {plot((stand_w[,i]),main=assets[i],
                     ylab = '',xlab = '',type='l')}

lrets = as.matrix(df[,-1])

RCov = RCor = array(NA,c(dm,dm,nnw))

for(t in 1:nnw){
    RCov[,,t] = t(lrets[which(Week==indw[t]),])%*%lrets[which(Week==indw[t]),]
    RCor[,,t] = cov2cor(RCov[,,t])
}


C = combinat::combn(1:dm,2)

par(mfrow=c(2,ceiling(ncol(C)/2)))
for(i in 1:ncol(C)){
    plot(indw, RCor[C[1,i],C[2,i],],type='l',ylim=c(-1,1),
         main=c(assets[C[1,i]],assets[C[2,i]]))
}



