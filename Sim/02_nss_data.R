rm(list=ls(all=TRUE))
library(readr)
library(moments)
library(OpenMx)
library(sfsmisc)
library(MASS)


data  = read_csv("nss_data/10_dim_daily_return.csv")
nn    = dim(data)[1]
dm    = 10 # choose the dimension
cc    = 1 # =1 open to close; =0 close to close
# ALWAYS use =1, after standardizing white noise
rets  = data.matrix(data[(dim(data)[1]-nn+1):(dim(data)[1]),(cc*10+2):(cc*10+1+dm)])
date  = as.Date(data$Date[(dim(data)[1]-nn+1):(dim(data)[1])],format="%Y-%m-%d")
names = c( "BAC", "JPM",  "IBM",  "MSFT",
           "XOM",  "AA",   "AXP",  "DD",
           "GE",   "KO" )


desc = matrix(NA,ncol=10,nrow=10)

for(i in 1:dm){
    desc[i,1] = mean(rets[,i])
    desc[i,2] = median(rets[,i])
    desc[i,3] = sd(rets[,i])
    desc[i,4] = skewness(rets[,i])
    desc[i,5] = kurtosis(rets[,i])
    t = shapiro.test(rets[,i])
    desc[i,6] = t$p.value
    t = ks.test(rets[,i],y='pnorm')
    desc[i,7] = t$p.value
    t = jarque.test(rets[,i])
    desc[i,8] = t$p.value
    t = Box.test (rets[,i], lag = 5, type = "Ljung")
    desc[i,9] = t$p.value
    t = Box.test (rets[,i], lag = 10, type = "Ljung")
    desc[i,10] = t$p.value
}

rownames(desc) = names
colnames(desc) = c('mean','median','sd','skew.','kurt.',
                   'SW','KS','JB','LB(5)','LB(10)')
round(desc,5)



##------------------------
## load Realized Measures
##------------------------

data  = read_csv("nss_data/10_dim_realized_covar.csv")
reaco = data.matrix(data[(dim(data)[1]-nn+1):(dim(data)[1]),-1])
RCov  = RCor = array(NA,c(dm,dm,nn))
stand = RVs = matrix(NA,ncol=dm,nrow=nn)

for(t in 1:nn){
    x         = vech2full(reaco[t,])
    RCov[,,t] = x[1:dm,1:dm]
    RCor[,,t] = nearcor(round(cov2cor(RCov[,,t]),8))$cor
    stand[t,] = rets[t,]/sqrt(diag(RCov[,,t]))
    RVs[t,]   = diag(RCov[,,t])
}

apply(stand,2,mean)
apply(stand,2,sd)

stand = stand-matrix(rep(apply(stand, 2, mean),nn),nrow=nn,byrow = TRUE)
stand = stand/matrix(rep(apply(stand, 2, sd),nn),nrow=nn,byrow = TRUE)

round(apply(stand,2,mean),4)
round(apply(stand,2,sd),4)


desc = matrix(NA,ncol=10,nrow=10)

for(i in 1:dm){
    desc[i,1] = mean(stand[,i])
    desc[i,2] = median(stand[,i])
    desc[i,3] = sd(stand[,i])
    desc[i,4] = skewness(stand[,i])
    desc[i,5] = kurtosis(stand[,i])
    t = shapiro.test(stand[,i])
    desc[i,6] = t$p.value
    t = ks.test(stand[,i],y='pnorm')
    desc[i,7] = t$p.value
    t = jarque.test(stand[,i])
    desc[i,8] = t$p.value
    t = Box.test (stand[,i], lag = 5, type = "Ljung")
    desc[i,9] = t$p.value
    t = Box.test (stand[,i], lag = 10, type = "Ljung")
    desc[i,10] = t$p.value
}

rownames(desc) = names
colnames(desc) = c('mean','median','sd','skew.','kurt.',
                   'SW','KS','JB','LB(5)','LB(10)')

round(desc,2)

dfs = rep(NA,dm)

for(i in 1:dm){
    model = fitdistr(stand[,i], "t")
    dfs[i]=model$estimate[3]
}




