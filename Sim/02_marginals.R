rm(list=ls(all=TRUE))
load('data.Rdata')

nn      = dim(data$ret)[1]
dm      = dim(data$ret)[2]
oos     = 5*12 # last 5 years

##----
## Estimate, predict and save output: ARMA-GARCH
##----

spec <- ugarchspec(variance.model = list(model = "sGARCH",
                                         garchOrder = c(1, 1),
                                         submodel = NULL,
                                         external.regressors = NULL,
                                         variance.targeting = FALSE),

                   mean.model     = list(armaOrder = c(1, 1),
                                         external.regressors = NULL,
                                         include.mean=TRUE),
                   distribution.model = "norm")

resid_arma_garch = u_arma_garch = hs_arma_garch=matrix(NA,ncol=dm,nrow=nn)
res_dist_arma_garch = list()

for(i in 1:dm){
    garch <- ugarchfit(spec = spec, data = data$ret[,i],
                       solver.control = list(trace=0))
    resid_arma_garch[,i]    = garch@fit$z
    hs_arma_garch[,i]    = garch@fit$sigma
    res_dist_arma_garch     = c(res_dist_arma_garch,list(edfun(resid_arma_garch[,i])))
    u_arma_garch[,i]        = (res_dist_arma_garch[[i]]$pfun(resid_arma_garch[,i]))*nn/(nn+1)
}

data = append(data,list(u_arma_garch=u_arma_garch))

par(mfrow=c(2,ceiling(dm/2)))
for(i in 1:dm){
    plot(data$date,hs_arma_garch[,i]^2,type='l',
         main=data$assets[i],ylab='',xlab='',lwd=2)
    lines(data$date,data$RV[,i],col=2)
}

##----
## Estimate, predict and save output: ARMA-RV
##----

# obtain u~U(0,1)data

resid_arma_rv = u_arma_rv = hs_arma_rv=matrix(NA,ncol=dm,nrow=nn)
res_dist_arma_rv = list()

for(i in 1:dm){
    e     = data$ret[,i]/sqrt(data$RV[,i])
    model = arima0(e,order=c(1,0,1))
    sigma = sd(model$residuals)
    resid_arma_rv[,i] =model$residuals/sigma
    res_dist_arma_rv = c(res_dist_arma_rv,list(edfun(resid_arma_rv[,i])))
    u_arma_rv[,i] = (res_dist_arma_rv[[i]]$pfun(resid_arma_rv[,i]))*nn/(nn+1)
}

par(mfrow=c(2,ceiling(dm/2)))
for(i in 1:dm){
    plot(density(resid_arma_garch[,i]),lwd=2)
    lines(density(resid_arma_rv[,i]),col=2)
}

data = append(data,list(u_arma_rv=u_arma_rv))

save(data,file = 'data_u.Rdata')

# Dynamics of RV via HAR-type model with X1, X3 and X6


RVs  = matrix(NA,ncol=3,nrow=nn)
lags = c(1,3,6)
lrv  = log(data$RV[,1])
RVs[1,] = mean(lrv)

for(j in 1:3){
    for(t in 1:(nn-1)){
        RVs[(t+1),j] = (1/length(lrv[max(0,(t-lags[j]+1)):t]))*
            sum(lrv[max(0,(t-lags[j]+1)):t])
    }
}

colnames(RVs) =lags
RVs=data.frame(RVs)

m1=lm(lrv~X1+X3+X6,data=RVs)
se <- sqrt(diag(vcov(m1)))
m1$coefficients
m1$coefficients/se

