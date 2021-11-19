rm(list=ls(all=TRUE))

# install.packages("remotes")
# remotes::install_github("jonathancornelissen/highfrequency")
library(highfrequency)
library("data.table")
library(readxl)

rcKernel <- rKernelCov(rData = sampleOneMinuteData, makeReturns = TRUE)
rcKernel

assets  = c( "SP500","FTSE100","SPLAC","NIKKEI","HSI")
dm      = length(assets)
df      = NULL
df$Date = seq(as.Date("1975/1/1"), as.Date(Sys.Date()), by = "day")

freq= 'month'

for (i in 1:dm){
    data  = read_excel(paste('asset_data/',assets[i],'.xlsx',sep=''))
    tmp   = data$PX_LAST
    Date  = as.Date(data$Date)
    tmpdf = data.frame(tmp,Date)
    colnames(tmpdf) = c(assets[i],"Date")
    df    = merge.data.frame(df,tmpdf,by='Date',all=TRUE)
}

df    = na.omit(df)
df$DT = as.Date(cut(df$Date, freq))
df$Date = NULL

setDT(df)
class(df)

df

rcKernel <- rKernelCov(rData = df, makeReturns = TRUE,kernelParam =0.9,
                       cor=TRUE)

indw      = unique(df$DT)



rck = array(as.numeric(unlist(rcKernel)), dim=c(dm, dm, length(rcKernel)))
rck = rck[,,-1]

load('data_month.Rdata')


C = combinat::combn(1:dm,2)
par(mfrow=c(2,ceiling(ncol(C)/2)))
for(i in 1:ncol(C)){
    plot(data$RCor[C[1,i],C[2,i],],type='l',ylim=c(-1.2,1.2),
         main=c(assets[C[1,i]],assets[C[2,i]]),ylab='',
         xlab='')
    lines(rck[C[1,i],C[2,i],],col=2,lwd=1)
}

data$rck = rck


save(data,file='data_month.Rdata')
