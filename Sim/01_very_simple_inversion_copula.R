rm(list=ls(all=TRUE))

set.seed(123)
nn     = 1000
dim    = 3
mu.tr  = round(rnorm(dim)*10)
Sig.tr = round(rWishart(1,10,diag(dim)))[,,1]
R.tr   = diag((diag(Sig.tr)^(-1/2)))%*%Sig.tr%*%diag((diag(Sig.tr)^(-1/2)))

mu.tr
Sig.tr
R.tr
cov2cor(Sig.tr)

# pseudo-data

x      = matrix(rnorm(nn*dim),ncol=dim)%*%chol(Sig.tr)+matrix(mu.tr,nrow=nn,ncol=dim,byrow=T)
apply(x,2,mean)
mu.tr
cov(x)
Sig.tr

# u-data

u      = matrix(NA,nrow=nn,ncol=dim)
for (i in 1:dim){
    u[,i] = pnorm(x[,i],mu.tr[i],sqrt(diag(Sig.tr))[i])
}

par(mfrow=c(1,dim))
for(i in 1:dim) {
    hist(u[,i],freq=F)
    abline(h=1)}

cor(u)
R.tr

# y-data
library(dst)




