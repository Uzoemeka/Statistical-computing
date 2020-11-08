set.seed(2)
y=rpois(20,0.2)
Lambda=seq(0.01,0.8,length.out = 100)
l=-20*Lambda+20*0.2*log(Lambda)
plot(Lambda,l,type = "l",xlab="parameters",ylab="LogL")
mm=Lambda[which.max(l)]
abline(v=0.2)
abline(v=mm,col="red")
lh0=max(l)-1.92
abline(h=lh0,col="red")
Lambda[lh0]
abline(v=0.062,col="blue")
abline(v=0.46484848,col="green")
l2=-20*0.1+20*0.2*log(0.1)
l2
1-pchisq(1.54, df=20, lower.tail=FALSE)
##########################################################################


#================= statistical computing =============================
#================= 2(a) 
set.seed(10)
n=2000;a=10;b=2;sig=30
t2=runif(n,20,100);e=rexp(n,1/sig)
y=a + b*t2 + e
dat=data.frame(y,t2)
hist(y)

#============== 2(b)i
#=== loglikelihood of the model
rm(list = ls())
min.RSS=function(theta,x){
  J = n*log(theta[3]) + (1/(2*theta[3]^2))*sum((y-theta[1]-theta[2]*t2-theta[3])^2) + (n/2)*log(2*pi)
  return(J) 
}

#=======loglikelihood of the model with respect to theta1
theta1=seq(2, 35, len = 10000) 
log.min=rep(NA, len=length(theta1)) 
for(i in 1:length(theta1)){
  log.min[i]=min.RSS(theta = c(theta1[i],2,30), x = dat)
}
plot(theta1, log.min, type='l',main="log-likelihood with respect to a",xlab="a",ylab="log-likelihood")
index.min=which(log.min == min(log.min))
theta1[index.min]#==initial value for a

#=======loglikelihood of the model with respect to theta2
theta2=seq(1, 5, len = 1000) 
log.min=rep(NA, len=length(theta2)) 
for(i in 1:length(theta2)){
  log.min[i]=min.RSS(theta = c(10,theta2[i],30), x = dat)
}
plot(theta2, log.min, type='l',main="log-likelihood with respect to b",xlab="b",ylab="log-likelihood")
index.min=which(log.min == min(log.min))
theta2[index.min]#==initial value for b

#=======loglikelihood of the model with respect to theta3
theta3=seq(18, 30, len = 10000) 
log.min=rep(NA, len=length(theta3)) 
for(i in 1:length(theta3)){
  log.min[i]=min.RSS(theta = c(10,2,theta3[i]), x = dat)
}
plot(theta3, log.min, type='l',main="log-likelihood with respect to sigma",xlab="sigma",ylab="log-likelihood")
index.min=which(log.min == min(log.min))
theta3[index.min]#==initial value for sigma


#======== estimation of the model parameters
theta=c(8.7,1.9,23)
result=optim(theta, fn= min.RSS, dat)
result
y1=result$par[1] + result$par[2]*t2 + result$par[3]

plot(t2,y,main="Regression line",xlab="predictor",ylab="response")
lines(t2,y1,type="l",col="red")


##################################################################
rm(list=ls())
B.RSS=function(theta,x){
  r=y-theta[1]-theta[2]*t2
  bs= -dnorm(y, mean=theta[1]-theta[2]*t2-theta[3], sd=sqrt(theta[3]), log=F)
  return(bs) 
}
#pnorm(q, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)
set.seed(11)
n=20;a=10;b=2;sig=30
t2=runif(n,20,100);e=rexp(n,1/sig)
y=a + b*t2 + e
dat=data.frame(y,t2)

theta=c(10,2,32)
result1=optim(theta, fn= B.RSS, dat)
result1
y2=result1$par[1] + result1$par[2]*t2 + result1$par[3]
lines(t2,y2,type="l",col="blue")


#===================2(c)==============
#===== for sig=1
rm(list = ls())
min.RSS=function(r,sig1){
  n=20
  J = n*log(sig1) + (1/(2*sig1^2))*sum((r-sig1)^2) + (n/2)*log(2*pi)
  return(J) 
}
rr=seq(-3, 3, len = 100) 
log.min=rep(NA, len=length(rr)) 
for(i in 1:length(rr)){
  log.min[i]=min.RSS(r=rr[i],sig1=1)
}
plot(rr, log.min, type='l')

#=======for r=1
sig.1=seq(1/5, 10, len = 100) 
log.min=rep(NA, len=length(sig.1)) 
for(i in 1:length(sig.1)){
  log.min[i]=min.RSS(r=1,sig1=sig.1[i])
}
plot(sig.1, log.min, type='l')



###
rm(list = ls())
B.RSS=function(r,sig1){
  n=20
  bs=mean(-dexp(r, rate = 1/sig1, log = FALSE) + 1/(4*sig1))
  return(bs) 
}
rr=seq(-3, 3, len = 100) 
log.min=rep(NA, len=length(rr)) 
for(i in 1:length(rr)){
  log.min[i]=B.RSS(r=rr[i],sig1=1)
}
plot(rr, log.min, type='l')



#########
gramschmidt <- function(x) {
  x <- as.matrix(x)
  # Get the number of rows and columns of the matrix
  n <- ncol(x)
  m <- nrow(x)
  
  # Initialize the Q and R matrices
  q <- matrix(0, m, n)
  r <- matrix(0, n, n)
  
  for (j in 1:n) {
    v = x[,j] # Step 1 of the Gram-Schmidt process v1 = a1
    # Skip the first column
    if (j > 1) {
      for (i in 1:(j-1)) {
        r[i,j] <- t(q[,i]) %*% x[,j] # Find the inner product (noted to be q^T a earlier)
        # Subtract the projection from v which causes v to become perpendicular to all columns of Q
        v <- v - r[i,j] * q[,i] 
      }      
    }
    # Find the L2 norm of the jth diagonal of R
    r[j,j] <- sqrt(sum(v^2))
    # The orthogonalized result is found and stored in the ith column of Q.
    q[,j] <- v / r[j,j]
  }
  
  # Collect the Q and R matrices into a list and return
  qrcomp <- list('Q'=q, 'R'=r)
  return(qrcomp)
}

#===============================================================================
#===============================================================================
rm(list=ls())
set.seed(10)
n=20;a=10;b=2;sig=30
t2=runif(n,20,100)
e=rexp(n,1/sig)
y=a + b*t2 + e


rss <- function(p,y,t2){
  a1=p[1];b1=p[2];log.s=p[3]
  rr=y - a1 - b1*t2 - log.s
  Log.p=-sum(dexp(rr,1/exp(log.s), log=TRUE))
  return(Log.p)
}

result1=optim(par=c(1,1,1), fn=rss, y=y, t2=t2)
y1=result1$par[1] + result1$par[2]*t2 + result1$par[3]

plot(t2,y,main="Regression line",xlab="predictor",ylab="response",ylim = c(60,219))
lines(t2,y1,type="l",col="red")


#===================2(c)==============
#===== for sig=1
rss1 = function(rr,s){
 J=-sum(dexp(rr,1/exp(s), log=TRUE))
  return(J)
}

rr=seq(-3, 3, len = 100) 
log.min1=rep(NA, len=length(rr)) 
for(i in 1:length(rr)){
  log.min1[i]=rss1(rr=rr[i],s=1)
}
plot(rr, log.min1, type='l')

#=======for r=1
sig=seq(1/5, 10, len = 100) 
log.min2=rep(NA, len=length(sig)) 
for(i in 1:length(sig)){
  log.min2[i]=rss1(r=1,s=sig[i])
}
plot(sig, log.min2, type='l')


#=================================================================================
Brier <- function(p,y,t2){
  a1=p[1];b1=p[2];log.s=p[3]
  r=y - a1 - b1*t2
  nrm=(norm(dexp(as.matrix(r),1/exp(log.s)),type="f"))^2
  Log.p = mean(-dexp(r,1/exp(log.s)) + 0.5*nrm )
  return(Log.p)
 }

result2=optim(par=c(10,2,30), fn=Brier, y=y, t2=t2)
y2=result2$par[1] + result2$par[2]*t2 + result2$par[3]

plot(t2,y,main="Regression line",xlab="predictor",ylab="response",ylim = c(60,219))
lines(t2,y2,type="l",col="blue")

#============================================================================
Brier1 = function(rr,s){
  nrm=(norm(dexp(as.matrix(rr),1/exp(s)),type="f"))^2
  B= mean(-dexp(rr,1/exp(s)) + 0.5*nrm )
  return(B)
}

rr=seq(-3, 3, len = 100) 
log.min3=rep(NA, len=length(rr)) 
for(i in 1:length(rr)){
  log.min3[i]=Brier1(rr=rr[i],s=1)
}
plot(rr, log.min3, type='l')


#=======for r=1
sig=seq(1/5, 10, len = 100) 
log.min4=rep(NA, len=length(sig)) 
for(i in 1:length(sig)){
  log.min4[i]=Brier1(r=1,s=sig[i])
}
plot(sig, log.min4, type='l')


# r=as.matrix(y - a - b*t2)
# nrm=(norm(dexp(r,1/exp(30)),type="f"))^2
# 
# bnrm <- function(p,y,t2){
#   a1=p[1];b1=p[2];log.s=p[3]
#   bnrm=(norm(dexp(as.matrix(y - a1 - b1*t2),1/exp(log.s)),type="f"))^2
#   return(bnrm)
# }
# nrm=optim(par=c(1,1,1), fn=bnrm, y=y, t2=t2)


# s <- 2
# x <- 2
# dexp(x,s)
# s*exp(-s*x)


# fcn=function(p,y,t2){
#   a1=p[1];b1=p[2];log.s=p[3]
#   e=y - a1 - b1*t2
#   # Log.p= n*exp(log.s) + (1/exp(log.s))*1/sum(e)
#   Log.p <- -sum(log(  (1*exp(log.s))  * exp(-e*exp(log.s))  )  )
#   return(Log.p)
# }
# 
# e
# phi <- 1/2
# dexp(e,1/phi, log=FALSE)
# (1/phi) * exp(-e/phi)
# 
# 
# 
# result=optim(par=c(10,2,30), fn=fcn, y=y, t2=t2)
# y1=result$par[1] + result$par[2]*t2 + result$par[3]
# 
# plot(t2,y,main="Regression line",xlab="predictor",ylab="response",ylim = c(50,240))
# lines(t2,y1,type="l",col="red")

# rss <- function(p,y,t2){
#   a1=p[1];b1=p[2];sig=p[3]
#   e=y - a1 - b1*t2
#   Log.p=-sum(log(dexp(e,sig)))
#   return(Log.p)
# }


# gd=function(p,y,t2){# Gradient of 'fn'
#     a1=p[1];b1=p[2];s=p[3]
#     c(-n/exp(s) + sum(y - a1 - b1*t2)/exp(s^2), n/exp(s), sum(t2)/exp(s))
# }
# optim(par=c(1,1,1), fcn, gd, y,t2, method = "CG")


#====================================================================
xx=c(1,2,3,4,5)
zz=c(3,4,2,2,1)
yy=c(30,40,22,33,40)

funk=function(param,x,y,z){
  a=rep(param[1],5)
  b=param[2]
  d=param[3]
  fit=sum((y-(a+b*x+z*d))^2)
  return(fit)
}