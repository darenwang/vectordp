rm(list = ls())

#####functions#######
source("functions.R")
source("dp_functions.R")
source("localscreen.R")
 
#####packages#######
library(glmnet)
library(gglasso)
 
set.seed(1)
 
#START define parameters
#rho can not be too large so that the time series is stable
#candidate of rho
 rho.can=c( 0.4)
 
 ll=length(rho.can)
 
#p is dimension 


p=20

#sigma is the variance of noise
eps_sigma=1

 
n=30
 
# minimal sample size for lasso and group lasso estiamtors; needed for any lasso estimators
delta1=10
delta.local=10
 
#RR is the number of repetition
RR=100
 
 

#parameters
v1= 2*(seq(1,p,1)%%2)-1
 
 v2=-v1
AA=matrix(0,nrow = p,ncol=p-2)


#records of change point estimation
dp.rec =  vector("list", length = ll)
lr.rec = vector("list", length = ll)
 
 
#checking the errors of DP and LR 
haus1=rep(0,ll)
haus2=rep(0,ll)
#
for( candidate in 1:ll){
  
  #generate transition matrices
  A1=cbind(v1,v2,AA)
  A2=cbind(v2,v1,AA)
  A3=A1
  
  print(paste("candidate=" ,candidate))
  rho=rho.can[candidate]
  A1=A1*rho
  A2=A2*rho
  A3=A3*rho
   

  
   
  
  
  
  
   
  
  for( rr in 1:RR){
    print(paste("rr=" ,rr))
    data=gen.var.data (eps_sigma,p,A1, 2*n+1,vzero=c( ))
    data=cbind(data,gen.var.data (eps_sigma,p,A2, 2*n,vzero=c(data[,ncol(data)])))
    data=cbind(data,gen.var.data (eps_sigma,p,A3, 2*n,vzero=c(data[,ncol(data)])))
    
    #For convenient due to data splitting, col(data) (or the length of the time series) is odd.
    #the change points are at 2*n and 4*n
    
    lambda.lasso.list=c(  3,4 )
    
    gam.list=seq(1, 50,2) 
    #START: DP
     start.time <- Sys.time()
     dp.estimate=dp.main.function (data, p, lambda.lasso.list ,gam.list, delta1)
     end.time <- Sys.time();time.taken <- end.time - start.time;print(time.taken)
    
     print(  "dp=")
     print(dp.estimate)
     dp.rec[[candidate]][[rr]]=list(dp.estimate[-c(1,length(dp.estimate))])
     haus1[candidate]=haus1[candidate]+hausdorff.distance(dp.estimate,c(0, 2*n,4*n,6*n) )
     print(haus1/rr)
    #END:  DP
    
    
    
     lambda.group.list =c( 0.5,1, 1.5 )
     
    #START  local refinement#include X and Y
     lr.estimate= local.main.function  (data,delta.local, lambda.group.list,dp.estimate)
     
     print( "lr=")
       print(lr.estimate)
       lr.rec[[candidate]][[rr]]=list(lr.estimate[-c(1,length(lr.estimate))])
       haus2[candidate]=haus2[candidate]+hausdorff.distance(lr.estimate,c(0, 2*n,4*n,6*n) )
       print(haus2/rr)
    #END local refinement
    
   
    
  }
}



 



 


