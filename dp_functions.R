
#compute res given interval [s,e]

residualvar=function(X, Y,s,e,p,lambda){
  options(warn=-1)
  estimate=matrix(0,nrow=p,ncol=p)
  
  
  for ( m in 1:p){
    out=glmnet(x=t(X[,  s  : e ]), y=Y[m,  s  :e], family=c("gaussian"),
               alpha = 1, lambda=lambda/sqrt(e-s) )#,intercept=F)
    estimate[m,]= 
      as.vector(out$beta )
    
    
  }
   norm(estimate-A1,type="F" )
   norm(A1,type="F")
  y.estimate= estimate%*%X[, s:e ]
  return(norm(y.estimate -Y[, s:e ], type="F" )^2)
}
###end of function



# pull the change point given partitions
dp.pull.change= function(N,  partition){
  cc =N
  estimate.change=cc
  while(cc[1]!=1){
    estimate.change=c(partition[cc],estimate.change)
    cc=partition[cc]
  }
  return(estimate.change[-length(estimate.change)]-1)
  
} 
###end of function

#dp matrix  function
dp.matrix.function=function(X.train, Y.train, p, lambda.lasso , delta1){
  N=ncol(X.train)
dp.matrix=matrix(Inf, N, N) 
for (s in 1:N){
  for (e in 1:N){
    if ( e>s+delta1  ){
      #print(c(i,j))
      dp.matrix[s,e]=residualvar(X=X.train, Y=Y.train,s ,e ,p,lambda.lasso)
    }
    
    
  }
  
}  
return(dp.matrix)
}

#end 

#compute res plus gam
res.and.gam=function(r,l ,gam,Best.value,dp.matrix   ){
  
  result=  ifelse(l==1, -gam, Best.value[l] ) +gam+   dp.matrix[l,r] 
  return(result)
}

###end of function


dp.gam.function=function(dp.matrix ,gam, delta1){
  N=ncol(dp.matrix)
 Best.value=rep(Inf, N)
partition= rep(1,N ) 
 

 
for( r in 2:N){
   
  b=sapply(1:r,function(l) 
    res.and.gam(r,l, gam,Best.value,dp.matrix  )     
  )
  if (min(b) <Best.value[r]){Best.value[r]=min(b)
  partition[r]=which.min(b)} 
  
  
}


dp.estimate=dp.pull.change(N, partition ) #pull out change points
 
return(c( dp.estimate ,N))
}

dp.gam.list=function(X.train, Y.train, p, lambda.lasso ,gam.list, delta1){
  
  dp.matrix =dp.matrix.function (X.train, Y.train, p, lambda.lasso , delta1) 
  dp.list=vector('list',length= length(gam.list))
  for ( i in 1: length(gam.list)){ 
  dp.list[[i]]=   dp.gam.function(dp.matrix ,gam.list[i], delta1)  
  }
  return(dp.list)
  
}

compute.test.errors=function(X.train, Y.train,  X.test,Y.test,lambda.lasso,estimate.changes){
  res=0
  for ( l in 2:length(estimate.changes)){
    
    s=estimate.changes[l-1]+1
    e=estimate.changes[l]
    estimate=matrix(0,nrow=p,ncol=p)
    for( m in 1:p){
    out=glmnet(x=t(X.train[,  s  : e ]), y=Y.train[m,  s  :e], family=c("gaussian"),
               alpha = 1, lambda=lambda.lasso/sqrt(e-s )) #,intercept=F)
    estimate[m,]= 
      as.vector(out$beta ) }
    res=res+ norm( Y.test[,s  :e] -estimate%*%X.test[, s:e ], type="F" )^2
   
  }  
  return(res)
  
}

dp.main=function(X.train, Y.train,X.test,Y.test, p, lambda.lasso.list ,gam.list, delta1){
  rec.temp=rep(0, nrow=length(lambda.lasso.list) )
  dp.lambda.list=vector('list', length=length(lambda.lasso.list))
  for (ll in 1:length(lambda.lasso.list)){
    #print(ll)
    dp.list= dp.gam.list (X.train, Y.train, p, lambda.lasso.list[ll] ,gam.list, delta1)  
    res.vec.temp= sapply(1:length(gam.list), function(j)
      compute.test.errors (X.train, Y.train,  X.test,Y.test,lambda.lasso.list[[ll]],dp.list[[j]]) )
    rec.temp[ll]=min(res.vec.temp)       
    dp.lambda.list[[ll]]= dp.list[[which.min(res.vec.temp)]]
  }
  
  
  return(dp.lambda.list[[which.min( rec.temp)]])
  
}

#end of dp main


#use in the simulations
dp.main.function=function(data, p, lambda.lasso.list ,gam.list, delta1){
 data.temp=data
   if (ncol(data)%%2 ==0){
    data.temp=data[,2:ncol(data)]
  }
  
  
  X=data.temp[,1:(ncol(data.temp)-1)]
  Y=data.temp[,2:ncol(data.temp)]
  X.train= X[,seq(1,ncol( X),2)]
  X.test= X[,seq(1,ncol( X)-1,2)+1]
  Y.train= Y[,seq(1,ncol( Y),2)]
  Y.test= Y[,seq(1,ncol( Y)-1,2)+1]
  return(dp.main(X.train, Y.train,X.test,Y.test, p, lambda.lasso.list ,gam.list, delta1))
  
}