#one change point


convert.design.matrix.one.change=function(X,t){
  ee=ncol(X)
  xx1=t(X)
  xx1[ (t+1):ee,]=0
  
  xx2=t(X)
  xx2[1:t,]=0

  xx=cbind(xx1/sqrt(ee),xx2/sqrt(ee))
  index=c()
  for ( pp in 1:p){
    index=c(index,pp,pp+p)
  }
  xxx=xx[,index]
  return(xxx)
  
  }

train.res.group.lasso =function(t,y, X   ,lambda.group){
  group=rep(1:p,each=2)
   out=gglasso(x=convert.design.matrix.one.change(X,t),y=y,group=group, loss="ls",
              lambda=lambda.group/ncol(X) ,intercept = FALSE,eps = 0.0001)
  alpha=as.vector(out$beta)
  #test
  #beta=c(alpha[1:p ]/sqrt(t),alpha[(p+1):(2*p)]/sqrt(ncol(X)-t))
  #beta
  res=sum(( y - convert.design.matrix.one.change(X,t) %*%alpha  )^2)
  return(res)
}

test.res.group.lasso =function(t,y, X, y.test, X.test   ,lambda.group){
 
  group=rep(1:p,each=2)
   out=gglasso(x=convert.design.matrix.one.change(X,t),y=y,group=group, loss="ls",
              lambda=lambda.group/ ncol(X) ,intercept = FALSE,eps = 0.001)
   #test
  alpha=as.vector(out$beta)
  beta=c(alpha[1:p ]/sqrt(t),alpha[(p+1):(2*p)]/sqrt(ncol(X)-t))
  beta
  res=sum(( y.test - convert.design.matrix.one.change(X.test,t) %*%alpha  )^2)
  return(res)
}



res.at.t=function(s,e,t,X.train,Y.train,lambda.group ){
return(sum(sapply( 1:p, function(m) 
  train.res.group.lasso(t , y=Y.train[m,   s :e], X = X.train[, s: e ]  ,lambda.group ) )))
}


find.one.change.grouplasso=function(s,e,X.train,Y.train,delta.local   ,lambda.group ){
  
  estimate= (s+e)/2
 
  if( e-s > 2*delta.local){
  can.vec=c((s+delta.local): (e-delta.local))
  can.vec=can.vec[which(can.vec%%  2==0)]
   res.seq=sapply(can.vec,function(t) 
    res.at.t(s,e,t-s+2  ,X.train,Y.train,lambda.group)) 
  #plot(can.vec,res.seq)
  estimate= can.vec[which.min(res.seq)]
  #lv= min(res.seq)
  }

    return(  estimate )
  
}


group.lasso.test.error=function(s,e,X.train, Y.train, X.test, Y.test, delta.local, lambda.group.list ){
  estimate.temp=rep(0,length(lambda.group.list) )
  test.errors.temp=rep(0,length(lambda.group.list) )
  for ( ll in 1:length(lambda.group.list)){
     
  estimate.temp[ll]=find.one.change.grouplasso(s,e,X.train,Y.train,delta.local,lambda.group.list[ll] )
  test.errors.temp[ll]=sum(sapply( 1:p, function(m) 
    test.res.group.lasso( estimate.temp[ll]-s+2 ,y=Y.train[m,s:e],X= X.train[, s: e  ], 
                          y.test=Y.test[m, s:e], X.test= X.test[, s: e  ],lambda.group.list[ll] ) ))
  #print(estimate.temp[ll])
  }

  return(lambda.group.list[which.min(test.errors.temp)])
}

local.group=function(X, Y,X.train, Y.train, X.test, Y.test,delta.local, lambda.group.list,estimated.changes ){
  result.changes= estimated.changes
  
  KK=length(estimated.changes)-2
  if (KK>0){  
    result.changes[1]=1
  for ( kk in 1:KK){
    #print(kk)
    lambda.temp=group.lasso.test.error (s=result.changes[kk],e=result.changes[kk+2],X.train, Y.train, X.test, Y.test, delta.local, lambda.group.list )
    temp.estimate=  
      find.one.change.grouplasso(s=2*result.changes[kk],e=2*result.changes[kk+2],X ,Y ,delta.local   ,lambda.group= lambda.temp )
    result.changes[kk+1]=round(temp.estimate/2)
    #print(temp.estimate)
  }
    }
  result.changes[1]=0
  return(result.changes)
 
}


local.main.function= function (data,delta.local, lambda.group.list,dp.estimate){
  data.temp=data
  if (ncol(data)%%2 ==0){
    data.temp=data[,2:ncol(data)]
  }
  
    temp.estimate=dp.estimate/2

  X=data.temp[,1:(ncol(data.temp)-1)]
  Y=data.temp[,2:ncol(data.temp)]
  X.train= X[,seq(1,ncol( X),2)]
  X.test= X[,seq(1,ncol( X)-1,2)+1]
  Y.train= Y[,seq(1,ncol( Y),2)]
  Y.test= Y[,seq(1,ncol( Y)-1,2)+1]
  return(2* local.group (X, Y,X.train, Y.train, X.test, Y.test,delta.local, lambda.group.list,temp.estimate))
}


 
