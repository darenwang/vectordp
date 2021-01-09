#generate var data

gen.var.data= function(eps_sigma,p,A,n,vzero){
X=matrix(0,nrow=p,ncol=n)
 
if(length(vzero)!=p){
X[,1]= rnorm(p, mean=0, sd= eps_sigma)
}else{
  X[,1]=vzero
}
for ( t in 2:n){
  
  
    X[,t]=rnorm(p,mean=A%*%X [,t-1],sd=eps_sigma)

  
}
return(X)}

#compuate residuals on given interval 
 

 

  



#hausdorff.distance
hausdorff.distance=function(v1,v2){
  p1=length(v1)
  p2=length(v2)
  distance.mat=matrix(0,nrow=p1,ncol=p2)
  for( i in 1: p1){
    for( j in 1: p2){
      distance.mat[i,j]=abs(v1[i]-v2[j])
      
    }
  }
  dis=max(max( apply(distance.mat, 1, min)),max( apply(distance.mat, 2, min)))
  return(dis)
  
}


colSD=function(matrix){
  Ncol=ncol(matrix)
  result=sapply(c(1:Ncol), function(t) sd(matrix[,t]))
  return(result)  
  
}

