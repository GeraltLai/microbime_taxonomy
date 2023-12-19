####################
##### function #####
####################
### manova
manova=function(data,meta,factor){
  ##### distance matrix 
  matrix=diag(x=1,nrow=nrow(data))
  for(i in 1:(nrow(data)-1))
    for(j in (i+1):nrow(data)){
      total=lapply(1:length(data[i,]), function(x) min(data[i,x],data[j,x]))
      total=do.call(rbind,total)
      s1=sum(data[i,])
      s2=sum(data[j,])
      ans=1-(2*sum(total)/(s1+s2))
      matrix[i,j]=ans
    }
  label=t(matrix)[lower.tri(matrix)]
  # After transposing, remove the triangle and fill the lower triangle
  matrix[lower.tri(matrix)]=t(matrix)[lower.tri(matrix)]
  ##### manova
  # total number of observations (sst)
  N=nrow(data)
  sst=sum(label^2)/N
  sst
  # within-group or residual sum of squares (ssw)
  a=length(unique(meta[[factor]]))
  k=0
  adj_dist=label
  number_factor = table(meta[[factor]])
  for(i in 1:(nrow(data)-1))
    for(j in (i+1):nrow(data)){
      k=k+1
      if(meta[i,factor]==meta[j,factor]){
        n = number_factor[meta[i,factor]]
        adj_dist[k]=label[k]/sqrt(n)
      }else{
        adj_dist[k]=0
      }
    }
  ssw=sum(adj_dist^2)
  ssw
  ### ssa 
  ssa=sst-ssw
  ssa
  ### F value
  f=(ssa/(a-1))/(ssw/(N-a))
  return(data.frame(ssa,ssw,sst,f))
}

### permutation to calculate p-value
permutation=function(f_value,data,meta,factor,permutation=999){
  count=0
  for(i in 1:permutation){
    new_meta=meta
    x=new_meta[[factor]]
    y=sample(x,length(x))
    new_meta[[factor]]=y
    ans=manova(data=data,meta=new_meta,factor=factor)
    f_per=ans[4]
    if(f_per>f_value){
      count=count+1
    }
  }
  return((count+1)/(permutation+1))
}