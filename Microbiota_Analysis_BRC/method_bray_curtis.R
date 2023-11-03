rm(list = ls())

a = c(4,1,5,10)
b = c(2,1,3,6)


a = c(100, 100, 140, 350, 275, 325)
b = c(100, 200, 350, 400, 200, 300)

a = c(7, 6 ,10, 300, 500, 250)
b = c(9, 4 ,8, 11, 13, 1)
##### basic
sapply(1:length(a), function(x) abs(a[x]-b[x])/sum(a,b))
sum(sapply(1:length(a), function(x) abs(a[x]-b[x])/sum(a,b)))
basic = sum(sapply(1:length(a), function(x) abs(a[x]-b[x])/sum(a,b)))

##### method 1
weight = sapply(1:length(a), function(x) 1/length(a))
sapply(1:length(a), function(x) weight[x] * abs(a[x]-b[x])/(a[x]+b[x]))
sum(sapply(1:length(a), function(x) weight[x] * abs(a[x]-b[x])/(a[x]+b[x])))
method_1 = sum(sapply(1:length(a), function(x) weight[x] * abs(a[x]-b[x])/(a[x]+b[x])))

##### Standard Bray-Curtis(4)
weight = sapply(1:length(a), function(x) (a[x]+b[x])/sum(a,b))
sum(weight)
sapply(1:length(a), function(x) weight[x] * abs(a[x]-b[x])/(a[x]+b[x]))
sum(sapply(1:length(a), function(x) weight[x] * abs(a[x]-b[x])/(a[x]+b[x])))
method_2 = sum(sapply(1:length(a), function(x) weight[x] * abs(a[x]-b[x])/(a[x]+b[x])))

##### Proposed method (9)
weight = sapply(1:length(a), function(x) 1/((a[x] + b[x])/sum(a,b)))
top = sapply(1:length(a), function(x) weight[x] * abs(a[x]-b[x]))
down = sapply(1:length(a), function(x) weight[x] * (a[x] + b[x]))
method_3 = sum(top)/sum(down)

##### Proposed method (10)
reciprocal = sapply(1:length(a), function(x) 1/((a[x]+b[x])/sum(a,b)))
weight = sapply(1:length(a), function(x) reciprocal[x]/sum(reciprocal))
sapply(1:length(a), function(x) weight[x] * abs(a[x]-b[x])/(a[x]+b[x]))
sum(sapply(1:length(a), function(x) weight[x] * abs(a[x]-b[x])/(a[x]+b[x])))
method_4 = sum(sapply(1:length(a), function(x) weight[x] * abs(a[x]-b[x])/(a[x]+b[x])))

#### Compare 5 method
data.frame(basic, method_1, method_2, method_3, method_4)

### basic
weighted_bray = function(data){
  matrix = diag(x = 1, nrow = nrow(data))
  for(ii in 1:(nrow(data)-1))
    for(jj in (ii+1):nrow(data)){
      A = data[ii,] ; B = data[jj,]
      value = sapply(1:length(A), function(x) abs(A[x]-B[x])/sum(A,B))
      value = do.call(cbind, value)
      matrix[ii,jj] = sum(value)
    }
  label = t(matrix)[lower.tri(matrix)]
  # After transposing, remove the triangle and fill the lower triangle
  matrix[lower.tri(matrix)]=t(matrix)[lower.tri(matrix)]
  return(matrix)
}

##### method 1
weighted_bray = function(data){
  matrix = diag(x = 1, nrow = nrow(data))
  for(ii in 1:(nrow(data)-1))
    for(jj in (ii+1):nrow(data)){
      A = data[ii,] ; B = data[jj,]
      weight = sapply(1:length(A), function(x) 1/length(A))
      value = sapply(1:length(A), function(x) weight[x] * abs(A[x]-B[x])/(A[x]+B[x]))
      matrix[ii,jj] = sum(value)
    }
  label = t(matrix)[lower.tri(matrix)]
  # After transposing, remove the triangle and fill the lower triangle
  matrix[lower.tri(matrix)]=t(matrix)[lower.tri(matrix)]
  return(matrix)
}

##### method 2
weighted_bray = function(data){
  matrix = diag(x = 1, nrow = nrow(data))
  for(ii in 1:(nrow(data)-1))
    for(jj in (ii+1):nrow(data)){
      A = data[ii,] ; B = data[jj,]
      weight = sapply(1:length(A), function(x) (A[x]+B[x])/sum(A,B))
      value = sapply(1:length(A), function(x) weight[x] * abs(A[x]-B[x])/(A[x]+B[x]))
      matrix[ii,jj] = sum(value)
    }
  label = t(matrix)[lower.tri(matrix)]
  # After transposing, remove the triangle and fill the lower triangle
  matrix[lower.tri(matrix)]=t(matrix)[lower.tri(matrix)]
  return(matrix)
}

##### method 3
weighted_bray = function(data){
  matrix = diag(x = 1, nrow = nrow(data))
  for(ii in 1:(nrow(data)-1))
    for(jj in (ii+1):nrow(data)){
      A = data[ii,] ; B = data[jj,]
      weight = sapply(1:length(A), function(x) 1/((A[x] + B[x])/sum(A,B)))
      top = sapply(1:length(A), function(x) weight[x] * abs(A[x]-B[x]))
      down = sapply(1:length(A), function(x) weight[x] * (A[x] + B[x]))
      top = do.call(cbind, top) ; down = do.call(cbind, down)
      matrix[ii,jj] = sum(top)/sum(down)
    }
  label = t(matrix)[lower.tri(matrix)]
  # After transposing, remove the triangle and fill the lower triangle
  matrix[lower.tri(matrix)]=t(matrix)[lower.tri(matrix)]
  return(matrix)
}

##### method 4
weighted_bray = function(data){
  matrix = diag(x = 1, nrow = nrow(data))
  for(ii in 1:(nrow(data)-1))
    for(jj in (ii+1):nrow(data)){
      A = data[ii,] ; B = data[jj,]
      reciprocal = sapply(1:length(A), function(x) 1/((A[x]+B[x])/sum(A,B)))
      reciprocal = do.call(cbind, reciprocal)
      weight = sapply(1:length(A), function(x) reciprocal[x]/sum(reciprocal))
      value = sapply(1:length(A), function(x) weight[x] * abs(A[x]-B[x])/(A[x]+B[x]))
      value = do.call(cbind, value)
      matrix[ii,jj] = sum(value)
    }
  label = t(matrix)[lower.tri(matrix)]
  # After transposing, remove the triangle and fill the lower triangle
  matrix[lower.tri(matrix)]=t(matrix)[lower.tri(matrix)]
  return(matrix)
}


