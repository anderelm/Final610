## Anderson Computing Final Q1

# direct inversion
r_squared_fun_1 = function(x,y){
  y = y - mean(y)
  x = sweep(x,2,colMeans(x),"-")
  SSE_intercept_only = sum(y^2)
  xty = crossprod(x,y)
  xtx = crossprod(x)
  beta = solve(xtx)%*%xty
  SST = crossprod(xty,beta)[1,1]
  r_squared = SST/SSE_intercept_only
  return(r_squared)
}

# solve linear system
r_squared_fun_2 = function(x,y){
  y = y - mean(y)
  x = sweep(x,2,colMeans(x),"-")
  SSE_intercept_only = sum(y^2)
  xty = crossprod(x,y)
  xtx = crossprod(x)
  beta = solve(xtx,xty)
  SST = crossprod(xty,beta)[1,1]
  r_squared = SST/SSE_intercept_only
  return(r_squared)
}

# invert using eigen system
r_squared_fun_3 = function(x,y){
  y = y - mean(y)
  x = sweep(x,2,colMeans(x),"-")
  SSE_intercept_only = sum(y^2)
  xty = crossprod(x,y)
  xtx = crossprod(x)
  eigen_xtx = eigen(xtx,symmetric=TRUE)
  u = eigen_xtx$vectors
  d_inv = diag(1/eigen_xtx$values)
  beta = u %*% d_inv %*% t(u) %*% xty
  SST = crossprod(xty,beta)[1,1]
  r_squared = SST/SSE_intercept_only
  return(r_squared)
}

# avoid inversion using svd
r_squared_fun_4 = function(x,y){
  y = y - mean(y)
  x = sweep(x,2,colMeans(x),"-")
  SSE_intercept_only = sum(y^2)
  svd_x = svd(x)
  u = svd_x$u
  uty = crossprod(u,y)
  SST = crossprod(uty)[1,1]
  r_squared = SST/SSE_intercept_only
  return(r_squared)
}

# avoid inversion using svd
# ignore right singular vectors
r_squared_fun_5 = function(x,y){
  y = y - mean(y)
  x = sweep(x,2,colMeans(x),"-")
  SSE_intercept_only = sum(y^2)
  svd_x = svd(x,nv=0)
  u = svd_x$u
  uty = crossprod(u,y)
  SST = crossprod(uty)[1,1]
  r_squared = SST/SSE_intercept_only
  return(r_squared)
}

setseed(1234567890)

N=20
n.vec=seq(N)
r.sq_times <- matrix(, nrow=5, ncol=20)
rownames(r.sq_times) = c("r_squared_fun_1", "r_squared_fun_2", "r_squared_fun_3", "r_squared_fun_4", "r_squared_fun_5")

p=seq(1, 20, by=1)
for (n in 1:length(p)) {
  utils::Rprof(interval=0.01)
  x=matrix(rnorm(1000*500), nrow=1000,  ncol=500)
  y=matrix(rnorm(1000*1), nrow=1000,  ncol=1)
  r_squared_fun_1(x, y)
  r_squared_fun_2(x, y)
  r_squared_fun_3(x, y)
  r_squared_fun_4(x, y)
  r_squared_fun_5(x, y)
  Rprof(NULL)
  Rprof_summ = summaryRprof()$by.total
  r.sq_times[1,n] = Rprof_summ['"r_squared_fun_1"', 'total.time']
  r.sq_times[2,n] = Rprof_summ['"r_squared_fun_2"', 'total.time']
  r.sq_times[3,n] = Rprof_summ['"r_squared_fun_3"', 'total.time']
  r.sq_times[4,n] = Rprof_summ['"r_squared_fun_4"', 'total.time']
  r.sq_times[5,n] = Rprof_summ['"r_squared_fun_5"', 'total.time']
  
}

## HAVE TO MAKE THE 5X6 MATRIX
summary_matrix=apply(r.sq_times, 1, summary)
summary_matrix

##The function that takes the least amount of time on average is function 2. 
# This function is likely faster than the others because it only calculates 
# the values that are needed to calculate R^2. For example, function 5 calculates
# the entire singular value decomposition and then only uses the u. 
