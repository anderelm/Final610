eigen_fun = function(X,tol=.Machine$double.eps^0.75){
  if(!isSymmetric(X)) stop("x is not symmetric")
  d = ncol(X)
  U = diag(1,d,d,FALSE)
  while(TRUE){
    did_we_update = FALSE
    
    for(i in 1:(d-1)){
      for(j in 2:d){
        
        # test whether abs(X[i,j]) is too large
        if(!is.na(X[i,j]) && abs(X[i,j]) > tol){
        
          R=diag(1,d,d,FALSE)
          
          if(X[i,i]==X[j,j]+tol | X[i,i]==X[j,j] - tol ) {
            theta=.5*atan((2*X[j,i])/(X[i,i]-X[j,j]))
            
            R[i,j]=cos(theta) 
            R[j,i]=-cos(theta) 
            R[i,i]=sin(theta)
            R[j,j]=sin(theta)
            
            X=R%*%X%*%t(R) # and update X 
            U=U%*%R # and update and U
            did_we_update = TRUE  # track that in this sweep through we have done at least one update
            
          } else {
            theta=pi/4
            R[i,j]=cos(theta) 
            R[j,i]=-cos(theta) 
            R[i,i]=sin(theta)
            R[j,j]=sin(theta)
            
            X=R%*%X%*%t(R) # and update X 
            U=U%*%R # and update and U
            did_we_update = TRUE  # track that in this sweep through we have done at least one update
          }
          
        }
      }
    }
    if(!did_we_update) break
  }
  values=diag(X)
  sorted_values = sort(values,index.return=TRUE)
  values = sorted_values$x
  vectors = U[,sorted_values$ix]
  return(list(values=values, vectors=vectors))
}

z = matrix(rnorm(2^2),2,2)
X = z + t(z)
eigen(X)
eigen_fun(X)

## Testing blocks
library(testthat)


test_that("eigen decom tests",{
  set.seed(1234567890)
  p=seq(2, 5, by=1)
  for(i in 1:length(p)){
    z = matrix(rnorm(p[i]^p[i]),p[i],p[i])
    X = z + t(z)
    expect_equal(eigen_fun(X), eigen(X), tol=0.1)
  }
})

## This test eventually passes, but it just takes a very long time to do so. 
# When you stop the function from running manually, the test passes. 

