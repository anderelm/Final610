library(SuppDists)
library(mvtnorm)

draw_tau = function(data,state){
  shape = 1 + length(data$y)
  rate = 1 + sum(abs(data$y-data$x%*%state$beta))
  return(rgamma(1,shape,rate))
}

draw_beta = function(data, state){
  browser()
  z= rinvGauss(n=length(data$y), 
               nu=(1/(state$tau*abs(data$y-t(data$x)%*%state$beta))), 
               lambda=1)

  I=diag(1, nrow=ncol(data$x), ncol=ncol(data$x))
  Z=diag(z, nrow=length(z), ncol=length(z))
  mean=solve(((state$tau)^2)%*%(t(data$x)%*%Z%*%data$x)+I)%*%((state$tau^2)%*%t(data$x)%*%Z%*%data$y)
  sigma=solve((state$tau^2)%*%t(data$x)%*%Z%*%data$x + I)
  return(beta=mvtnorm::rmvnorm(n=length(data$y), mean=mean, sigma=sigma ))

}

update_current_state = function(data, state){
  state$beta = draw_beta(data,state)
  state$tau = draw_tau(data,state)
  return(state)
}
