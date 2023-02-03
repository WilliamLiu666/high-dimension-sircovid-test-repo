library(MASS)
library(coda)
HMC = function (U, grad_U, epsilon, Li, current_q,M,invM)
{
  q = matrix(current_q,length(current_q),1)
  p = matrix(mvrnorm(1,rep(0,d),M),length(current_q),1) # independent standard normal variates
  current_p = p
  
  # Make a half step for momentum at the beginning
  p = p - epsilon * grad_U(q) / 2
  
  # Alternate full steps for position and momentum
  for (i in 1:Li)
  {
    # Make a full step for the position
    q = q + epsilon * invM%*%p
    # Make a full step for the momentum, except at end of trajectory
    if (i!=Li) p = p - epsilon * grad_U(q)
  }
  
  # Make a half step for momentum at the end.
  p = p - epsilon * grad_U(q) / 2
  
  # Negate momentum at end of trajectory to make the proposal symmetric
  p = -p
  
  # Evaluate potential and kinetic energies at start and end of trajectory
  current_U = U(current_q)
  current_K = t(current_p) %*%invM%*% current_p /2
  proposed_U = U(q)
  proposed_K = t(p) %*% invM%*%p /2
  
  # Accept or reject the state at end of trajectory, returning either
  # the position at the end of the trajectory or the initial position
  if (runif(1) < exp(current_U-proposed_U+current_K-proposed_K))
  {
    return (q) # accept
  }
  else
  {
    return (current_q) # reject
  }
  
}

acc_rate <- function(x){
  num <- dim(x)[2]-1
  count <- 0
  for (i in 1:num){
    if (x[1,i+1]==x[1,i]){
      count<-count+1
    }
  }
  return((num-count)/num)
}

f = function(x) t(x) %*% invtrueM %*% x/2-log(det(M))/2
grad_ <- function(x) invtrueM %*% x
grad_f = function(x,eps = 1e-4){
  n <- length(x)
  theta_map <- cbind(rep(0, n), diag(eps, n)) + matrix(rep(x, n + 1), ncol = n + 1)
  LP <- apply(theta_map,2,f)
  (LP[2:11]-LP[1])/eps
} 

N <- 1000
Li <- 1
d  <- 30
M <-matrix(runif(d^2)*2-1,nrow=d)
M <- t(M) %*% M
trueM <- M
invtrueM <- solve(trueM)
theta <- mvrnorm(1,rep(0,d),trueM)

#################################True M#################################

epsilon <- .021
output.true <- matrix(0,d,N)
output.true[,1] <- theta

for (i in 1:(N-1)){
  if (i%%100 == 0){
    print(i)
  }
  output.true[,i+1] <- HMC (f, grad_, epsilon, Li, output.true[,i],trueM,invtrueM)
}

acc_rate(output.true)
sum(effectiveSize(t(output.true)))

#################################Diagonal M##########################################


M.diag <- diag(diag(trueM),d)
M.diag <- solve(M.diag)
invM.diag <- solve(M.diag)

epsilon <- .065
output.diag <- matrix(0,d,N)
output.diag[,1] <- theta

for (i in 1:(N-1)){
  if (i%%100 == 0){
    print(i)
  }
  output.diag[,i+1] <- HMC (f, grad_, epsilon, Li, output.diag[,i],M.diag,invM.diag)
}

acc_rate(output.diag)
sum(effectiveSize(t(output.diag)))

#################################Diagonal M 2##########################################


M.diag2 <- solve(trueM)
M.diag2 <- diag(diag(M.diag2),d)
invM.diag2 <- solve(M.diag2)

epsilon <- .45
output.diag <- matrix(0,d,N)
output.diag[,1] <- theta

for (i in 1:(N-1)){
  if (i%%100 == 0){
    print(i)
  }
  output.diag[,i+1] <- HMC (f, grad_, epsilon, Li, output.diag[,i],M.diag2,invM.diag2)
}

acc_rate(output.diag)
sum(effectiveSize(t(output.diag)))


##########################################################
L <- chol(trueM)
L <- t(L)

f.c = function(x) t(x) %*% x/2-log(det(M))/2
grad_.c <- function(x)  x


epsilon <- .0045
output.chol <- matrix(0,d,N)
output.chol[,1] <- theta

for (i in 1:(N-1)){
  if (i%%100 == 0){
    print(i)
  }
  output.chol[,i+1] <- HMC (f.c, grad_.c, epsilon, Li, output.chol[,i],diag(1,d,d),diag(1,d,d))
  if (output.chol[1,i] != output.chol[1,i+1]){
    print(sprintf('%s,%s',output.chol[1,i],output.chol[1,i+1]) )
    output.chol[,i+1] <- L %*% output.chol[,i+1]
  }
}

acc_rate(output.chol)
sum(effectiveSize(t(output.chol)))


