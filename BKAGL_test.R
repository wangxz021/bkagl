BKAGL_test <- function(Kmt, state) {
  N <- dim(Kmt[[1]])[1]
  D <- dim(Kmt[[1]])[2]
  T <- length(Kmt)
  P <- c(1:T)
  for(i in 1:T){
    P[i]<-dim(Kmt[[i]])[3]
  }
  
  G <- list()
  for(t in 1:T){
    G[[t]]<-list(mu = matrix(rnorm(P[t] * D ), P[t], D), sigma = diag(1, P[t], P[t])) 
  }
  for(t in 1:T){
    for (m in 1:P[t]){
      G[[t]]$mu[m,] <- crossprod(state$a$mu, Kmt[[t]][,,m])
    }
  }
  
  L <- list(mu = matrix(0, T, D), sigma = matrix(0, T, D))
  for (t in 1:T) {
    L$mu[t,] <- t(state$b[[t]]$mu)%*% G[[t]]$mu
  }
  
  f <- list(mu = matrix(0, D, 1), sigma = matrix(0, D, 1))
  f$mu <- crossprod(rbind(matrix(1, 1, D), L$mu), state$ec$mu)
  f$sigma <- 1 + diag(crossprod(rbind(matrix(1, 1, D), L$mu), state$ec$sigma) %*% rbind(matrix(1, 1, D), L$mu))
  
  pos <- 1 - pnorm((+state$parameters$margin - f$mu) / f$sigma)
  neg <- pnorm((-state$parameters$margin - f$mu) / f$sigma)
  p <- pos / (pos + neg)
  
  prediction <- list(f = f, p = p)
}