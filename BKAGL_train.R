BKAGL_train <- function(Km, y, parameters) {    
  set.seed(parameters$seed)
  N <- dim(Km[[1]])[1]
  T <- length(Km)
  P <- c(1:T)
  for(i in 1:T){
    P[i]<-dim(Km[[i]])[3]
  }
  
  sigma_g <- parameters$sigma_g
  
  lambda <- list(alpha = matrix(parameters$alpha_lambda + 0.5, N, 1), beta = matrix(parameters$beta_lambda, N, 1))
  a <- list(mu = matrix(rnorm(N), N, 1), sigma = diag(1, N, N))
  G <- list()
  for(i in 1:T){
    G[[i]]<-list(mu = matrix(rnorm(P[i] * N ), P[i], N), sigma = diag(1, P[i], P[i])) 
  }
  
  eta <- list()
  for(i in 1:T){
    eta[[i]]<-list(alpha = matrix(parameters$alpha_eta + 0.5, P[i], 1), beta = matrix(parameters$beta_eta, P[i], 1)) 
  }
  
  b <- list()
  for(i in 1:T){
    b[[i]]<-list(mu = matrix(rnorm(P[i]), P[i], 1), sigma = diag(1, P[i], P[i])) 
  }
  
  L <- list(mu = (abs(matrix(rnorm(T * N), T, N)) + parameters$margin) * sign(matrix(y, T, N, byrow = TRUE)), sigma = diag(1, T, T))
  gamma <- list(alpha = matrix(parameters$alpha_gamma + 0.5, T, 1), beta = matrix(parameters$beta_gamma, T, 1))
  omega <- list(alpha = parameters$alpha_omega + 0.5, beta = parameters$beta_omega)
  ec <- list(mu = rbind(0, matrix(1, T, 1)), sigma = diag(1, T + 1, T + 1))
  f <- list(mu = (abs(matrix(rnorm(N), N, 1)) + parameters$margin) * sign(y), sigma = matrix(1, N, 1))
  
  KmKmT <- matrix(0,N,N)
  for(i in 1:T){
    
    KmKmT <-  KmKmT +  tcrossprod(matrix(Km[[i]],N,N*P[i]))
    
  }
  tep <- matrix(0,N,1) 
  KmtimesGT.mu <- matrix(0,N,1)
  for(i in 1:T){
    KmtimesGT.mu <-  tep  +  matrix(Km[[i]],N,N*P[i])%*% matrix(t(G[[i]]$mu),N*P[i],1)
    tep <- KmtimesGT.mu
  }
  
  lower <- matrix(-1e40, N, 1)
  lower[which(y > 0)] <- +parameters$margin
  upper <- matrix(+1e40, N, 1)
  upper[which(y < 0)] <- -parameters$margin
  
  btimesbT.mu <- list()
  for(t in 1:T){
    btimesbT.mu[[t]] <- b[[t]]$mu %*% t(b[[t]]$mu) + b[[t]]$sigma
  }
  
  atimesaT.mu <- tcrossprod(a$mu, a$mu) + a$sigma
  
  GtimesGT.mu <- list()
  for(t in 1:T){
    GtimesGT.mu[[t]] <- tcrossprod(G[[t]]$mu) + N * G[[t]]$sigma
  }
  btimesGT.mu <- matrix(1, T, N)
  LtimesLT.mu <- tcrossprod(L$mu, L$mu) + N * L$sigma
  ctimescT.mu <- tcrossprod(ec$mu[2:(T + 1)], ec$mu[2:(T + 1)]) + ec$sigma[2:(T + 1), 2:(T + 1)]
  etimeseT.mu <- ec$mu[1]^2 + ec$sigma[1, 1]
  ctimese.mu <- ec$mu[2:(T + 1)] * ec$mu[1] + ec$sigma[2:(T + 1), 1]
  
  for (iter in 1:parameters$iteration) {
    
    # update lambda
    lambda$beta <- 1 / (1 / parameters$beta_lambda + 0.5 * diag(atimesaT.mu))
    
    # update a
    a$sigma <- chol2inv(chol(diag(as.vector(lambda$alpha * lambda$beta), N, N) + KmKmT))
    a$mu <- a$sigma %*% KmtimesGT.mu / sigma_g^2
    atimesaT.mu <- tcrossprod(a$mu, a$mu) + a$sigma
    
    # update G
    for(t in 1:T){
      G[[t]]$sigma <- chol2inv(chol(diag(1, P[t], P[t]) + btimesbT.mu[[t]]))
      for(j in 1:N){
        G[[t]]$mu[,j] <- G[[t]]$sigma %*% (matrix(t(a$mu) %*% matrix(Km[[t]],N,N*P[t]),P[t],N,byrow=T)[,j] + L$mu[t,j] * b[[t]]$mu)
      }
    }
    tep <- matrix(0,N,1) 
    for(i in 1:T){
      KmtimesGT.mu <-  tep  +  matrix(Km[[i]],N,N*P[i])%*% matrix(t(G[[i]]$mu),N*P[i],1)
      tep <- KmtimesGT.mu
    }
    
    # update eta
    for(t in 1:T){
      eta[[t]]$beta <- 1 / (1 / parameters$beta_eta + 0.5 * diag(btimesbT.mu[[t]]))
    }
    
    # update b
    for(t in 1:T){
      b[[t]]$sigma <- chol2inv(chol(diag(as.vector(eta[[t]]$alpha * eta[[t]]$beta), P[t], P[t]) +((tcrossprod(G[[t]]$mu) + N * G[[t]]$sigma)/ sigma_g^2)))
      b[[t]]$mu <- b[[t]]$sigma %*% (matrix(G[[t]]$mu,P[t],N) %*% t(matrix(L$mu[t,],1,N))) / sigma_g^2
      btimesbT.mu[[t]] <- b[[t]]$mu %*% t(b[[t]]$mu) + b[[t]]$sigma
    }
    
    # update L
    L$sigma <- chol2inv(chol(diag(1, T, T) / sigma_g^2 + ctimescT.mu))
    for(t in 1:T){
      for(j in 1:N){
        btimesGT.mu[t,j] <- t(b[[t]]$mu) %*% G[[t]]$mu[,j]
      }
    }
    L$mu <- L$sigma %*% (btimesGT.mu / sigma_g^2 + tcrossprod(ec$mu[2:(T + 1)], f$mu) - matrix(ctimese.mu, T, N, byrow = FALSE))
    LtimesLT.mu <- tcrossprod(L$mu, L$mu) + N * L$sigma
    
    # update gamma
    gamma$beta <- 1 / (1 / parameters$beta_gamma + 0.5 * diag(ctimescT.mu))
    
    # update omega
    omega$beta <- 1 / (1 / parameters$beta_omega + 0.5 * etimeseT.mu)
    
    # update e and c
    ec$sigma <- chol2inv(chol(rbind(cbind(omega$alpha * omega$beta + N, t(rowSums(L$mu))), cbind(rowSums(L$mu), diag(as.vector(gamma$alpha * gamma$beta), T, T) + LtimesLT.mu))))
    ec$mu <- ec$sigma %*% (rbind(matrix(1, 1, N), L$mu) %*% f$mu)
    etimeseT.mu <- ec$mu[1]^2 + ec$sigma[1, 1]
    ctimescT.mu <- tcrossprod(ec$mu[2:(T + 1)], ec$mu[2:(T + 1)]) + ec$sigma[2:(T + 1), 2:(T + 1)]
    ctimese.mu <- ec$mu[2:(T + 1)] * ec$mu[1] + ec$sigma[2:(T + 1), 1]
    
    # update f
    output <- crossprod(rbind(matrix(1, 1, N), L$mu), ec$mu)
    alpha_norm <- lower - output
    beta_norm <- upper - output
    normalization <- pnorm(beta_norm) - pnorm(alpha_norm)
    normalization[which(normalization == 0)] <- 1
    f$mu <- output + (dnorm(alpha_norm) - dnorm(beta_norm)) / normalization
    f$sigma <- 1 + (alpha_norm * dnorm(alpha_norm) - beta_norm * dnorm(beta_norm)) / normalization - (dnorm(alpha_norm) - dnorm(beta_norm))^2 / normalization^2
  }
  
  state <- list(lambda = lambda, a = a, eta = eta, b = b, gamma = gamma, omega = omega, ec = ec, parameters = parameters)
  
}