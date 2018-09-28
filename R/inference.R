library(pbivnorm)
library(mvtnorm)

##### calculate posterior #####
# calculate posterior based on one trait
post <- function(data, X = NULL, id, LPMfit){
  
  if (length(data) != length(id))
    stop("data should be a list of the same length of id.")
  
  if(length(id) == 1){
    
    if (is.null(X)){
      Phi <- pnorm(LPMfit$beta[id])
      
      comp.pi1 <- Phi*LPMfit$alpha[id]
      comp.pi0 <- (1 - Phi)/(data[[1]]$p^(LPMfit$alpha[id] - 1))
      comp.L <- comp.pi1 + comp.pi0
      post <- comp.pi1/comp.L
      
      posterior <- data.frame(SNP = data[[1]]$SNP, posterior = post)
      
    }
    else{
      D <- ncol(X)
      
      posterior <- NULL
      
      data_X <- merge(data[[1]], X, by = "SNP")
      current_X <- cbind(rep(1, nrow(data_X)), data_X[, 3:(D+1)])
      
      Xbeta <- as.vector(LPMfit$beta[id, ]%*%t(current_X))
      Phi <- pnorm(Xbeta)
      
      comp.pi1 <- Phi*LPMfit$alpha[id]
      comp.pi0 <- (1 - Phi)/(data_X[, 2]^(LPMfit$alpha[id] - 1))
      comp.L <- comp.pi1 + comp.pi0
      post <- comp.pi1/comp.L
      
      posterior <- data.frame(SNP = data_X[, 1], posterior = post)
      
    }
  }
  
  if(length(id) == 2){
    data_pair <- merge(data[[1]], data[[2]], by = "SNP")
    
    if (is.null(X)){
      alpha <- LPMfit$alpha[c(id[1], id[2])]
      beta  <- LPMfit$beta[c(id[1], id[2])]
      rho   <- LPMfit$R[id[1], id[2]]
      
      Phi <- pnorm(beta)
      Phi11 <- pbivnorm(beta[1], beta[2], rho)
      Phi10 <- -Phi11 + Phi[1]
      Phi01 <- -Phi11 + Phi[2]
      Phi00 <- 1 + Phi11 - Phi[1] - Phi[2]
      
      comp.pi11 <- Phi11*alpha[1]
      comp.pi10 <- Phi10*alpha[1]/(data_pair[, 3]^(alpha[2] - 1))
      comp.pi01 <- Phi01*alpha[2]/(data_pair[, 2]^(alpha[1] - 1))
      comp.pi00 <- Phi00/(data_pair[, 2]^(alpha[1] - 1)*alpha[2]*data_pair[, 3]^(alpha[2] - 1))
      comp.pi   <- comp.pi11 + comp.pi10 + comp.pi01 + comp.pi00
      
      SNP <- data_pair$SNP
      
    }
    else{
      alpha <- LPMfit$alpha[c(id[1], id[2])]
      beta  <- LPMfit$beta[c(id[1], id[2]), ]
      rho   <- LPMfit$R[id[1], id[2]]
      
      X <- data.frame(intercept = rep(1, nrow(X)), X)
      data_X <- merge(data_pair, X, by = "SNP")
      
      Xbeta <- as.matrix(data_X[,-(1:3)]) %*% t(beta)
      
      Phi <- pnorm(Xbeta)
      Phi11 <- pbivnorm(Xbeta[, 1], Xbeta[, 2], rho)
      Phi10 <- -Phi11 + Phi[, 1]
      Phi01 <- -Phi11 + Phi[, 2]
      Phi00 <- 1 + Phi11 - Phi[, 1] - Phi[, 2]
      
      comp.pi11 <- Phi11*alpha[1]*alpha[2]
      comp.pi10 <- Phi10*alpha[1]/(data_X[, 3]^(alpha[2] - 1))
      comp.pi01 <- Phi01*alpha[2]/(data_X[, 2]^(alpha[1] - 1))
      comp.pi00 <- Phi00/(data_X[, 2]^(alpha[1] - 1)*data_X[, 3]^(alpha[2] - 1))
      comp.pi   <- comp.pi11 + comp.pi10 + comp.pi01 + comp.pi00
      
      SNP <- data_X$SNP
    }
    
    post.joint <- comp.pi11/comp.pi
    post.marginal1 <- (comp.pi11 + comp.pi10)/comp.pi
    post.marginal2 <- (comp.pi11 + comp.pi01)/comp.pi
    
    posterior <- data.frame(SNP = SNP, 
                            post.joint = post.joint,
                            post.marginal1 = post.marginal1, 
                            post.marginal2 = post.marginal2)
  }
  
  if(length(id) == 3){
    data_triple <- merge(data[[1]], data[[2]], by = "SNP")
    data_triple <- merge(data_triple, data[[3]], by = "SNP")
    names(data_triple)[2:4] <- names(data)[id]
    
    if (is.null(X)){
      M <- nrow(data_triple)
      
      current_beta <- LPMfit$beta[id]
      
      Phi111 <- pmvnorm(upper = current_beta, corr = LPMfit$R[id, id])
      Phi110 <- pmvnorm(upper = current_beta*c(1, 1, -1),
                        corr = LPMfit$R[id, id]*matrix(c(1, 1, -1, 1, 1, -1, -1, -1, 1), 3, 3))
      Phi101 <- pmvnorm(upper = current_beta*c(1, -1, 1),
                        corr = LPMfit$R[id, id]*matrix(c(1, -1, 1, -1, 1, -1, 1, -1, 1), 3, 3))
      Phi100 <- pmvnorm(upper = current_beta*c(1, -1, -1),
                        corr = LPMfit$R[id, id]*matrix(c(1, -1, -1, -1, 1, 1, -1, 1, 1), 3, 3))
      Phi011 <- pmvnorm(upper = current_beta*c(-1, 1, 1),
                        corr = LPMfit$R[id, id]*matrix(c(1, -1, -1, -1, 1, 1, -1, 1, 1), 3, 3))
      Phi010 <- pmvnorm(upper = current_beta*c(-1, 1, -1),
                        corr = LPMfit$R[id, id]*matrix(c(1, -1, 1, -1, 1, -1, 1, -1, 1), 3, 3))
      Phi001 <- pmvnorm(upper = current_beta*c(-1, -1, 1),
                        corr = LPMfit$R[id, id]*matrix(c(1, 1, -1, 1, 1, -1, -1, -1, 1), 3, 3))
      Phi000 <- pmvnorm(upper = current_beta*c(-1, -1, -1), corr = LPMfit$R[id, id])
      
      alpha_P1 <- LPMfit$alpha[id[1]]*data_triple[, 2]^(LPMfit$alpha[id[1]] - 1)
      alpha_P2 <- LPMfit$alpha[id[2]]*data_triple[, 3]^(LPMfit$alpha[id[2]] - 1)
      alpha_P3 <- LPMfit$alpha[id[3]]*data_triple[, 4]^(LPMfit$alpha[id[3]] - 1)
      
      SNP <- data_triple$SNP
    }
    else{
      data_X <- merge(data_triple, X, by = "SNP")
      
      M <- nrow(data_X)
      
      current_X <- as.matrix(cbind(rep(1, M), data_X[, 5:(ncol(data_X))]))
      
      Xbeta <- current_X %*% t(LPMfit$beta[id, ])
      
      Phi111 <- numeric(M)
      Phi110 <- numeric(M)
      Phi101 <- numeric(M)
      Phi100 <- numeric(M)
      Phi011 <- numeric(M)
      Phi010 <- numeric(M)
      Phi001 <- numeric(M)
      Phi000 <- numeric(M)
      
      for (j in 1:M){
        Phi111[j] <- pmvnorm(upper = Xbeta[j, ], corr = LPMfit$R[id, id])
        Phi110[j] <- pmvnorm(upper = Xbeta[j, ]*c(1, 1, -1),
                             corr = LPMfit$R[id, id]*matrix(c(1, 1, -1, 1, 1, -1, -1, -1, 1), 3, 3))
        Phi101[j] <- pmvnorm(upper = Xbeta[j, ]*c(1, -1, 1),
                             corr = LPMfit$R[id, id]*matrix(c(1, -1, 1, -1, 1, -1, 1, -1, 1), 3, 3))
        Phi100[j] <- pmvnorm(upper = Xbeta[j, ]*c(1, -1, -1),
                             corr = LPMfit$R[id, id]*matrix(c(1, -1, -1, -1, 1, 1, -1, 1, 1), 3, 3))
        Phi011[j] <- pmvnorm(upper = Xbeta[j, ]*c(-1, 1, 1),
                             corr = LPMfit$R[id, id]*matrix(c(1, -1, -1, -1, 1, 1, -1, 1, 1), 3, 3))
        Phi010[j] <- pmvnorm(upper = Xbeta[j, ]*c(-1, 1, -1),
                             corr = LPMfit$R[id, id]*matrix(c(1, -1, 1, -1, 1, -1, 1, -1, 1), 3, 3))
        Phi001[j] <- pmvnorm(upper = Xbeta[j, ]*c(-1, -1, 1),
                             corr = LPMfit$R[id, id]*matrix(c(1, 1, -1, 1, 1, -1, -1, -1, 1), 3, 3))
        Phi000[j] <- pmvnorm(upper = Xbeta[j, ]*c(-1, -1, -1), corr = LPMfit$R[id, id])
        
      }
      
      alpha_P1 <- LPMfit$alpha[id[1]]*data_X[, 2]^(LPMfit$alpha[id[1]] - 1)
      alpha_P2 <- LPMfit$alpha[id[2]]*data_X[, 3]^(LPMfit$alpha[id[2]] - 1)
      alpha_P3 <- LPMfit$alpha[id[3]]*data_X[, 4]^(LPMfit$alpha[id[3]] - 1)
      
      SNP <- data_X$SNP
    }
    
    comp_pi111 <- Phi111
    comp_pi110 <- Phi110/alpha_P3
    comp_pi101 <- Phi101/alpha_P2
    comp_pi100 <- Phi100/(alpha_P2*alpha_P3)
    comp_pi011 <- Phi011/alpha_P1
    comp_pi010 <- Phi010/(alpha_P1*alpha_P3)
    comp_pi001 <- Phi001/(alpha_P1*alpha_P2)
    comp_pi000 <- Phi000/(alpha_P1*alpha_P2*alpha_P3)
    
    comp_pi <- comp_pi111 + comp_pi110 + comp_pi101 + comp_pi100 + comp_pi011 +
      comp_pi010 + comp_pi001 + comp_pi000
    
    post.joint <- comp_pi111/comp_pi
    post.marginal1 <- (comp_pi111 + comp_pi110 + comp_pi101 + comp_pi100)/comp_pi
    post.marginal2 <- (comp_pi111 + comp_pi110 + comp_pi011 + comp_pi010)/comp_pi
    post.marginal3 <- (comp_pi111 + comp_pi101 + comp_pi011 + comp_pi001)/comp_pi
    post.marginal12 <- (comp_pi111 + comp_pi110)/comp_pi
    post.marginal13 <- (comp_pi111 + comp_pi101)/comp_pi
    post.marginal23 <- (comp_pi111 + comp_pi011)/comp_pi
    
    posterior <- data.frame(SNP = SNP, 
                            post.joint = post.joint,
                            post.marginal1 = post.marginal1, 
                            post.marginal2 = post.marginal2,
                            post.marginal3 = post.marginal3, 
                            post.marginal12 = post.marginal12,
                            post.marginal13 = post.marginal13, 
                            post.marginal23 = post.marginal23)
  }
  
  return(posterior)
}

##### inference for associated SNPs #####
assoc <- function(post, FDRset = 0.1, fdrControl){
  K <- ncol(post)

  eta <- post

  for (k in 2:K){
    
    eta[, k] <- 0
    
    if (fdrControl == "global"){
      eta[which(post2FDR(post[, k]) <= FDRset), k] <- 1
    }
    if (fdrControl == "local"){
      eta[which((1 - post[, k]) <= FDRset), k] <- 1
    }
  }
  
  if (K == 2){
    names(eta)[2] <- "eta"
  }
  if (K == 4){
    names(eta)[2:4] <- c("eta.joint", "eta.marginal1", "eta.marginal2")
  }
  if (K == 8){
    names(eta)[2:8] <- c("eta.joint", "eta.marginal1", "eta.marginal2",  "eta.marginal3",
                         "eta.marginal12",  "eta.marginal13",  "eta.marginal23")
  }

  return(eta)
}

###### transform posterior to FDR #####
post2FDR <- function(posterior){

  M          <- length(posterior)
  fdr        <- 1 - posterior
  rank.fdr   <- rank(fdr)
  sort.fdr   <- sort(fdr)
  cumsum.fdr <- cumsum(sort.fdr)
  sort.FDR   <- cumsum.fdr/seq(1, M, 1)
  FDR        <- sort.FDR[rank.fdr]

  return(FDR)
}

##### test rho=0 #####
test_rho <- function(bLPMfit){

  Npairs <- length(bLPMfit$rho)
  K <- nrow(bLPMfit$R)

  p_value <- matrix(0, K, K)

  current_pair <- 0
  for (i in 1:(K-1)){
    for (j in (i+1):K){
      current_pair <- current_pair + 1
      lambda <- 2*(bLPMfit$L_stage3_List[[current_pair]][length(bLPMfit$L_stage3_List[[current_pair]])] -
              bLPMfit$L_stage2_List[[current_pair]][length(bLPMfit$L_stage2_List[[current_pair]])])
      p_value[i, j] <- 1 - pchisq(lambda, 1)
    }
  }

  p_value <- p_value + t(p_value)
  diag(p_value) <- 0

  rownames(p_value) <- rownames(bLPMfit$R)
  colnames(p_value) <- colnames(bLPMfit$R)

  return(p_value)
}

##### test beta=0 #####
test_beta <- function(data, X, id, LPMfit){
  
  data_X <- merge(data[[id]], X, by = "SNP")
  
  D <- ncol(X)
  M <- nrow(data_X)
  
  current_X <- as.matrix(cbind(rep(1, nrow(data_X)), data_X[, 3:(D+1)]))
  alpha <- LPMfit$alpha[id]
  beta <- LPMfit$beta[id, ]
  Pvalue <- data_X[, 2]
  
  Xbeta <- as.vector(current_X%*%(as.matrix(beta)))
  Phi <- pnorm(Xbeta)
  phi <- dnorm(Xbeta)
  
  comp <- Phi*alpha*Pvalue^(alpha- 1) + 1 - Phi
  comp_g_alpha <- Phi*Pvalue^(alpha - 1)*(1 + alpha*log(Pvalue))
  comp_g_beta <- (alpha*Pvalue^(alpha - 1) - 1)*phi
  comp_H_alpha <- Phi*Pvalue^(alpha-1)*(2+alpha*log(Pvalue))*log(Pvalue)
  comp_H_alpha_beta <- Pvalue^(alpha-1)*(1+alpha*log(Pvalue))*phi
  
  H_alpha <- sum(-(comp_g_alpha/comp)^2 + comp_H_alpha/comp)
  H_alpha_beta <- t(current_X)%*%(-comp_g_alpha/comp*comp_g_beta/comp + comp_H_alpha_beta/comp)
  H_beta <- (t(current_X)*matrix(rep(-(comp_g_beta/comp)^2 - comp_g_beta*Xbeta/comp, each = D), D, M))%*%current_X
  H <- matrix(0, D+1, D+1)
  H[1, 1] <- H_alpha
  H[2:(D+1), 1] <- H_alpha_beta
  H[1, 2:(D+1)] <- t(H_alpha_beta)
  H[2:(D+1), 2:(D+1)] <- H_beta
  
  inv_H <- solve(-H)
  
  W <- beta^2/diag(inv_H)[-1]      
  se <- sqrt(diag(inv_H)[-1])
  p_value <- 1 - pchisq(W, 1)
  
  return(list(p_value = p_value, se = se))
}

# Louise method
test_beta_louise <- function(data, X, id, LPMfit){

  data_X <- merge(data[[id]], X, by = "SNP")

  D <- ncol(X)

  current_X <- as.matrix(cbind(rep(1, nrow(data_X)), data_X[, 3:(D+1)]))
  alpha <- LPMfit$alpha[id]
  beta <- LPMfit$beta[id, ]
  Pvalue <- data_X[, 2]

  Xbeta <- as.vector(current_X%*%(as.matrix(beta)))
  Phi <- pnorm(Xbeta)
  phi <- dnorm(Xbeta)

  comp.pi1 <- Phi*alpha*Pvalue^(alpha- 1)
  comp.pi0 <- 1 - Phi
  comp.L <- comp.pi1 + comp.pi0
  E_eta <- comp.pi1/comp.L

  Ic <- matrix(0, D+1, D+1)
  Ic[1,1] <- sum(E_eta)/alpha^2
  Ic[2:(D+1), 2:(D+1)] <- t(current_X)%*%current_X

  E_SS <- matrix(0, (D+1), (D+1))
  u <- 1/alpha + log(Pvalue)
  v <- phi*(alpha*Pvalue^(alpha-1) - 1)/comp.L
  E_Z <- Xbeta + v
  E_Z2 <- Xbeta*E_Z + 1
  V3 <- E_eta*(Xbeta + phi/Phi - E_Z)
  E_SS[1, 1] <- sum(u^2*(E_eta - E_eta^2))
  E_SS[2:(D+1), 2:(D+1)] <- t(current_X)%*%((E_Z2 - E_Z^2)*current_X)
  E_SS[1, 2:(D+1)] <- (u*V3)%*%current_X
  E_SS[2:(D+1), 1] <- t(E_SS[1, 2:(D+1)])

  I <- Ic - E_SS

  invI <- solve(I)

  W <- LPMfit$beta[id, ]^2/diag(invI)[-1]
  se <- sqrt(diag(invI)[-1])
  p_value <- 1 - pchisq(W, 1)

  return(list(p_value = p_value, se = se))
}
