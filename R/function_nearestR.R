#################   This R code is designed by Ying Cui at National University of Singapore to solve  ##
#################                          min 0.5*<X-Sigma, X-Sigma>
#################                          s.t. X_ii =b_i, i=1,2,...,n
#################                               X>=tau*I (symmetric and positive semi-definite)  ########
#################
#################                           based on the algorithm  in                   #################
#################                 ``A Quadratically Convergent Newton Method fo           r#################
#################                    Computing the Nearest Correlation Matrix             #################
#################                           By Houduo Qi and Defeng Sun                   #################
#################                     SIAM J. Matrix Anal. Appl. 28 (2006) 360--385.
#################
#################                       Last modified date: March 19, 2016                #################
#################                                                                       #################
#################  The  input arguments  Sigma, b>0, tau>=0, and tol (tolerance error)     #################
#################                                                                      #################
#################                                                                        #################
#################      For the correlation matrix problem, set b = rep(1,n)             #################
#################                                                                       #################
#################      For a positive definite matrix                                     #################
#################        set tau = 1.0e-5 for example                                    #################
#################        set tol = 1.0e-6 or lower if no very high accuracy required     #################
#################        If the algorithm terminates with the given tol, then it means    ##################
#################         the relative gap between the  objective function value of the obtained solution  ###
#################         and the objective function value of the global solution is no more than tol or  ###
#################         the violation of the feasibility is smaller than tol*(1+norm(b,'2'))    ############
#################
#################        The outputs include the optimal primal solution, the dual solution and others   ########
#################                 Diagonal Preconditioner is used                                  #################
#################
#################      Send your comments to hdqi@soton.ac.uk  or matsundf@nus.edu.sg              #################
#################
################# Warning:  Though the code works extremely well, it is your call to use it or not. #################


CorrelationMatrix <- function(Sigma,b,tau=0,tol=1e-6)
{

  Sigma<- as.matrix(Sigma)
  n <- dim(Sigma)[1]

  if (missing(b) )
  {
    b = rep(1,n)
  }




  ptm <- proc.time()    #### time:start

  b <- as.array(b)

  Sigma <- (Sigma+t(Sigma))/2   #### make Sigma symmetric
  b0 <- array(1,c(n,1))
  error_tol <- max(1e-12,tol)

  printyes <- 1
  printsubyes <- 1
  Res_b <- array(0,c(300,1))
  y <- array(0,c(n,1))                #### initial dual value
  Fy <- array(0,c(n,1))



  Iter_whole <- 200
  Iter_inner <- 20
  Iter_CG <- 200

  f_eval <- 0
  iter_k <- 0
  iter_linesearch <- 0
  tol_CG <- 1e-2               #### relative accuracy of the conjugate gradient method for solving the Newton direction
  sigma_1 <- 1e-4              #### tolerance in the line search of the Newton method

  x0 <- y
  c <- array(1,c(n,1))
  d <- array(0,c(n,1))
  val_Sigma <- 0.5*norm(Sigma,'F')^2

  if (n == 1)
  {
    X <- Sigma + y
  }
  else
  {
    X <- Sigma + diag(c(y))
  }
  X <- (X + t(X))/2
  Eigen_X <- Myeigen(X,n)
  Eigen_P <- Eigen_X$vectors
  Eigen_lambda <- Eigen_X$values

  gradient <- Corrsub_gradient(y,Eigen_X,b0,n)
  f0 <- gradient$f
  Fy <- gradient$Fy
  val_dual <- val_Sigma - f0
  X <- PCA(X,Eigen_X,b0,n)
  val_obj <- 0.5*norm(X-Sigma,'F')^2
  gap <- (val_obj - val_dual)/(1+abs(val_dual) + abs(val_obj))

  f <- f0
  f_eval <- f_eval + 1
  b <- b0 - Fy
  norm_b <- norm(b,'2')
  norm_b0 <- norm(b0,'2')+1
  norm_b_rel <- norm_b/norm_b0

  Omega12 <- omega_mat(Eigen_lambda,n)

  # if (c == 1)
  # {
  #   cat(' iter      grad          rel_grad        p_obj           d_obj           rel_gap     CG_iter\n')
  #   cat(iter_k, format(norm_b,scientific = T,digits = 4),format(norm_b_rel,scientific = T,digits = 4),format(val_obj,scientific= T,digits = 4),  format(val_dual, scientific= T,digits = 4), format(gap,scientific= T,digits = 4),'\n',sep="\t")
  # }

  x0 <- y
  while (gap > error_tol & norm_b_rel > error_tol & iter_k< Iter_whole)
  {
    c <- precond_matrix(Omega12,Eigen_P,n)
    CG_result<- pre_cg(b,tol_CG,Iter_CG,c,Omega12,Eigen_P,n)

    d <- CG_result$dir
    if (printsubyes == 1)
    {
      iter_CG <- CG_result$iter
      if (CG_result$flag != 0)
      {
        cat('\n Warning: Not a completed Newton-CG step\n ')
      }

    }


    slope <- sum((Fy - b0)*d)

    y <- x0+d

    if (n == 1)
    {
      X <- Sigma + y
    }
    else
    {
      X <- Sigma + diag(c(y))
    }

    X <- (X+t(X))/2
    Eigen_X <- Myeigen(X,n)
    Eigen_P <- Eigen_X$vectors
    Eigen_lambda <- Eigen_X$values

    gradient <- Corrsub_gradient(y,Eigen_X,b0,n)
    f <- gradient$f
    Fy <- gradient$Fy

    k_inner <- 0
    while(k_inner <= Iter_inner & f>f0 + sigma_1*0.5^k_inner*slope +1e-6)
    {
      k_inner <- k_inner+1
      y <- x0 + 0.5^k_inner*d

      X <- Sigma + diag(c(y))

      X <- (X + t(X))/2
      Eigen_X <- Myeigen(X,n)
      Eigen_P <- Eigen_X$vectors
      Eigen_lambda <- Eigen_X$values

      gradient <- Corrsub_gradient(y,Eigen_X,b0,n)
      f <- gradient$f
      Fy <- gradient$Fy
    }

    f_eval <- f_eval + k_inner + 1
    x0 <- y
    f0 <- f
    val_dual <- val_Sigma  - f0
    X <- PCA(X,Eigen_X,b0,n)
    val_obj <- 0.5*norm(X-Sigma,'F')^2
    gap <- (val_obj - val_dual)/(1+abs(val_dual)+abs(val_obj))

    iter_k <- iter_k + 1
    b <- b0 - Fy
    norm_b <- norm(b,'2')
    Res_b[iter_k] <-  norm_b
    norm_b_rel <- norm_b/norm_b0



    Omega12 <- omega_mat(Eigen_lambda,n)
    # if (printsubyes == 1)
    # {
    #   iter_CG <- CG_result$iter
    #   cat(iter_k, format(norm_b,scientific = T,digits = 4), format(norm_b_rel,scientific = T,digits = 4),format(val_obj,scientific= T,digits = 4),  format(val_dual, scientific= T,digits = 4), format(gap,scientific= T,digits = 4),iter_CG,'\n',sep="\t")
    # }

  }



  X <- X + tau
  Final_f <- val_Sigma - f
  rank_X <- sum(Eigen_lambda >= 0)

  info <- list(CorrMat = X, dual_sol = y,rel_grad = norm_b_rel, rank = rank_X, gap = gap,iterations = iter_k)

  CPU_time<- proc.time() - ptm      ##### time:end


  # if (printyes == 1)
  # { #cat('\n ------Final computational details: the rank of the optimal solution : ', rank_X)
  #   cat('\n -------Computational time: ',CPU_time)
  #
  # }
  return(info)

}

Corrsub_gradient <- function(y,Eigen_X,b,n)
{

  f <- 0
  Fy <- array(0,c(n,1))

  P <- t(Eigen_X$vectors)
  lambda <- Eigen_X$values

  i<- 1
  while ( i<= n)
  {
    P[i,] <- (max(lambda[i],0))^0.5*P[i,]
    f <- f + (max(lambda[i],0))^2
    i <- i+1
  }

  i<-1
  while ( i<= n)
  {
    Fy[i] <- crossprod(P[,i])
    i <- i+1
  }

  f <- 0.5*f - sum(b*y)
  list(f = f, Fy = Fy)
}

PCA <- function(X,Eigen_X,b,n)
{

  if (n ==1)
  {
    X <- b
  }else{
    P <- Eigen_X$vectors
    lambda <- Eigen_X$values
    r <- sum(lambda >0)

    if (r==0)
    {
      X <- mat.or.vec(n, n)
    }
    else if (r == n)
    {
      X <- X
    }
    else if (r <= n/2)
    {
      lambda1 <- lambda[1:r]
      lambda1 <- t(lambda1^0.5)
      P1 <- P[,1:r]
      if (r>1)
      {
        lambda1 <- matrix(rep(lambda1,each = n),nrow = n)
        P1 <- lambda1*P1
        X <- tcrossprod(P1)
      }
      else
      {
        X <- lambda1^2*tcrossprod(P1)
      }

    }
    else
    {
      lambda2 <- -lambda[(r+1):n]
      lambda2 <- t(lambda2^0.5)
      lambda2 <- matrix(rep(lambda2,each = n),nrow = n)
      P2 <- P[,(r+1):n]
      P2 <- lambda2*P2
      X <- X + tcrossprod(P2)
    }




    d <- diag(X)
    d <- pmax(d,b)
    X <- X - diag(diag(X)) + diag(d)

    d <- (b/d)^0.5
    d2 <- tcrossprod(d)
    X <- X*d2
  }

  return(X)
}

omega_mat <- function(lambda,n){
  r <- sum(lambda >0)

  if (r>0)
  {
    if (r==n)
    {
      Omega12 <- matrix(1,n,n)
    }
    else
    {
      s <- n-r
      Omega12 <- matrix(1,nrow = r,ncol = s)
      for (i in 1:r)
      {
        for (j in 1:s)
        {
          Omega12[i,j] <- lambda[i]/(lambda[i]-lambda[r+j]);
        }
      }
    }
  }
  else
  {
    Omega12 <- matrix(nrow=0, ncol=0)
  }
  return(Omega12)
}

pre_cg <- function(b,tol,maxit,c,Omega12,P,n){
  r <- b  ### initial value for CG: zero
  n2b <- norm(b,'2')
  tolb <- tol*n2b
  p <- array(0,c(n,1))
  flag <- 1
  iterk <- 0


  relres <- 1000
  z <- r/c  ####z = M\r; here M =diag(c); if M is not the identity matrix
  rz1 <- sum(r*z)
  rz2 <- 1
  d <- z

  for (k in 1:maxit)
  {
    if (k > 1)
    {
      beta <- rz1/rz2
      d <- z + beta*d
    }
    w <- Jacobian_matrix(d,Omega12,P,n)
    denom <- sum(d*w)
    iterk <- k
    relres <- norm(r,'2')/n2b
    if (denom <= 0 )
    {
      sssss <- 0
      p <- d/norm(d)
      break
    }
    else
    {
      alpha <- rz1/denom
      p <- p + alpha*d
      r <- r - alpha*w
    }
    z <- r/c
    if (norm(r,'2') <= tolb)
    {
      iterk <- k
      relres <- norm(r,'2')/n2b
      flag <- 0
      break
    }
    rz2 <- rz1
    rz1 <- sum(r*z)
  }

  list(dir = p, flag = flag, iter = iterk)

}

Jacobian_matrix <- function(x,Omega12,P,n)
{

  Ax <- array(0,c(n,1))

  Dim <-  dim(Omega12)
  r <- Dim[1]
  s <- Dim[2]

  if (r>0)
  {
    H1 <- P[,1:r]
    if (r == 1)
    {
      H1<- matrix(H1,n,1)
    }
    if (r< n/2)
    {
      i <- 1
      while (i <= n)
      {
        H1[i,] <- x[i] *H1[i,]
        i <- i+1
      }

      Omega12 <- Omega12* (t(H1)%*% P[,(r+1):n])
      HH1 <- t(H1)%*%tcrossprod(P[,1:r])+ Omega12%*%t(P[,(r+1):n])

      HH2 <- t(P[,1:r]%*%Omega12)

      H <- rbind(HH1,HH2)
      i <- 1
      while (i<=n){
        Ax[i] <- P[i,]%*%H[,i];
        Ax[i] <- Ax[i] + 1e-10*x[i]; ### add a small perturbation
        i <- i+1;
      }
    }
    else
    {
      if (r == n )
      {Ax <- (1+1e-10)*x}
      else
      {
        H2 <- P[,(r+1):n]
        if (n-r == 1)
        {
          H2<- matrix(H2,n,1)
        }
        i <- 1
        while (i<=n)
        {
          H2[i,] <- x[i]*H2[i,]
          i <- i+1
        }
        Omega12 <- 1-Omega12
        Omega12 <- Omega12 * (t(P[,1:r])%*%H2)

        HH1 <- Omega12 %*% t(P[,(r+1):n])
        HH2 <- t(Omega12) %*% t(P[,1:r]) + (t(P[,(r+1):n])%*%H2)%*%t(P[,(r+1):n])
        H <- rbind(HH1,HH2)

        i <- 1
        while (i<=n)
        {
          Ax[i] <- -P[i,]%*%H[,i]
          Ax[i] <- x[i] + Ax[i] + 1e-10*x[i]
          i <- i+1
        }
      }
    }
  }
  return(Ax)
}

precond_matrix <- function(Omega12,P,n)
{

  Dim <-  dim(Omega12)
  r <- Dim[1]
  s <- Dim[2]
  c <- array(1,c(n,1))
  if (r >0)
  {
    if ( r<n/2 )
    {
      H <- t(P)
      H <- H*H

      H_temp <- H[1:r,]
      if (r == 1)
      {
        H_temp <- matrix(H_temp,1,n)
      }
      H12 <- t(H_temp) %*% Omega12

      d <- array(1,c(r,1))
      for (i in 1:n)
      {
        c[i] <- sum(H[1:r,i])%*% (t(d)%*%H[1:r,i])
        c[i] <- c[i] + 2*(H12[i,]%*%H[(r+1):n,i])
        if (c[i] <1e-8)
        {
          c[i] <- 1e-8
        }

      }
    }
    else
    {
      if (r<n)
      {
        H <- t(P)
        H <- H*H
        Omega12 <- 1 - Omega12
        H12 <- Omega12 %*% H[(r+1):n,]
        d <- array(1,c(s,1))
        dd <- array(1,c(n,1))

        for (i in 1:n)
        {


          c[i] <- sum(H[(r+1):n,i])%*% (t(d)%*%H[(r+1):n,i])
          c[i] <- c[i] + 2*(t(H[1:r,i])%*%H12[,i])
          alpha <- sum(H[,i])
          c[i] <- alpha * (t(H[,i])%*%dd) - c[i]
          if (c[i] < 1e-8)
          {
            c[i] <-  1e-8
          }

        }
      }
    }
  }


  return(c)
}



Myeigen <- function(X,n)
{
  Eigen_X <- eigen(X,symmetric = TRUE)
  P <- Eigen_X$vectors
  lambda <- Eigen_X$values

  if (is.unsorted(-lambda))   #### arrange the eigenvalues in non-increasing order
  {
    lambda_sort <- sort(lambda,decreasing = TRUE,index.return = TRUE)
    lambda <- lambda_sort$x
    idx <- lambda_new$ix
    P <- P[,idx]
  }

  list(values = lambda, vectors = P)

}




