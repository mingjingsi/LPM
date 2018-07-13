bLPM <- function(data, X = NULL, alpha = NULL, pi1_ = NULL, maxiter = 1e4, tol = 1e-8, coreNum = 1, verbose = FALSE){

  K    <- length(data) # No. of GWASs
  name <- names(data)  # name of GWASs
  
  # quality control for p-values
  No_small_p <- 0
  for(i in 1:K){
    # check whether p-values are in [0, 1]
    if (any(data[[i]]$p < 0 | data[[i]]$p > 1)) {
      stop("Some p-values are smaller than zero or larger than one. p-value should be ranged between zero and one." )
    }
    
    # set zero p-values to small values to avoid log(0)
    if (any(data[[i]]$p < 1e-30)) {
      No_small_p <- No_small_p + sum(data[[i]]$p < 1e-30)
      
      data[[i]]$p[data[[i]]$p < 1e-30] <- 1e-30
    }
  }
  if (No_small_p != 0){
    message("Info: Some SNPs have p-values close to zero." )
    message("Info: Number of SNPs with p-values close to zero: ", No_small_p)
    message("Info: p-values for these SNPs are set to ", 1e-30 )
  }
  
  # initial values of alpha and pi1_
  if (is.null(alpha)) {
    alpha <- rep(0.1, K)
  }
  else{
    check <- (length(alpha)==K)
    if (check == 0)
      stop("alpha should be a vector of the same length of data.")
  }
  if (is.null(pi1_)) {
    pi1_ <- rep(0.1, K)
  }
  else{
    check <- (length(pi1_)==K)
    if (check == 0)
      stop("pi1_ should be a vector of the same length of data.")
  }
  
  # names and initial values for pair-wise GWAS
  Npairs     <- K*(K-1)/2            # No. of pairs
  name_pair  <- numeric(Npairs)      # name of pairs
  alpha_pair <- matrix(0, 2, Npairs) # initial values of alpha for the pairs
  pi1_pair   <- matrix(0, 2, Npairs) # intial values of pi1_ for the pairs
  pair_id    <- matrix(0, 2, Npairs) # 
  temp_pair <- 0
  for (i in 1:(K-1)){
    for (j in (i+1):K){
      temp_pair <- temp_pair + 1
      name_pair[temp_pair]    <- paste(name[i], name[j], sep = "_")
      alpha_pair[, temp_pair] <- c(alpha[i], alpha[j])
      pi1_pair[, temp_pair]   <- c(pi1_[i], pi1_[j])
      pair_id[, temp_pair]    <- c(i, j)
    }
  }
  
  # report information
  message("Info: Number of GWASs: ", K)
  
  # extract p-values
  Pvalue <- NULL
  for (i in 1:K){
    Pvalue <- c(Pvalue, list(data[[i]]$p))
  }

  # check whether covariates are inputted
  if (is.null(X)) {
    D <- 0                   # No. of covariates
    Pvalue_pair_id  <- NULL  # SNP id for pair data
    for (i in 1:(K-1)){
      for (j in (i+1):K){
        temp1 <- data[[i]]
        temp1$p <- 1:nrow(temp1)
        temp2 <- data[[j]]
        temp2$p <- 1:nrow(temp2)
        temp <- merge(temp1, temp2, by = "SNP")
        Pvalue_pair_id <- c(Pvalue_pair_id, list(as.matrix(temp[, 2:3])))
      }
    }
  }
  else {
    D <- ncol(X) - 1                                # No. of covariates
    X <- data.frame(intercept = rep(1, nrow(X)), X) # add a colume of 1 as the intercept
    Pvalue_X_pair_id  <- NULL                       # SNP id for pair data and X
    for (i in 1:(K-1)){
      for (j in (i+1):K){
        temp1 <- data[[i]]
        temp1$p <- 1:nrow(temp1)
        temp2 <- data[[j]]
        temp2$p <- 1:nrow(temp2)
        temp3 <- data.frame(SNP = X$SNP, id = 1:nrow(X))
        temp <- merge(temp1, temp2, by = "SNP")
        temp <- merge(temp, temp3, by = "SNP")
        Pvalue_X_pair_id <- c(Pvalue_X_pair_id, list(as.matrix(temp[, 2:4])))
      }
    }
    X$SNP <- NULL
    X_name <- names(X) # name of covariates
    X <- as.matrix(X)
  }
  
  rm(data)

  # report information
  message("Info: Number of GWAS pairs: ", Npairs)
  message("Info: Number of covariates: ", D)

  # fit bLPM
  if (D == 0) {
    # If no covarites
    fit <- bLPM_noX_Rcpp(Pvalue, Pvalue_pair_id, alpha_pair, pi1_pair, pair_id, maxiter, tol, coreNum)
    
    colnames(fit$beta0) <- name_pair
    colnames(fit$beta0_stage2) <- name_pair
    
    # compute R
    if (pair == FALSE){
      fit$R <- matrix(0, K, K)
      temp_pair <- 0
      for (i in 1:(K-1)){
        for (j in (i+1):K){
          temp_pair <- temp_pair + 1
          fit$R[i, j] <- fit$rho[temp_pair]
        }
      }
      fit$R <- fit$R + t(fit$R)
      diag(fit$R) <- 1
      rownames(fit$R) <- name
      colnames(fit$R) <- name
    }
    
    if(verbose == FALSE){
      fit$beta0_stage2 <- NULL
    }
  }
  else{
    fit <- bLPM_Rcpp(Pvalue, X, Pvalue_X_pair_id, alpha_pair, pi1_pair, pair_id, maxiter, tol, coreNum)
    
    dimnames(fit$beta)[[2]] <- X_name
    dimnames(fit$beta_stage2)[[2]] <- X_name
    if (Npairs != 1){
      dimnames(fit$beta)[[3]] <- name_pair
      dimnames(fit$beta_stage2)[[3]] <- name_pair
    }
    
    # compute R
    fit$R <- matrix(0, K, K)
    temp_pair <- 0
    for (i in 1:(K-1)){
      for (j in (i+1):K){
        temp_pair <- temp_pair + 1
        fit$R[i, j] <- fit$rho[temp_pair]
      }
    }
    fit$R <- fit$R + t(fit$R)
    diag(fit$R) <- 1
    rownames(fit$R) <- name
    colnames(fit$R) <- name
    
    if(verbose == FALSE){
      fit$beta_stage2 <- NULL
    }
  }

  colnames(fit$alpha) <- name_pair
  colnames(fit$alpha_stage1) <- name_pair
  colnames(fit$alpha_stage2) <- name_pair
  colnames(fit$pi1_stage1) <- name_pair
  rownames(fit$rho) <- name_pair
  names(fit$L_stage1_List) <- name_pair
  names(fit$L_stage2_List) <- name_pair
  names(fit$L_stage3_List) <- name_pair
  
  if(verbose == FALSE){
    fit$alpha_stage1 <- NULL
    fit$alpha_stage2 <- NULL
    fit$pi1_stage1 <- NULL
  }

  return(fit)

}

bLPM_add <- function(data, data_add, X = NULL, fit, maxiter = 1e4, tol = 1e-8, coreNum = 1, verbose = FALSE){
  
  K     <- length(data)     # No. of original GWASs
  K_add <- length(data_add) # No. of added GWASs
  K_all <- K + K_add        # total No. of GWASs
  
  name     <- names(data)       # name of original GWASs
  name_add <- names(data_add)   # name of added GWASs
  name_all <- c(name, name_add) # name of all GWASs
  
  Npairs_all <- (K_all)*(K_all-1)/2     # No. of all pairs
  name_pair_all  <- numeric(Npairs_all) # name of all pairs
  temp_pair <- 0
  for (i in 1:(K_all-1)){
    for (j in (i+1):K_all){
      temp_pair <- temp_pair + 1
      name_pair_all[temp_pair] <- paste(name_all[i], name_all[j], sep = "_")
    }
  }
  
  # quality control for p-values
  No_small_p <- 0
  for(i in 1:K){
    # check whether p-values are in [0, 1]
    if (any(data[[i]]$p < 0 | data[[i]]$p > 1)) {
      stop("Some p-values are smaller than zero or larger than one. p-value should be ranged between zero and one." )
    }
    
    # set zero p-values to small values to avoid log(0)
    if (any(data[[i]]$p < 1e-30)) {
      No_small_p <- No_small_p + sum(data[[i]]$p < 1e-30)
      
      data[[i]]$p[data[[i]]$p < 1e-30] <- 1e-30
    }
  }
  for(i in 1:K_add){
    # check whether p-values are in [0, 1]
    if (any(data_add[[i]]$p < 0 | data_add[[i]]$p > 1)) {
      stop("Some p-values are smaller than zero or larger than one. p-value should be ranged between zero and one." )
    }
    
    # set zero p-values to small values to avoid log(0)
    if (any(data_add[[i]]$p < 1e-30)) {
      No_small_p <- No_small_p + sum(data_add[[i]]$p < 1e-30)
      
      data_add[[i]]$p[data_add[[i]]$p < 1e-30] <- 1e-30
    }
  }
  if (No_small_p != 0){
    message("Info: Some SNPs have p-values close to zero." )
    message("Info: Number of SNPs with p-values close to zero: ", No_small_p)
    message("Info: p-values for these SNPs are set to ", 1e-30 )
  }
  
  # convert to pair-wise GWAS
  Npairs     <- K_all*(K_all-1)/2 - K*(K-1)/2 # No. of added pairs
  alpha_pair <- matrix(0.1, 2, Npairs)        # initial values of alpha for the pairs
  pi1_pair   <- matrix(0.1, 2, Npairs)        # intial values of pi1_ for the pairs
  pair_id    <- matrix(0.1, 2, Npairs)        # added pair id
  temp_pair <- 0
  for (i in 1:K){
    for (j in 1:K_add){
      temp_pair <- temp_pair + 1
      pair_id[, temp_pair] <- c(i, j+K)
    }
  }
  if (K_add > 1){
    for (i in 1:(K_add-1)){
      for (j in (i+1):K_add){
        temp_pair <- temp_pair + 1
        pair_id[, temp_pair] <- c(i+K, j+K)
      }
    }
  }
  
  # report information
  message("Info: Number of original GWASs: ", K)
  message("Info: Number of added GWASs: ", K_add)
  
  # extract p-values
  Pvalue <- NULL
  for (i in 1:K){
    Pvalue <- c(Pvalue, list(data[[i]]$p))
  }
  for (i in 1:K_add){
    Pvalue <- c(Pvalue, list(data_add[[i]]$p))
  }
  
  # check whether covariates are inputted
  if (is.null(X)) {
    D <- 0                 # No. of covariates
    Pvalue_pair_id <- NULL # SNP id for the added pairs
    for (i in 1:K){
      for (j in 1:K_add){
        temp1 <- data[[i]]
        temp1$p <- 1:nrow(temp1)
        temp2 <- data_add[[j]]
        temp2$p <- 1:nrow(temp2)
        temp <- merge(temp1, temp2, by = "SNP")
        Pvalue_pair_id <- c(Pvalue_pair_id, list(as.matrix(temp[, 2:3])))
      }
    }
    if (K_add > 1){
      for (i in 1:(K_add-1)){
        for (j in (i+1):K_add){
          temp1 <- data_add[[i]]
          temp1$p <- 1:nrow(temp1)
          temp2 <- data_add[[j]]
          temp2$p <- 1:nrow(temp2)
          temp <- merge(temp1, temp2, by = "SNP")
          Pvalue_pair_id <- c(Pvalue_pair_id, list(as.matrix(temp[, 2:3])))
        }
      }
    }
  }
  else {
    D <- ncol(X) - 1                                # No. of covariates
    X <- data.frame(intercept = rep(1, nrow(X)), X) # add a colume of 1 as the intercept
    Pvalue_X_pair_id <- NULL # SNP id for the added pairs and X
    for (i in 1:K){
      for (j in 1:K_add){
        temp1 <- data[[i]]
        temp1$p <- 1:nrow(temp1)
        temp2 <- data_add[[j]]
        temp2$p <- 1:nrow(temp2)
        temp3 <- data.frame(SNP = X$SNP, id = 1:nrow(X))
        temp <- merge(temp1, temp2, by = "SNP")
        temp <- merge(temp, temp3, by = "SNP")
        Pvalue_X_pair_id <- c(Pvalue_X_pair_id, list(as.matrix(temp[, 2:4])))
      }
    }
    if (K_add > 1){
      for (i in 1:(K_add-1)){
        for (j in (i+1):K_add){
          temp1 <- data_add[[i]]
          temp1$p <- 1:nrow(temp1)
          temp2 <- data_add[[j]]
          temp2$p <- 1:nrow(temp2)
          temp3 <- data.frame(SNP = X$SNP, id = 1:nrow(X))
          temp <- merge(temp1, temp2, by = "SNP")
          temp <- merge(temp, temp3, by = "SNP")
          Pvalue_X_pair_id <- c(Pvalue_X_pair_id, list(as.matrix(temp[, 2:4])))
        }
      }
    }
    X$SNP <- NULL
    X_name <- names(X) # name of covariates
    X <- as.matrix(X)
  }
  
  # report information
  message("Info: Number of added GWAS pairs: ", Npairs)
  message("Info: Number of covariates: ", D)
  
  # fit bLPM
  if (D == 0) {
    # If no covarites
    fit_add <- bLPM_noX_Rcpp(Pvalue, Pvalue_pair_id, alpha_pair, pi1_pair, pair_id, maxiter, tol, coreNum)
    
    # combine original fit with added fit
    if(verbose == FALSE){
      fit_all <- list(alpha = matrix(0, 2, Npairs_all),
                      beta0 = matrix(0, 2, Npairs_all),
                      rho = numeric(Npairs_all),
                      L_stage1_List = NULL,
                      L_stage2_List = NULL,
                      L_stage3_List = NULL)
      temp_pair_all <- 0
      temp_pair <- 0
      temp_pair_add <- 0
      for (i in 1:(K-1)){
        for (j in (i+1):K){
          temp_pair_all <- temp_pair_all + 1
          temp_pair <- temp_pair + 1
          fit_all$beta0[, temp_pair_all] <- fit$beta0[, temp_pair]
        }
        for (j in 1:K_add){
          temp_pair_all <- temp_pair_all + 1
          temp_pair_add <- temp_pair_add + 1
          fit_all$beta0[, temp_pair_all] <- fit_add$beta0[, temp_pair_add]
        }
      }
      for (i in K:(K_all-1)){
        for (j in (i+1):K_all){
          temp_pair_all <- temp_pair_all + 1
          temp_pair_add <- temp_pair_add + 1
          fit_all$beta0[, temp_pair_all] <- fit_add$beta0[, temp_pair_add]
        }
      }
      
      colnames(fit_all$beta0) <- name_pair_all
    }
    else{
      fit_all <- list(alpha = matrix(0, 2, Npairs_all),
                      alpha_stage1 = matrix(0, 2, Npairs_all),
                      alpha_stage2 = matrix(0, 2, Npairs_all),
                      beta0 = matrix(0, 2, Npairs_all),
                      pi1_stage1 = matrix(0, 2, Npairs_all),
                      beta0_stage2 = matrix(0, 2, Npairs_all),
                      rho = numeric(Npairs_all),
                      L_stage1_List = NULL,
                      L_stage2_List = NULL,
                      L_stage3_List = NULL)
      temp_pair_all <- 0
      temp_pair <- 0
      temp_pair_add <- 0
      for (i in 1:(K-1)){
        for (j in (i+1):K){
          temp_pair_all <- temp_pair_all + 1
          temp_pair <- temp_pair + 1
          fit_all$beta0[, temp_pair_all] <- fit$beta0[, temp_pair]
          fit_all$beta0_stage2[, temp_pair_all] <- fit$beta0_stage2[, temp_pair]
        }
        for (j in 1:K_add){
          temp_pair_all <- temp_pair_all + 1
          temp_pair_add <- temp_pair_add + 1
          fit_all$beta0[, temp_pair_all] <- fit_add$beta0[, temp_pair_add]
          fit_all$beta0_stage2[, temp_pair_all] <- fit_add$beta0_stage2[, temp_pair_add]
        }
      }
      for (i in K:(K_all-1)){
        for (j in (i+1):K_all){
          temp_pair_all <- temp_pair_all + 1
          temp_pair_add <- temp_pair_add + 1
          fit_all$beta0[, temp_pair_all] <- fit_add$beta0[, temp_pair_add]
          fit_all$beta0_stage2[, temp_pair_all] <- fit_add$beta0_stage2[, temp_pair_add]
        }
      }
      
      colnames(fit_all$beta0) <- name_pair_all
      colnames(fit_all$beta0_stage2) <- name_pair_all
    }

  }
  else{
    fit_add <- bLPM_Rcpp(Pvalue, X, Pvalue_X_pair_id, alpha_pair, pi1_pair, pair_id, maxiter, tol, coreNum)
    
    # combine original fit with added fit
    if(verbose == FALSE){
      fit_all <- list(alpha = matrix(0, 2, Npairs_all),
                      beta = array(0, c(2, D+1, Npairs_all)),
                      rho = numeric(Npairs_all),
                      L_stage1_List = NULL,
                      L_stage2_List = NULL,
                      L_stage3_List = NULL)
    }
    else{
      fit_all <- list(alpha = matrix(0, 2, Npairs_all),
                      alpha_stage1 = matrix(0, 2, Npairs_all),
                      alpha_stage2 = matrix(0, 2, Npairs_all),
                      beta = array(0, c(2, D+1, Npairs_all)),
                      pi1_stage1 = matrix(0, 2, Npairs_all),
                      beta_stage2 = array(0, c(2, D+1, Npairs_all)),
                      rho = numeric(Npairs_all),
                      L_stage1_List = NULL,
                      L_stage2_List = NULL,
                      L_stage3_List = NULL)
    }
    
    temp_pair_all <- 0
    temp_pair <- 0
    temp_pair_add <- 0
    for (i in 1:(K-1)){
      for (j in (i+1):K){
        temp_pair_all <- temp_pair_all + 1
        temp_pair <- temp_pair + 1
        fit_all$beta[, , temp_pair_all] <- fit$beta[, , temp_pair]
        if(verbose == TRUE)
          fit_all$beta_stage2[, , temp_pair_all] <- fit$beta_stage2[, , temp_pair]
      }
      for (j in 1:K_add){
        temp_pair_all <- temp_pair_all + 1
        temp_pair_add <- temp_pair_add + 1
        fit_all$beta[, , temp_pair_all] <- fit_add$beta[, , temp_pair_add]
        if(verbose == TRUE)
          fit_all$beta_stage2[, , temp_pair_all] <- fit_add$beta_stage2[, , temp_pair_add]
      }
    }
    for (i in K:(K_all-1)){
      for (j in (i+1):K_all){
        temp_pair_all <- temp_pair_all + 1
        temp_pair_add <- temp_pair_add + 1
        fit_all$beta[, , temp_pair_all] <- fit_add$beta[, , temp_pair_add]
        if(verbose == TRUE)
          fit_all$beta_stage2[, , temp_pair_all] <- fit_add$beta_stage2[, , temp_pair_add]
      }
    }
    
    dimnames(fit_all$beta)[[2]] <- X_name
    dimnames(fit_all$beta)[[3]] <- name_pair_all
    
    if(verbose == TRUE){
      dimnames(fit_all$beta_stage2)[[2]] <- X_name
      dimnames(fit_all$beta_stage2)[[3]] <- name_pair_all
    }
 
  }

  # combine original fit with added fit
  temp_pair_all <- 0
  temp_pair <- 0
  temp_pair_add <- 0
  for (i in 1:(K-1)){
    for (j in (i+1):K){
      temp_pair_all <- temp_pair_all + 1
      temp_pair <- temp_pair + 1
      fit_all$alpha[, temp_pair_all] <- fit$alpha[, temp_pair]
      fit_all$rho[temp_pair_all] <- fit$rho[temp_pair]
      fit_all$L_stage1_List <- c(fit_all$L_stage1_List, fit$L_stage1_List[temp_pair])
      fit_all$L_stage2_List <- c(fit_all$L_stage2_List, fit$L_stage2_List[temp_pair])
      fit_all$L_stage3_List <- c(fit_all$L_stage3_List, fit$L_stage3_List[temp_pair])
      
      if(verbose == TRUE){
        fit_all$alpha_stage1[, temp_pair_all] <- fit$alpha_stage1[, temp_pair]
        fit_all$alpha_stage2[, temp_pair_all] <- fit$alpha_stage2[, temp_pair]
        fit_all$pi1_stage1[, temp_pair_all] <- fit$pi1_stage1[, temp_pair]
      }
    }
    for (j in 1:K_add){
      temp_pair_all <- temp_pair_all + 1
      temp_pair_add <- temp_pair_add + 1
      fit_all$alpha[, temp_pair_all] <- fit_add$alpha[, temp_pair_add]
      fit_all$rho[temp_pair_all] <- fit_add$rho[temp_pair_add]
      fit_all$L_stage1_List <- c(fit_all$L_stage1_List, fit_add$L_stage1_List[temp_pair_add])
      fit_all$L_stage2_List <- c(fit_all$L_stage2_List, fit_add$L_stage2_List[temp_pair_add])
      fit_all$L_stage3_List <- c(fit_all$L_stage3_List, fit_add$L_stage3_List[temp_pair_add])
      
      if(verbose == TRUE){
        fit_all$alpha_stage1[, temp_pair_all] <- fit_add$alpha_stage1[, temp_pair_add]
        fit_all$alpha_stage2[, temp_pair_all] <- fit_add$alpha_stage2[, temp_pair_add]
        fit_all$pi1_stage1[, temp_pair_all] <- fit_add$pi1_stage1[, temp_pair_add]
      }
    }
  }
  for (i in K:(K_all-1)){
    for (j in (i+1):K_all){
      temp_pair_all <- temp_pair_all + 1
      temp_pair_add <- temp_pair_add + 1
      fit_all$alpha[, temp_pair_all] <- fit_add$alpha[, temp_pair_add]
      fit_all$rho[temp_pair_all] <- fit_add$rho[temp_pair_add]
      fit_all$L_stage1_List <- c(fit_all$L_stage1_List, fit_add$L_stage1_List[temp_pair_add])
      fit_all$L_stage2_List <- c(fit_all$L_stage2_List, fit_add$L_stage2_List[temp_pair_add])
      fit_all$L_stage3_List <- c(fit_all$L_stage3_List, fit_add$L_stage3_List[temp_pair_add])
      
      if(verbose == TRUE){
        fit_all$alpha_stage1[, temp_pair_all] <- fit_add$alpha_stage1[, temp_pair_add]
        fit_all$alpha_stage2[, temp_pair_all] <- fit_add$alpha_stage2[, temp_pair_add]
        fit_all$pi1_stage1[, temp_pair_all] <- fit_add$pi1_stage1[, temp_pair_add]
      }
    }
  }

  # compute R
  fit_all$R <- matrix(0, K_all, K_all)
  temp_pair_all <- 0
  for (i in 1:(K_all-1)){
    for (j in (i+1):K_all){
      temp_pair_all <- temp_pair_all + 1
      fit_all$R[i, j] <- fit_all$rho[temp_pair_all]
    }
  }
  fit_all$R <- fit_all$R + t(fit_all$R)
  diag(fit_all$R) <- 1
  rownames(fit_all$R) <- name_all
  colnames(fit_all$R) <- name_all
  
  colnames(fit_all$alpha) <- name_pair_all
  names(fit_all$rho) <- name_pair_all
  names(fit_all$L_stage1_List) <- name_pair_all
  names(fit_all$L_stage2_List) <- name_pair_all
  names(fit_all$L_stage3_List) <- name_pair_all
  
  if(verbose == TRUE){
    colnames(fit_all$alpha_stage1) <- name_pair_all
    colnames(fit_all$alpha_stage2) <- name_pair_all
    colnames(fit_all$pi1_stage1) <- name_pair_all
  }

  return(fit_all)
}

LPM <- function(bLPMfit){
  K <- ncol(bLPMfit$R)
  
  LPMfit <- NULL
  current_pair <- 0
  alpha <- numeric(K)
  
  if (is.null(bLPMfit$beta)){
    beta <- numeric(K)
    for (i in 1:(K-1)){
      for (j in (i+1):K){
        current_pair <- current_pair + 1
        alpha[i] <- alpha[i] + bLPMfit$alpha[1, current_pair]
        alpha[j] <- alpha[j] + bLPMfit$alpha[2, current_pair]
        beta[i] <- beta[i] + bLPMfit$beta0[1, current_pair]
        beta[j] <- beta[j] + bLPMfit$beta0[2, current_pair]
      }
    }
    
    names(beta) <- colnames(bLPMfit$R)
    
  }
  else{
    beta <- matrix(0, K, length(bLPMfit$beta[1, , 1]))
    for (i in 1:(K-1)){
      for (j in (i+1):K){
        current_pair <- current_pair + 1
        alpha[i] <- alpha[i] + bLPMfit$alpha[1, current_pair]
        alpha[j] <- alpha[j] + bLPMfit$alpha[2, current_pair]
        beta[i, ] <- beta[i, ] + bLPMfit$beta[1, , current_pair]
        beta[j, ] <- beta[j, ] + bLPMfit$beta[2, , current_pair]
      }
    }
    
    colnames(beta) <- colnames(bLPMfit$beta[, , 1])
    rownames(beta)<- names(alpha)
    
  }
  names(alpha) <- colnames(bLPMfit$R)
  
  LPMfit$alpha <- alpha/(K-1)
  LPMfit$beta <- beta/(K-1)
  LPMfit$R <- CorrelationMatrix(bLPMfit$R, 1)$CorrMat
  
  colnames(LPMfit$R) <- colnames(bLPMfit$R)
  rownames(LPMfit$R) <- rownames(bLPMfit$R)
  
  return(LPMfit)
}
