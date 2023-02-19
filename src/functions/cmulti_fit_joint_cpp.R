
cmulti_fit_joint <- function (Yarray, # Array with dimensions (nsurvey x nrint x ntint)
                              rarray, # distance intervals for each point count
                              tarray, # time intervals for each point count
                              X1 = NULL,     # Design matrix for distance modeling
                              X2 = NULL,     # Design matrix for removal modeling
                              maxdistint = 10, # Max distance for numerical integration
                              tau_inits = NULL, 
                              phi_inits = NULL,
                              method = "Nelder-Mead", ...) {
  
  
  
  input_data <- list(Yarray = Yarray,
                     rarray = rarray,
                     tarray = tarray,
                     X1 = X1,
                     X2 = X2,
                     maxdistint = maxdistint)
  
  # Only conduct analysis on point counts with non-zero total counts
  Ysum <- apply(Yarray,1,sum,na.rm = TRUE)
  Ykeep <- which(Ysum > 0)
  if (length(Ykeep) != length(Ysum)){
    Yarray <- Yarray[Ykeep, , ]
    rarray<- rarray[Ykeep, ]
    tarray<- tarray[Ykeep, ]
  }
  nsurvey <- dim(Yarray)[1] # Number of surveys
  nrint <- apply(rarray,1,function(x)length(na.omit(x))) # Number of distance bins for each point count
  ntint <- apply(tarray,1,function(x)length(na.omit(x))) # Number of time bins for each point count
  
  # Format parameters and check they are named via column names in design matrix
  if (!is.null(X1)){
    X1 <- X1[Ykeep, ]
    tau_params <- colnames(X1)
  } else {
    X1 <- matrix(1,nrow = nsurvey,ncol=1)
    tau_params <- colnames(X1)[1] <- "log_tau"
  }
  
  if (!is.null(X2)){
    X2 <- X2[Ykeep, ]
    phi_params <- colnames(X2)
  } else {
    X2 <- matrix(1,nrow = nsurvey,ncol=1)
    phi_params <- colnames(X2)[1] <- "log_phi"
  }
  
  # Calculate maximum distance for integration at each point count
  max_r <- apply(rarray,1,max,na.rm = TRUE)
  max_r[max_r == Inf] <- maxdistint
  
  # Initial values
  if (length(tau_inits) != ncol(X1)) tau_inits <- NULL
  if (is.null(tau_inits)) {
    tau_inits <- rep(0, ncol(X1))
    names(tau_inits) <- tau_params
  }
  
  if (length(phi_inits) != ncol(X2)) phi_inits <- NULL
  if (is.null(phi_inits)) {
    phi_inits <- rep(0, ncol(X2))
    names(phi_inits) <- phi_params
  }
  
  inits <- c(tau_inits,phi_inits)
  
  
  # browser()
  nlimit <- c(.Machine$double.xmin, .Machine$double.xmax)^(1/3)
  Yarray_x = aperm(Yarray, c(2,3,1))
  res <- optim(inits, nll_fun, method = method, hessian = TRUE,
               X1=X1, X2=X2, tau_params=tau_params, nsurvey=nsurvey, 
               Yarray = Yarray_x,tarray= tarray, rarray=rarray, 
               nrint= nrint, ntint= ntint, max_r= max_r, Ysum=Ysum,nlimit= nlimit)
  
  rval <- list(input_data = input_data,
               coefficients = res$par, 
               #vcov = try(.solvenear(res$hessian)), 
               loglik = -res$value)
  
  #if (inherits(rval$vcov, "try-error")) rval$vcov <- matrix(NA, length(rval$coefficients), length(rval$coefficients))
  rval$coefficients <- unname(rval$coefficients)
  rval$vcov <- unname(rval$vcov)
  rval$results <- res
  rval
}
