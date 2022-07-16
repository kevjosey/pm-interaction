# pm1.5 contrast
pm_contrast <- function(model, pm0, pm1) {
  
  V <- model$var2
  
  X <- data.frame(CVD = c(0,1,0,1), PM = c(pm0, pm0, pm1, pm1))
  X$pm_nocvd <- ifelse(X$CVD == 0, X$PM, 23)
  X$pm_cvd <- ifelse(X$CVD == 1, X$PM, 23)
  
  W.tmp <- model.matrix(delete.response(model$terms), X)[,-1]
  W <- W.tmp - matrix(rep(model$means, each = nrow(W.tmp)), nrow = nrow(W.tmp))
  
  lp <- predict(model, newdata = X, se = FALSE, type = "lp")
  
  h2 <- c(lp[2] - lp[1])
  h3 <- c(lp[3] - lp[1])
  h4 <- c(lp[4] - lp[1])
  
  w2 <- c(W[2,] - W[1,])
  w3 <- c(W[3,] - W[1,])
  w4 <- c(W[4,] - W[1,])
  
  est <- c(h2, h3, h4)
  estvar <- diag(rbind(w2,w3,w4) %*% V %*% cbind(w2,w3,w4))
  lower <- exp(est - 1.96*sqrt(estvar))
  upper <- exp(est + 1.96*sqrt(estvar))
  
  cbind(log.hr = est, log.var = estvar, hr = exp(est),
        lower = lower, upper = upper,
        pm0 = rep(pm0, 3), pm1 = rep(pm1,3))
  
}

# Cox model additive interaction with Splines
add_interact_cox <- function(model, pm0, pm1, conf.level = 0.95) {
  
  lvl <- 1 - ((1 - conf.level)/2)
  z <- qnorm(lvl, mean = 0, sd = 1)
  
  X <- data.frame(CVD = c(0,1,0,1), PM = c(pm0, pm0, pm1, pm1))
  X$pm_nocvd <- ifelse(X$CVD == 0, X$PM, 23)
  X$pm_cvd <- ifelse(X$CVD == 1, X$PM, 23)
  
  lp <- predict(model, newdata = X, se = FALSE, type = "lp")
  V <- model$var2
  
  W.tmp <- model.matrix(delete.response(model$terms), X)[,-1]
  W <- W.tmp #- matrix(rep(model$means, each = nrow(W.tmp)), nrow = nrow(W.tmp))
  
  col1 <- length(model$coefficients)
  col0 <- (col1 - 1)/2 + 1
  
  # these indexes change based on df of pspline()
  a2 <- W[2,2:col0] - W[1,2:col0]
  a3 <- W[3,2:col0] - W[1,2:col0]
  a4 <- W[4,2:col0] - W[1,2:col0]
  b2 <- W[2,(col0 + 1):col1] - W[1,(col0 + 1):col1]
  b3 <- W[3,(col0 + 1):col1] - W[1,(col0 + 1):col1]
  b4 <- W[4,(col0 + 1):col1] - W[1,(col0 + 1):col1]
  
  h0 <- c(exp(lp[4] - lp[1]) - exp(lp[2] - lp[1]))
  h1 <- c(a4*exp(lp[4] - lp[1]) - a3*exp(lp[3] - lp[1]) - a2*exp(lp[2] - lp[1]))
  h2 <- c(b4*exp(lp[4] - lp[1]) - b3*exp(lp[3] - lp[1]) - b2*exp(lp[2] - lp[1]))
  
  reri.var <- t(c(h0,h1,h2)) %*% V %*% c(h0,h1,h2) 
  reri.se <- sqrt(reri.var)
  reri.p <- c(exp(lp[4] - lp[1]) - exp(lp[3] - lp[1]) - exp(lp[2] - lp[1])) + 1
  reri.l <- reri.p - (z * reri.se)
  reri.u <- reri.p + (z * reri.se)
  rval <- c(est = unname(reri.p), lower = reri.l, upper = reri.u, pm0 = pm0, pm1 = pm1)
  
  return(rval)
  
}

# Aalen model additive interaction
add_interact_aalen <- function(model, pm0, pm1, idx = NULL, conf.level = 0.95, monotone = TRUE, ...) {
  
  lvl <- 1 - ((1 - conf.level)/2)
  z <- qnorm(lvl, mean = 0, sd = 1)
  
  if (class(model)[1] != "aalen") 
    stop("Error: model must be aalen object")
  
  design <- rbind(c(1,0,pm0,0), 
                  c(1,1,pm0,pm0),
                  c(1,0,pm1,0), 
                  c(1,1,pm1,pm1))
  
  if (monotone) {
    
    out <- sapply(1:length(idx), function(i, ...) {
      
      beta <- model$cum[idx[i],-1]
      time <- unname(model$cum[idx[i],1])
      lp <- c(design%*%beta)
      V <- model$covariance[[idx[i]]]
      
      reri.p <- c(exp(-lp[2]) + exp(-lp[3]) - exp(-lp[4]) - exp(-lp[1]))
      
      h0 <- -c(reri.p)
      h1 <- c(exp(-lp[4]) - exp(-lp[2]))
      h2 <- c(pm1*exp(-lp[4]) - pm1*exp(-lp[3]) - pm0*exp(-lp[2]) + pm0*exp(-lp[1]))
      h3 <- c(pm1*exp(-lp[4]) - pm0*exp(-lp[2]))
      
      reri.var <- t(c(h0,h1,h2,h3)) %*% V %*% c(h0,h1,h2,h3) 
      reri.se <- sqrt(reri.var)
      reri.l <- reri.p - (z * reri.se)
      reri.u <- reri.p + (z * reri.se)
      rval <- c(time = time, pm0 = pm0, pm1 = pm1, 
                est = reri.p, lower = reri.l, upper = reri.u)
      
      return(rval)
      
    })
    
  } else {
    
    out <- sapply(1:length(idx), function(i, ...) { 
      
      beta <- model$cum[idx[i],-1]
      time <- unname(model$cum[idx[i],1])
      lp <- c(design%*%beta)
      V <- model$covariance[[idx[i]]]
      
      reri.p <- c(exp(-lp[2]) + exp(-lp[3]) - exp(-lp[4]) - 1)
      
      h0 <- -c(reri.p + 1)
      h1 <- c(exp(-lp[4]) - exp(-lp[2]))
      h2 <- c(pm1*exp(-lp[4]) - pm1*exp(-lp[3]) - pm0*exp(-lp[2]) )
      h3 <- c(pm1*exp(-lp[4]) - pm0*exp(-lp[2]))
      
      reri.var <- t(c(h0,h1,h2,h3)) %*% V %*% c(h0,h1,h2,h3) 
      reri.se <- sqrt(reri.var)
      reri.l <- reri.p - (z * reri.se)
      reri.u <- reri.p + (z * reri.se)
      rval <- c(time = time, pm0 = pm0, pm1 = pm1,
                est = reri.p, lower = reri.l, upper = reri.u)
      
      return(rval)
      
    })
  }
  
  return(out)
  
}
