library(data.table)
library(ranger)
library(xgboost)
library(gee)

# Function for IPW-TM using random forest
ipwtm_rf <- function(exposure, numerator = NULL, denominator, id, timevar, data, trunc = NULL,
                     continuous = FALSE, censor = FALSE, num.trees = 100, max.depth = 8, num.threads = 12, ...) {
  
  # Order the data based on id and timevar
  order.orig <- order(data[[id]], data[[timevar]])
  data <- data.table(data)[order.orig]
  tempdat <- data.table(
    id = data[[id]],
    timevar = data[[timevar]],
    exposure = data[[exposure]]
  )
  
  if (!continuous) {
    
    # Calculate numerator probabilities
    if (is.null(numerator)) {
      mod1.pred <- mean(tempdat$exposure, na.rm = TRUE)
    } else {
      mod1 <- ranger(update.formula(numerator, formula(paste0("factor(", exposure,") ~ ."))), probability = TRUE,
                     data = data, num.trees = num.trees, max.depth = max.depth, num.threads = num.threads,
                     respect.unordered.factors = "partition")
      mod1.pred <- mod1$predictions[, 2]
    }
    
    tempdat$p.numerator <- 1 - mod1.pred
    
    if (!censor)
      tempdat$p.numerator[tempdat$exposure == 1] <- mod1.pred[tempdat$exposure == 1]
    
    # Calculate denominator probabilities
    mod2 <- ranger(update.formula(denominator, formula(paste0("factor(", exposure,") ~ ."))), probability = TRUE,
                   data = data, num.trees = num.trees, max.depth = max.depth, num.threads = num.threads,
                   respect.unordered.factors = "partition")
    mod2.pred <- mod2$predictions[, 2]
    tempdat$p.denominator <- 1 - mod2.pred
    
    if (!censor)
      tempdat$p.denominator[tempdat$exposure == 1] <- mod2.pred[tempdat$exposure == 1]
    
  } else if (continuous) {
    
    # Calculate numerator probabilities for continuous exposure
    if (is.null(numerator)) {
      mod1.pred <- mean(tempdat$exposure, na.rm = TRUE)
      mod1.sd <- sd(tempdat$exposure, na.rm = TRUE)
    } else {
      mod1 <- ranger(update.formula(numerator, formula(paste0(exposure," ~ ."))),
                     data = data, num.trees = num.trees, max.depth = max.depth,
                     num.threads = num.threads, respect.unordered.factors = "partition")
      mod1.pred <- mod1$predictions
      mod1.sd <- sd(tempdat$exposure - mod1.pred)
    }
    
    mod2 <- ranger(update.formula(denominator, formula(paste0(exposure," ~ ."))),
                   data = data, num.trees = num.trees, max.depth = max.depth,
                   num.threads = num.threads, respect.unordered.factors = "partition")
    mod2.pred <- mod2$predictions
    mod2.sd <- sd(tempdat$exposure - mod2.pred)
    
    # numerator
    a.num <- (tempdat$exposure - mod1.pred) / mod1.sd
    dens.num <- density(a.num)
    tempdat$p.numerator <- approx(x = dens.num$x, y = dens.num$y, xout = a.num)$y / mod1.sd
    
    # denominator
    a.denom <- (tempdat$exposure - mod2.pred) / mod2.sd
    dens.denom <- density(a.denom)
    tempdat$p.denominator <- approx(x = dens.denom$x, y = dens.denom$y, xout = a.denom)$y / mod2.sd
    
  }
  
  # Calculate inverse probability weights
  tempdat$ipw.weights <- tempdat$p.numerator / tempdat$p.denominator
  
  if (sum(is.na(tempdat$ipw.weights)) > 0)
    stop("NA's in weights!")
  
  if (!is.null(trunc)) {
    # Truncate weights if specified
    tempdat$weights.trunc <- tempdat$ipw.weights
    tempdat$weights.trunc[tempdat$ipw.weights < quantile(tempdat$ipw.weights, 0 + trunc)] <- quantile(tempdat$ipw.weights, 0 + trunc)
    tempdat$weights.trunc[tempdat$ipw.weights > quantile(tempdat$ipw.weights, 1 - trunc)] <- quantile(tempdat$ipw.weights, 1 - trunc)
  }
  
  if (is.null(numerator)) {
    if (is.null(trunc))
      return(list(ipw.weights = tempdat$ipw.weights, call = match.call(), den.mod = mod2))
    else
      return(list(ipw.weights = tempdat$ipw.weights, weights.trunc = tempdat$weights.trunc,
                  call = match.call(), den.mod = mod2))
  } else {
    if (is.null(trunc))
      return(list(ipw.weights = tempdat$ipw.weights, call = match.call(), num.mod = mod1, den.mod = mod2))
    else
      return(list(ipw.weights = tempdat$ipw.weights, weights.trunc = tempdat$weights.trunc,
                  call = match.call(), num.mod = mod1, den.mod = mod2))
  }
}

# Function for IPW-TM using XGBoost
ipwtm_xgb <- function(exposure, numerator = NULL, denominator, id, timevar, data, trunc = NULL,
                      continuous = FALSE, censor = FALSE, num.trees = 100, max.depth = 8, num.threads = 12, ...) {
  
  order.orig <- order(data[[id]], data[[timevar]])
  data <- data.table(data)[order.orig]
  tempdat <- data.table(
    id = data[[id]],
    timevar = data[[timevar]],
    exposure = data[[exposure]]
  )
  
  if (!continuous) {
    
    if (is.null(numerator)) {
      mod1.pred <- mean(tempdat$exposure, na.rm = TRUE)
    } else {
      mat <- model.matrix(numerator, data = data)
      mod1 <- xgboost(data = mat, label = tempdat$exposure, objective = "binary:logistic", 
                      nrounds = num.trees, max_depth = max.depth, nthread = num.threads)
      mod1.pred <- predict(mod1, newdata = mat)
    }
    
    tempdat$p.numerator <- 1 - mod1.pred
    
    if (!censor)
      tempdat$p.numerator[tempdat$exposure == 1] <- mod1.pred[tempdat$exposure == 1]
    
    mat <- model.matrix(denominator, data = data)
    mod2 <- xgboost(data = mat, label = tempdat$exposure, objective = "binary:logistic", 
                    nrounds = num.trees, max_depth = max.depth, nthread = num.threads)
    mod2.pred <- predict(mod2, newdata = mat)
    tempdat$p.denominator <- 1 - mod2.pred
    
    if (!censor)
      tempdat$p.denominator[tempdat$exposure == 1] <- mod2.pred[tempdat$exposure == 1]
    
  } else if (continuous) {
    
    if (is.null(numerator)) {
      mod1.pred <- mean(tempdat$exposure, na.rm = TRUE)
      mod1.sd <- sd(tempdat$exposure, na.rm = TRUE)
    } else {
      mat <- model.matrix(numerator, data = data)
      mod1 <- xgboost(data = mat, label = tempdat$exposure, objective = "reg:squarederror", 
                      nrounds = num.trees, max_depth = max.depth, nthread = num.threads)
      mod1.pred <- predict(mod1, newdata = mat)
      mod1.sd <- sd(tempdat$exposure - mod1.pred)
    }
    
    mat <- model.matrix(denominator, data = data)
    mod2 <- xgboost(data = mat, label = tempdat$exposure, objective = "reg:squarederror", 
                    nrounds = num.trees, max_depth = max.depth, nthread = num.threads)
    mod2.pred <- predict(mod2, newdata = mat)
    mod2.sd <- sd(tempdat$exposure - mod2.pred)
    
    # numerator
    a.num <- (tempdat$exposure - mod1.pred) / mod1.sd
    dens.num <- density(a.num)
    tempdat$p.numerator <- approx(x = dens.num$x, y = dens.num$y, xout = a.num)$y / mod1.sd
    
    # denominator
    a.denom <- (tempdat$exposure - mod2.pred) / mod2.sd
    dens.denom <- density(a.denom)
    tempdat$p.denominator <- approx(x = dens.denom$x, y = dens.denom$y, xout = a.denom)$y / mod2.sd
    
  }
  
  tempdat$ipw.weights <- tempdat$p.numerator / tempdat$p.denominator
  
  if (sum(is.na(tempdat$ipw.weights)) > 0)
    stop("NA's in weights!")
  
  if (!is.null(trunc)) {
    tempdat$weights.trunc <- tempdat$ipw.weights
    tempdat$weights.trunc[tempdat$ipw.weights < quantile(tempdat$ipw.weights, 0 + trunc)] <- quantile(tempdat$ipw.weights, 0 + trunc)
    tempdat$weights.trunc[tempdat$ipw.weights > quantile(tempdat$ipw.weights, 1 - trunc)] <- quantile(tempdat$ipw.weights, 1 - trunc)
  }
  
  if (is.null(numerator)) {
    if (is.null(trunc))
      return(list(ipw.weights = tempdat$ipw.weights, call = match.call(), den.mod = mod2))
    else
      return(list(ipw.weights = tempdat$ipw.weights, weights.trunc = tempdat$weights.trunc,
                  call = match.call(), den.mod = mod2))
  } else {
    if (is.null(trunc))
      return(list(ipw.weights = tempdat$ipw.weights, call = match.call(), num.mod = mod1, den.mod = mod2))
    else
      return(list(ipw.weights = tempdat$ipw.weights, weights.trunc = tempdat$weights.trunc,
                  call = match.call(), num.mod = mod1, den.mod = mod2))
  }
}

# Function for IPW-TM using GEE
ipwtm_gee <- function(exposure, numerator = NULL, denominator, id, timevar, data, trunc = NULL, ...) {
  
  order.orig <- order(data[[id]], data[[timevar]])
  data <- data.table(data)[order.orig]
  tempdat <- data.table(
    id = data[[id]],
    timevar = data[[timevar]],
    exposure = data[[exposure]]
  )
  
  # Calculate numerator probabilities using GEE
  if (is.null(numerator)) {
    fmla.num <- update.formula(numerator, formula(paste0(exposure," ~ 1")))
  } else {
    fmla.num <- update.formula(numerator, formula(paste0(exposure," ~ .")))
  }
  
  mod1 <- geeglm(formula = fmla.num, data = data, id = id, corstr = "ar1", ...)
  mod1.pred <- as.numeric(predict(mod1))
  mod1.sd <- sqrt(as.numeric(summary(mod1)$dispersion[1]))
  
  a.num <- (tempdat$exposure - mod1.pred) / mod1.sd
  dens.num <- density(a.num)
  tempdat$p.numerator <- approx(x = dens.num$x, y = dens.num$y, xout = a.num)$y / mod1.sd
  
  # Calculate denominator probabilities using GEE
  mod2 <- geeglm(update.formula(denominator, formula(paste0(exposure," ~ ."))),
                 data = data, id = id, corstr = "ar1", ...)
  mod2.pred <- as.numeric(predict(mod2))
  mod2.sd <- sqrt(as.numeric(summary(mod2)$dispersion[1]))
  
  a.denom <- (tempdat$exposure - mod2.pred) / mod2.sd
  dens.denom <- density(a.denom)
  tempdat$p.denominator <- approx(x = dens.denom$x, y = dens.denom$y, xout = a.denom)$y / mod2.sd
  
  tempdat$ipw.weights <- tempdat$p.numerator / tempdat$p.denominator
  
  if (sum(is.na(tempdat$ipw.weights)) > 0)
    stop("NA's in weights!")
  
  if (!is.null(trunc)) {
    tempdat$weights.trunc <- tempdat$ipw.weights
    tempdat$weights.trunc[tempdat$ipw.weights < quantile(tempdat$ipw.weights, 0 + trunc)] <- quantile(tempdat$ipw.weights, 0 + trunc)
    tempdat$weights.trunc[tempdat$ipw.weights > quantile(tempdat$ipw.weights, 1 - trunc)] <- quantile(tempdat$ipw.weights, 1 - trunc)
  }
  
  if (is.null(numerator)) {
    if (is.null(trunc))
      return(list(ipw.weights = tempdat$ipw.weights, call = match.call(), den.mod = mod2))
    else
      return(list(ipw.weights = tempdat$ipw.weights, weights.trunc = tempdat$weights.trunc,
                  call = match.call(), den.mod = mod2))
  } else {
    if (is.null(trunc))
      return(list(ipw.weights = tempdat$ipw.weights, call = match.call(), num.mod = mod1, den.mod = mod2))
    else
      return(list(ipw.weights = tempdat$ipw.weights, weights.trunc = tempdat$weights.trunc,
                  call = match.call(), num.mod = mod1, den.mod = mod2))
  }
}
