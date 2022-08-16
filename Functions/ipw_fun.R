
ipwtm_ranger <- function(exposure, numerator = NULL, denominator, id, timevar, data, trunc = NULL,
                         continuous = FALSE, death = FALSE, num.trees = 100, max.depth = 8, num.threads = 12, ...) {
  
  tempcall <- match.call()
  
  if (!("exposure" %in% names(tempcall))) 
    stop("No exposure variable specified")
  if (!("denominator" %in% names(tempcall))) 
    stop("No denominator model specified")
  if (!is.null(tempcall$numerator) & !is(eval(tempcall$numerator), "formula")) 
    stop("Invalid numerator formula specified")
  if (!is.null(tempcall$denominator) & !is(eval(tempcall$denominator), "formula")) 
    stop("Invalid denominator formula specified")
  if (!("id" %in% names(tempcall))) 
    stop("No patient id specified")
  if (!("timevar" %in% names(tempcall))) 
    stop("No timevar specified")
  if (!("data" %in% names(tempcall))) 
    stop("No data specified")
  if (continuous == TRUE & death == TRUE)
    stop("How do you have a continuous death?")
  if (!is.null(tempcall$trunc)) {
    if (tempcall$trunc < 0 | tempcall$trunc > 0.5) 
      stop("Invalid truncation percentage specified (0-0.5)")
  }
  
  order.orig <- 1:nrow(data)
  order.orig <- order.orig[order(eval(parse(text = paste("data$", deparse(id, width.cutoff = 500), sep = ""))), 
                                 eval(parse(text = paste("data$", deparse(timevar, width.cutoff = 500), sep = ""))))]
  data <- data[order(eval(parse(text = paste("data$", deparse(id, width.cutoff = 500), sep = ""))), 
                     eval(parse(text = paste("data$", deparse(timevar, width.cutoff = 500), sep = "")))), ]
  tempdat <- data.frame(id = data[, as.character(id)],
                        timevar = data[, as.character(timevar)], 
                        exposure = data[, as.character(exposure)])
  
  if (!continuous){
  
    if (is.null(tempcall$numerator)) {
      
      tempdat$p.numerator <- 1 - mean(tempdat$exposure, na.rm = TRUE)
      
      if (!death)
        tempdat$p.numerator[tempdat$exposure == 1] <- mean(tempdat$exposure, na.rm = TRUE)
      
    } else {
      
      mod1 <- ranger(update.formula(numerator, formula(paste0("factor(",exposure,") ~ ."))), probability = TRUE,
                     data = data, num.trees = num.trees, max.depth = max.depth, num.threads = num.threads,
                     respect.unordered.factors = "partition")
      mod1.pred <- mod1$predictions[,2]
      tempdat$p.numerator <- 1 - c(mod1.pred)
      
      if (!death)
        tempdat$p.numerator[tempdat$exposure == 1] <- c(mod1.pred)[tempdat$exposure == 1]
      
    }
    
    mod2 <- ranger(update.formula(denominator, formula(paste0("factor(",exposure,") ~ ."))), probability = TRUE,
                   data = data, num.trees = num.trees, max.depth = max.depth, num.threads = num.threads,
                   respect.unordered.factors = "partition")
    mod2.pred <- mod2$predictions[,2]
    tempdat$p.denominator <- 1 - c(mod2.pred)
    
    if (!death)
      tempdat$p.denominator[tempdat$exposure == 1] <- c(mod2.pred)[tempdat$exposure == 1]
    
    tempdat$ipw.weights <- with(tempdat, p.numerator/p.denominator)

    
  } else if (continuous) {
    
    if (is.null(tempcall$numerator)) {
      
      tempdat$p.numerator <- dnorm(tempdat$exposure, mean(tempdat$exposure, na.rm = TRUE), 
                                   sd(tempdat$exposure, na.rm = TRUE))
      
    } else {
      
      mod1 <- ranger(update.formula(numerator, formula(paste0(exposure," ~ ."))),
                     data = data, num.trees = num.trees, max.depth = max.depth, 
                     num.threads = num.threads, respect.unordered.factors = "partition")
      mod1.pred <- mod1$predictions
      mod1.sd <- sd(tempdat$exposure - mod1$predictions)
      tempdat$p.numerator <- dnorm(tempdat$exposure, mod1.pred, mod1.sd)
      
    }
    
    mod2 <- ranger(update.formula(denominator, formula(paste0(exposure," ~ ."))),
                   data = data, num.trees = num.trees, max.depth = max.depth, 
                   num.threads = num.threads,respect.unordered.factors = "partition")
    mod2.pred <- mod2$predictions
    mod2.sd <- sd(tempdat$exposure - mod2$predictions)
    tempdat$p.denominator <- dnorm(tempdat$exposure, mod2.pred, mod2.sd)
    
    tempdat$ipw.weights <- with(tempdat, p.numerator/p.denominator)
    
  }
  
  if (sum(is.na(tempdat$ipw.weights)) > 0) 
    stop("NA's in weights!")
  
  if (!(is.null(tempcall$trunc))) {
    tempdat$weights.trunc <- tempdat$ipw.weights
    tempdat$weights.trunc[tempdat$ipw.weights <= quantile(tempdat$ipw.weights, 0 + trunc)] <- quantile(tempdat$ipw.weights, 0 + trunc)
    tempdat$weights.trunc[tempdat$ipw.weights > quantile(tempdat$ipw.weights, 1 - trunc)] <- quantile(tempdat$ipw.weights, 1 - trunc)
  }
  
  if (is.null(tempcall$trunc)) {
    
    if (is.null(tempcall$numerator)) 
      return(list(ipw.weights = tempdat$ipw.weights[order(order.orig)], 
                  call = tempcall, den.mod = mod2))
    else 
      return(list(ipw.weights = tempdat$ipw.weights[order(order.orig)], 
                  call = tempcall, num.mod = mod1, den.mod = mod2))
    
  } else {
    
    if (is.null(tempcall$numerator)) 
      return(list(ipw.weights = tempdat$ipw.weights[order(order.orig)], 
                  weights.trunc = tempdat$weights.trunc[order(order.orig)], 
                  call = tempcall, den.mod = mod2))
    else 
      return(list(ipw.weights = tempdat$ipw.weights[order(order.orig)], 
                  weights.trunc = tempdat$weights.trunc[order(order.orig)], 
                  call = tempcall, num.mod = mod1, den.mod = mod2))
    
  }
  
}

ipwtm_gee <- function(exposure, numerator = NULL, denominator, id, timevar, data, trunc = NULL, ...) {
  
  tempcall <- match.call()
  
  if (!("exposure" %in% names(tempcall))) 
    stop("No exposure variable specified")
  if (!("denominator" %in% names(tempcall))) 
    stop("No denominator model specified")
  if (!is.null(tempcall$numerator) & !is(eval(tempcall$numerator), "formula")) 
    stop("Invalid numerator formula specified")
  if (!is.null(tempcall$denominator) & !is(eval(tempcall$denominator), "formula")) 
    stop("Invalid denominator formula specified")
  if (!("id" %in% names(tempcall))) 
    stop("No patient id specified")
  if (!("timevar" %in% names(tempcall))) 
    stop("No timevar specified")
  if (!("data" %in% names(tempcall))) 
    stop("No data specified")
  if (!is.null(tempcall$trunc)) {
    if (tempcall$trunc < 0 | tempcall$trunc > 0.5) 
      stop("Invalid truncation percentage specified (0-0.5)")
  }
  
  order.orig <- 1:nrow(data)
  order.orig <- order.orig[order(eval(parse(text = paste("data$", deparse(id, width.cutoff = 500), sep = ""))), 
                                 eval(parse(text = paste("data$", deparse(timevar, width.cutoff = 500), sep = ""))))]
  data <- data[order(eval(parse(text = paste("data$", deparse(id, width.cutoff = 500), sep = ""))), 
                     eval(parse(text = paste("data$", deparse(timevar, width.cutoff = 500), sep = "")))), ]
  tempdat <- data.frame(id = data[, as.character(id)],
                        timevar = data[, as.character(timevar)], 
                        exposure = data[, as.character(exposure)])
  
  if (is.null(tempcall$numerator)) {
    
    tempdat$p.numerator <- dnorm(tempdat$exposure, 
                                 mean(tempdat$exposure, na.rm = TRUE),
                                 sd(tempdat$exposure, na.rm = TRUE))
    
  } else {
    
    mod1 <- geeglm(formula = update.formula(numerator, formula(paste0(exposure," ~ ."))), data = data, 
                   id = zip, corstr = "ar1", waves = ssn_time,  ...)
    mod1.pred <- as.numeric(predict(mod1))
    mod1.var <- as.numeric(summary(mod1)$dispersion[1])
    tempdat$p.numerator <- dnorm(tempdat$exposure, c(mod1.pred), sqrt(mod1.var))
    
  }
  
  mod2 <- geeglm(update.formula(denominator, formula(paste0(exposure," ~ ."))), 
                 data = data, id = zip, corstr = "ar1", waves = ssn_time,  ...)
  mod2.pred <- as.numeric(predict(mod2)) 
  mod2.var <- as.numeric(summary(mod2)$dispersion[1])
  tempdat$p.denominator <- dnorm(tempdat$exposure, c(mod2.pred), sqrt(mod2.var))
  
  tempdat$ipw.weights <- with(tempdat, p.numerator/p.denominator)
  
  if (sum(is.na(tempdat$ipw.weights)) > 0) 
    stop("NA's in weights!")
  
  if (!(is.null(tempcall$trunc))) {
    tempdat$weights.trunc <- tempdat$ipw.weights
    tempdat$weights.trunc[tempdat$ipw.weights <= quantile(tempdat$ipw.weights, 0 + trunc)] <- quantile(tempdat$ipw.weights, 0 + trunc)
    tempdat$weights.trunc[tempdat$ipw.weights > quantile(tempdat$ipw.weights, 1 - trunc)] <- quantile(tempdat$ipw.weights, 1 - trunc)
  }
  
  if (is.null(tempcall$trunc)) {
    
    if (is.null(tempcall$numerator)) 
      return(list(ipw.weights = tempdat$ipw.weights[order(order.orig)], 
                  call = tempcall, den.mod = mod2))
    else 
      return(list(ipw.weights = tempdat$ipw.weights[order(order.orig)], 
                  call = tempcall, num.mod = mod1, den.mod = mod2))
    
  } else {
    
    if (is.null(tempcall$numerator)) 
      return(list(ipw.weights = tempdat$ipw.weights[order(order.orig)], 
                  weights.trunc = tempdat$weights.trunc[order(order.orig)], 
                  call = tempcall, den.mod = mod2))
    else 
      return(list(ipw.weights = tempdat$ipw.weights[order(order.orig)], 
                  weights.trunc = tempdat$weights.trunc[order(order.orig)], 
                  call = tempcall, num.mod = mod1, den.mod = mod2))
    
  }
  
}

