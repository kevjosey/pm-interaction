dlmm_pred <- function(basis, model = NULL, coef = NULL, vcov = NULL,
                      xscale = c(), lag = list(), bylag = c(),
                      marginal = TRUE, xmarg = c(),
                      ci_level = 0.95)
{
  if (!any(class(basis) %in% "dlmm_basis")) stop("basis must be from dlmm_basis")
  var_name <- deparse(substitute(basis))
  z <- qnorm(1 - (1 - ci_level) / 2)

  x_names <- attributes(basis)$x_names
  M <- length(x_names)


  # Create basis
  seqlag <- function(range, bylag) seq(range[1], range[2], by = bylag)
  b <- list()
  len <- rep(0, M)
  lag_vals <- list()
  for (i in 1:M) {
    n <- x_names[i]
    attributes(basis)$basis[[n]]$x <- seqlag(attributes(basis)$lag[[n]], bylag[i])
    b[[n]] <- do.call(onebasis, attributes(basis)$basis[[n]])
    len[i] <- length(attributes(basis)$basis[[n]]$x)
    lag_vals[[n]] <- attributes(basis)$basis[[n]]$x
  }
  # Create crossbasis
  if (M > 1) {
    cb <- list()
    for (i in 1:M) {
      for (j in i:M) {
        n1 <- x_names[i]
        n2 <- x_names[j]
        n <- paste0(n1, "-", n2)
        cb[[n]] <- b[[n2]] %x% b[[n2]]
      }
    }
  }

  # Get coefficients and covariance matrix
  out <- list(ci_level = ci_level, lag_vals = lag_vals, x_names = x_names)
  if (!is.null(model)) {
    coef <- coef(model)
    vcov <- vcov(model)
  }
# stop()
  # Calculate fitted
  out$main <- list()
  out$interaction <- list()
  out$marginal <- list()
  for (i in 1:M) {
    n <- x_names[i]
    idx <- grep(paste0(var_name, "e", i, "\\.[0-9]{1,2}"), names(coef))
    out$main[[n]] <- list(coef = coef[idx], vcov = vcov[idx, idx],
                          fit = b[[n]] %*% coef[idx] * xscale[i],
                          se = sqrt(pmax(0, diag(b[[n]] %*% (vcov[idx, idx] %*% t(b[[n]]))))) * xscale[i])
    out$main[[n]]$low <- out$main[[n]]$fit - z * out$main[[n]]$se
    out$main[[n]]$high <- out$main[[n]]$fit + z * out$main[[n]]$se
    out$marginal[[n]] <- out$main[[n]]
  }

  if (M > 1) {
    for (i in 1:M) {
      for (j in i:M) {
        n1 <- x_names[i]
        n2 <- x_names[j]
        n <- paste0(n1, "-", n2)
        idx <- grep(paste0(var_name, "m", i, "\\.[0-9]{1,2}\\.m", j, "\\.[0-9]{1,2}"), names(coef))
        if (length(idx) > 0) {
          out$interaction[[n]] <- list(coef = coef[idx], vcov = vcov[idx, idx])
          out$interaction[[n]]$fit <-
            matrix(cb[[n]] %*% out$interaction[[n]]$coef * xscale[i] * xscale[j], len[j], len[i])
          out$interaction[[n]]$se <-
            matrix(pmax(0, sqrt(diag(cb[[n]] %*% (vcov[idx, idx] %*% t(cb[[n]])))) * xscale[i] * xscale[j]),
                   len[j], len[i])

          out$marginal[[n1]]$fit <- out$marginal[[n1]]$fit + colSums(out$interaction[[n]]$fit) * xmarg[j]
          out$marginal[[n2]]$fit <- out$marginal[[n2]]$fit + rowSums(out$interaction[[n]]$fit) * xmarg[i]
          out$marginal[[n1]]$se <- out$marginal[[n1]]$se + colSums(out$interaction[[n]]$se)
          out$marginal[[n2]]$se <- out$marginal[[n2]]$se + rowSums(out$interaction[[n]]$se)
          out$marginal[[n1]]$low <- out$marginal[[n1]]$fit - z * out$marginal[[n1]]$se
          out$marginal[[n1]]$high <- out$marginal[[n1]]$fit + z * out$marginal[[n1]]$se
          out$marginal[[n2]]$low <- out$marginal[[n2]]$fit - z * out$marginal[[n2]]$se
          out$marginal[[n2]]$high <- out$marginal[[n2]]$fit + z * out$marginal[[n2]]$se
        }
      }
    }
  }


  class(out) <- "dlmm_pred"
  return(out)
}
