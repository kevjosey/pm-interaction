dlmm_pen <- function(basis, smooth = TRUE, shrink = TRUE)
{
  x_names <- attributes(basis)$x_names
  M <- length(x_names)
  # Create basis
  seqlag <- function(range) seq(range[1], range[2], by = 1)
  b <- list()
  smoothS <- list()
  shrinkS <- list()
  for (i in 1:M) {
    n <- x_names[i]
    b[[n]] <- do.call(onebasis, attributes(basis)$basis[[n]])
    smoothS[[n]] <- attributes(b[[n]])$S
    e <- eigen(smoothS[[n]], symmetric = TRUE)
    # if (shrink & any(e$values < 1 / 10))
    #   e[which(e$values < 1 / 10)] <- 1 / 10
    smoothS[[n]] <- smoothS[[n]] / e$values[1]#e$vectors %*% diag(e$values) %*% t(e$vectors)
    shrinkS[[n]] <- tcrossprod(e$vectors[, which(e$values < 1/10)])
  }
  # Create penalization matrices
  pen_smooth <- lapply(1:M, function(i) matrix(0, ncol(basis), ncol(basis)))
  pen_shrink <- lapply(1:M, function(i) matrix(0, ncol(basis), ncol(basis)))
  for (i in 1:M) {
    n1 <- x_names[i]
    idx <- grep(paste0("e", i, "\\.[0-9]{1,2}"), colnames(basis))
    pen_smooth[[i]][idx, idx] <- smoothS[[n1]]
    pen_shrink[[i]][idx, idx] <- shrinkS[[n1]]

    for (j in i:M) {
      n2 <- x_names[j]
      idx <- grep(paste0("m", i, "\\.[0-9]{1,2}\\.m", j, "\\.[0-9]{1,2}"), colnames(basis))
      if (length(idx) > 0) {
        pen_smooth[[i]][idx, idx] <- smoothS[[n1]] %x% diag(ncol(b[[n2]]))
        pen_shrink[[i]][idx, idx] <- shrinkS[[n1]] %x% diag(ncol(b[[n2]]))
        pen_smooth[[j]][idx, idx] <- diag(ncol(b[[n1]])) %x% smoothS[[n2]]
        pen_shrink[[j]][idx, idx] <- diag(ncol(b[[n1]])) %x% shrinkS[[n2]]
      }
    }
  }

  if (smooth & shrink) {
    return(c(pen_smooth, pen_shrink))
  } else if (smooth) {
    return(pen_smooth)
  } else if (shrink) {
    return(pen_shrink)
  }
}
