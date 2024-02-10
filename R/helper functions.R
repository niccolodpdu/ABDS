#' Calculate the cosine between vector A and vector B
mycosine <- function(A, B) {
  return(sum(A * B) / sqrt(sum(A ^ 2) * sum(B ^ 2)))
}

#' Calculate the cosine between data and reference vector
cos_ref <- function(data, ref) {
  gene <- dim(data)[1]
  smp <- dim(data)[2]
  
  res <- array(0, dim = gene)
  for (i in 1:gene) {
    res[i] <- mycosine(ref, data[i, ])
  }
  return(res)
}

#' Calculate the cosine between data and CEG reference, e.g. c(1, 1, 1)
cos_iCEG <- function(data) {
  gene <- dim(data)[1]
  smp <- dim(data)[2]
  
  ref <- rep(1, smp)
  
  res <- array(0, dim = gene)
  for (i in 1:gene) {
    res[i] <- mycosine(ref, data[i, ])
  }
  return(res)
}

#' Calculate the cosine between data and closest iDEG reference, e.g. c(1, 0, 0)
cos_iDEG <- function(data) {
  gene <- dim(data)[1]
  smp <- dim(data)[2]
  
  ref <- vector("list", length = smp)
  for (i in 1:smp) {
    ref[[i]] <- rep(0, smp)
    ref[[i]][i] <- 1
  }
  
  res <- array(0, dim = gene)
  
  for (i in 1:gene) {
    temp <- NULL
    for (j in 1:smp) {
      temp <- c(temp, mycosine(ref[[j]], data[i, ]))
    }
    res[i] <- max(temp)
  }
  return(res)
}

#' total count normalization
totalcount <- function(data) {
  scalar <- colSums(data) / mean(colSums(data))
  data <- apply(data, 2, function(x)
    x / sum(x)) * mean(colSums(data))
  
  return(list(data = data, norm_factor = scalar))
}


#' This function computes cosine values to select markers.
cotMG <- function(data=NULL, Sest, thres.low=0.05, thres.high=1, cos.thres=1, 
                  top=NULL, per=NULL) {
  
  if (thres.low >= thres.high) {
    stop("thres.low must be smaller than thres.high!")
  }
  Valid <- NULL
  nonzero.ind <- NULL
  if (!is.null(data)) {
    if (is(data, "data.frame")) {
      data <- as.matrix(data)
    } else if (is(data, "SummarizedExperiment")) {
      data <- SummarizedExperiment::assay(data)
    } else if (is(data, "ExpressionSet")) {
      data <- Biobase::exprs(data)
    } else if (is(data, "matrix") == FALSE) {
      stop("Only matrix, dataframe, SummarizedExperiment and ExpressionSet
            object are supported for expression data!")
    }
    if (nrow(data) != nrow(Sest)){
      stop("The row number of data and S matrix should be the same!")
    }
    if (sum(is.na(data)) > 0) {
      stop("Data with missing values are not supported!")
    }
    if (sum(data<0) > 0) {
      stop("Only non-negative data are supported!")
    }
    if (is.null(rownames(data))) {
      rownames(data) <- seq_len(nrow(data))
    }
    nonzero.ind <- rowSums(data) > 0
    data <- data[nonzero.ind, ]
    X <- t(data)
    
    sigNorm <- apply(X, 2, function(x) norm(matrix(x),"F") )
    Valid <- sigNorm >= quantile(sigNorm, thres.low) &
      sigNorm <= quantile(sigNorm, thres.high)
    X <- X[, Valid]
  }
  
  if (is(Sest, "data.frame")) {
    Sest <- as.matrix(Sest)
  } else if (is(Sest, "matrix") == FALSE) {
    stop("Only matrix and data frame are supported for S matrix!")
  }
  if (sum(is.na(Sest)) > 0) {
    stop("S matrix with missing values are not supported!")
  }
  if (sum(Sest<0) > 0) {
    stop("Only non-negative S matrix are supported!")
  }
  if (!is.null(nonzero.ind)){
    Sest <- Sest[nonzero.ind, ]
  }
  if (is.null(rownames(Sest))) {
    rownames(Sest) <- seq_len(nrow(Sest))
  }
  S <- t(Sest)
  if (is.null(Valid)) {
    sigNorm <- apply(S, 2, function(x) norm(matrix(x),"F") )
    Valid <- sigNorm >= quantile(sigNorm, thres.low) &
      sigNorm <= quantile(sigNorm, thres.high)
  }
  S <- S[, Valid]
  S.norm <- apply(S, 2, function(x) x / norm(matrix(x),"2"))
  mg.cos <- matrix(NA, ncol(S), 2)
  rownames(mg.cos) <- colnames(S)
  colnames(mg.cos) <- c('cos', 'source')
  for (i in seq_len(ncol(S))){
    mg.cos[i, 1] <- max(S.norm[, i])
    mg.cos[i, 2] <- which.max(S.norm[, i])
  }
  
  if (!is.null(cos.thres)){
    mg.cos <- mg.cos[mg.cos[, 1] >= cos.thres, ]
  }
  
  mg.cos.sorted <- mg.cos[sort.int(mg.cos[, 1], decreasing=TRUE, 
                                   index.return=TRUE)$ix, ]
  if (!is.null(top)){
    mg.cos.sorted <- mg.cos.sorted[1:top, ]
  }
  
  mg.list <- vector("list", nrow(S))
  for (i in seq_along(mg.list)){
    mg.list[[i]] <- rownames(mg.cos.sorted[mg.cos.sorted[, 2] == i, ])
  }
  
  if (!is.null(per)){
    for (i in seq_along(mg.list)){
      if (length(mg.list[[i]]) > per){
        mg.list[[i]] <- mg.list[[i]][1:per]
      }
    }
  }
  
  return(list(mg.list = mg.list, mg.cos = mg.cos))
}

