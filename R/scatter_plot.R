#' @name scatter_plot
#' @title A function for illustrating data points in scatter plot with optional marker color-coding
#' @param data With samples on the columns and features on the rows.
#' @param mg.info Do you have the info for marker genes?
#' @param mg A list of marker genes
#' @param mg.col A vector of colors


scatter_plot<-function(data, mg.info = FALSE, mg = NULL, 
                       mg.col = c('red','orange','dodgerblue','green4','purple')){
  
  library(corpcor)
  
  X <- data
  Xproj <- X
  
  Xproj <- apply(Xproj, 1, function(x) x / sum(x))

  #colnames(Xproj) <- data[,1]
  A <- diag(dim(Xproj)[1])
  
  K = dim(A)[2]
  
  PS <- t(matrix(c(cos((seq_len(K) - 1) * 2 * pi / K),
                   sin((seq_len(K) - 1) * 2 * pi / K)), K))
  

  tmp <- PS %*% pseudoinverse(A)
  #tmp[1,] <- tmp[1,] / c(sqrt(tmp[1,] %*% tmp[1,]))
  #tmp[2,] <- tmp[2,] - c(tmp[2,] %*% tmp[1,]) * tmp[1,]
  #tmp[2,] <- tmp[2,] / c(sqrt(tmp[2,] %*% tmp[2,]))
  Xp <- tmp %*% Xproj
  
  plot(t(PS), xlab = NA, ylab = NA, asp = 1, pch = '.')
  #plot(Xp[1,], Xp[2,], col = rgb(200, 200, 200, max = 255, alpha = 125),
  #     xlab = NA, ylab = NA, asp = 1)
  points(Xp[1,], Xp[2,], col = rgb(0, 0, 0, max = 255, alpha = 125),
         xlab = NA, ylab = NA, asp = 1)
  
  #### If you want to add some marker genes, here is your place ####
  
  if (mg.info == TRUE) {
    for(i in seq_along(mg)){
      points(Xp[1,mg[[i]]], Xp[2,mg[[i]]], col=mg.col[i])
    }
  }

  
  points(t(PS), pch = 17)
}