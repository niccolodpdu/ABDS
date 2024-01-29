#' @name eCOT
#' @title Extended Cosine-based One-sample Test
#' @description This function computes cosine values to select down-regulated signature genes (DSGs).
#' @param input The input matrix with features on rows and samples on columns. The input matrix should have only one sample in each 
#'     group, or use 'super_sample' function first to create a super-sample that only has one sample in each group.
#' @param cos.thres The cosine threshold for DSGs.
#' @param top The upper bound of the number of total DSGs.
#' @param per The upper bound of the number of DSGs of each group.
#' @export result A matrix of: the highest cosine values with the reference vectors (column 1), group indexes for the DSG (column 2), and their feature names (column 3) 
#' @examples
#' output<-eCOT(input)
#' output<-eCOT(input,cos.thres=0.99)
#' output<-eCOT(input,top=100)
#' output<-eCOT(input,per=20)
#' # If you have multiple samples in one group, please use the 'super_sample' function to shrink the group size to one sample per group.
#' # See the end of this script.
#' # Please use only one of the filtering methods ('cos.thres', 'top', 'per') at one time.

eCOT<-function(input, cos.thres = NULL, top = NULL, per = NULL){
  
  if (is.null(row.names(input))){
    row.names(input)<-1:dim(input)[1]
  }
  
  result<-cos_sDEG(input)
  
  if (!is.null(cos.thres)){
    result<-result[(result[,1]>=cos.thres),]
  }
  
  if (!is.null(top)){
    result<-result[order(result[,1],decreasing = TRUE),]
    result<-result[1:min(dim(result)[1],top),]
  }
  
  if (!is.null(per)){
    temp<-NULL
    output<-NULL
    result<-result[order(result[,2],result[,1],decreasing = TRUE),]
    for (i in 1:max(result[,2])){
      temp<-result[result[,2]==i,]
      temp<-temp[1:min(dim(temp)[1],per),]
      output<-rbind(output,temp)
    }
    result<-output
  }
  
  
  return(result)
}


#### Helper Functions ####

#' Calculate the cosine between vector A and vector B
mycosine <- function(A, B) {
  return(sum(A * B) / sqrt(sum(A ^ 2) * sum(B ^ 2)))
}

#' Calculate the cosine between a data matrix (features on rows and samples on columns)
## Applied to data matrix with only one sample in each group, or use 'super_sample' function first.
cos_sDEG <- function(data) {
  gene <- dim(data)[1]
  smp <- dim(data)[2]
  
  ref <- vector("list", length = smp)
  for (i in 1:smp) {
    ref[[i]] <- rep(1, smp)
    ref[[i]][i] <- 0
  }
  
  cosine_value <- array(0, dim = gene)
  which_group <- array(0, dim = gene)
  
  for (i in 1:gene) {
    temp <- NULL
    for (j in 1:smp) {
      temp <- c(temp, mycosine(ref[[j]], data[i, ]))
    }
    cosine_value[i] <- max(temp)
    which_group[i] <- which.max(temp)
  }
  
  row_name<-row.names(data)
  
  return(cbind(cosine_value,which_group,row_name))
}

#' Super-sample, if you have several samples in one group & reduce them into one sample per group.
#' @param Input Input matrix, with features on rows and samples on columns.
#' @param nRep Number of samples in each group. E.g., Group 1 has 4 samples, Group 2 has 3 samples, Group 3 has 3 samples.
#' Then nRep = c(4,3,3)
super_sample<-function(input,nRep){
  output <- NULL
  start <- 1
  for(rep in nRep){
    output <- cbind(output, rowMeans(input[, start : (start + rep - 1)])) 
    start <- start + rep
  }
  return(output)
}
