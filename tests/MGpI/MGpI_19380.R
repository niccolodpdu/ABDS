library(Ternary)
library(hrbrthemes)
library(ggplot2)
library(DirichletReg)
library(DescTools)
library(psych)
library(scales)
library(devEMF)
library(rgl)
library(magick)
library(ggtern)
library(plyr)
library(directlabels)
library(ggplot2) 
library(devEMF)
library(nnls)
library(pcaMethods) 
library(lsa)
library(tilting)
library(limSolve)
library(scImpute)
library(GEOquery)
source('imputation/kNN.R')
source('imputation/gbmImpute.R')
source('imputation/lmImpute.R')
source('imputation/meanImpute.R')
source('imputation/SVD.R')
source('imputation/SVT.R')
source('imputation/SVTApprox.R')
source('imputation/tsImpute.R')
source('imputation/utils.R')
source('imputation/imputation_functions_wrappers.R')
source('MGpI.R')
source('helper functions.R')

gsm <- getGEO("GSE19380")
pheno <- pData(phenoData (gsm[[1]]))
data <- 2^exprs(gsm[[1]])

test_function<-function(data,rate,Total_missing,missing_rate_threshold=0.6,
                        MAR_option = 'All',pram = 12,sd_scaler=4,nRep=c(4,4,4,4)){
  
  # pram: parameter for SVD, PPCA, NIPALS, KNN
  # rate: MAR rate
  ################# Read-in Data #################
  data_grnd_full <- data [,1:16]
  
  ################# Introducing Missing Values #################
  
  MAR=rate*Total_missing
  
  ### MAR
  removed<-data_grnd_full
  
  if (MAR_option == 'SG only'){
    index<-sample(which(data_grnd_full[1:(dim(data_grnd_full)[2]*dim(data_grnd_full)[1]),]>0),MAR)
  } else {
    index<-sample(which(data_grnd_full>0),MAR)
  }
  
  for (i in 1:MAR){
    removed[index]<-0
  }
  sum(removed == 0)  
  
  ### LLOD
  removed[removed < quantile(removed,(Total_missing/(dim(data_grnd_full)[1]*dim(data_grnd_full)[2])))]<-0
  sum(removed == 0)  
  
  
  input_data <- removed
  before_norm<-input_data
  after_norm<-before_norm
  before_remove<-input_data
  after_remove<-NULL
  
  ## Mask-out (inactive)
  missing_rate<-apply(after_norm,1,function(x)
    sum(x==0)/length(x))
  
  # after_remove<-after_norm[missing_rate<missing_rate_threshold,]
  # dim(after_remove)
  # Please note that due to randomness, number of masked-out genes may not be identical in each experiment!
  
  ## Removing features with all missing value
  after_norm<-after_norm[missing_rate<1,]
  data_grnd_full<-data_grnd_full[missing_rate<1,]
  
  ################# Imputation & Evaluation #################
  
  ## Some methods only recognize NA instead of 0
  after_norm_na<-after_norm
  after_norm_na[after_norm_na==0]<-NA
  after_norm_na<-log(after_norm_na,2)
  
  ## Functions for evaluation
  rmse <- function(ground_truth,imputed) sqrt(mean((log(ground_truth,2) - log(imputed,2))^2))
  nrmse <- function(ground_truth,imputed) rmse(ground_truth,imputed)/sd(log(ground_truth))
  euclidean <- function(a, b) sqrt(sum((a - b)^2))
  # specifically for GSE110554 where data is distributed within [0,1), we use original space for RMSE and NRMSE instead of log
  
  nmiss<-sum(any(apply(after_norm_na,1,function(x) is.na(x))))
  
  ####  Min/2 Imputed  ####
  Half_min<-min(after_norm[after_norm!=0])/2
  Half_min_imputed<-after_norm
  Half_min_imputed[Half_min_imputed==0]<-Half_min

  rmse_Half_min_all <- rmse(data_grnd_full,
                            Half_min_imputed)
  

  nrmse_Half_min_all <- nrmse(data_grnd_full,
                              Half_min_imputed)
  
  
  ####  Mean Imputed  ####
  Mean<-mean(after_norm_na[!is.na(after_norm_na)])
  Mean_imputed<-after_norm_na
  Mean_imputed[is.na(after_norm_na)]<-Mean
  Mean_imputed<-2^Mean_imputed

  rmse_Mean_all <- rmse(data_grnd_full,
                        Mean_imputed)
  
  nrmse_Mean_all <- nrmse(data_grnd_full,
                          Mean_imputed)
  
  
  
  ####  ABDS Imputed  ####
  ABDS_imputed<-MGpI(after_norm,nRep = nRep,'global',missing_rate_threshold=1, sd_scaler=sd_scaler)
  rmse_ABDS_all <- rmse(data_grnd_full,
                        ABDS_imputed)
  

  nrmse_ABDS_all <- nrmse(data_grnd_full,
                          ABDS_imputed)
  
  
  
  
  ####  swKNN Imputed ####
  
  swKNN_imputed<-KNN(t(after_norm_na),pram)
  swKNN_imputed<-2^swKNN_imputed

  rmse_swKNN_all <- rmse(data_grnd_full,
                         swKNN_imputed)
  
  nrmse_swKNN_all <- nrmse(data_grnd_full,
                           swKNN_imputed)
  
  ####  NIPALS Imputed ####
  
  NIPALS_imputed<-completeObs(pca(
    object = after_norm_na,
    method = "nipals",
    nPcs = pram,
    center = TRUE,
    scale = "uv"
  ))
  
  NIPALS_imputed<-2^NIPALS_imputed
  rmse_NIPALS_all <- rmse(data_grnd_full,
                          NIPALS_imputed)
  
  nrmse_NIPALS_all <- nrmse(data_grnd_full,
                            NIPALS_imputed)
  
  
  
  
  ####  SVD Imputed ####
  
  SVD_imputed<-completeObs(pca(
    object = after_norm_na,
    method = "svdImpute",
    nPcs = pram,
    center = FALSE,
    scale = "none"
  ))
  
  SVD_imputed<-2^SVD_imputed

  rmse_SVD_all <- rmse(data_grnd_full,
                       SVD_imputed)

  nrmse_SVD_all <- nrmse(data_grnd_full,
                         SVD_imputed)
  
  
  ####  SVT Imputed ####
  
  SVT_imputed<-SVTImpute(after_norm_na, lambda = pram,
                         max.iters = 3000)$x
  
  SVT_imputed<-2^SVT_imputed

  rmse_SVT_all <- rmse(data_grnd_full,
                       SVT_imputed)
  
  nrmse_SVT_all <- nrmse(data_grnd_full,
                         SVT_imputed)
  
  
  ####  PPCA Imputed  ####
  
  PPCA_imputed<-completeObs(pca(
    object = after_norm_na,
    method = "ppca",
    nPcs = pram,
    center = TRUE,
    scale = "uv"
  ))
  
  PPCA_imputed<-2^PPCA_imputed

  rmse_PPCA_all <- rmse(data_grnd_full,
                        PPCA_imputed)

  nrmse_PPCA_all <- nrmse(data_grnd_full,
                          PPCA_imputed)
  
  
  
  ######## SG/MG Evaluation ########
  before_missing_super_sample <- NULL
  start <- 1
  for(rep in nRep){
    before_missing_super_sample <- cbind(before_missing_super_sample,
                                         rowMeans(data_grnd_full[, start : (start + rep - 1)]))
    start <- start + rep
  }
  
  SGlist<-cotMG(Sest=before_missing_super_sample,cos.thres = 0.9,per = 50)
  index<-as.vector(unlist(SGlist$mg.list))
  index_bool<-row.names(data_grnd_full)%in%index
  
  data_grnd_full_SG<-data_grnd_full[index_bool,]
  
  ABDS_imputed_SG<-ABDS_imputed[index_bool,]
  Half_min_imputed_SG<-Half_min_imputed[index_bool,]
  Mean_imputed_SG<-Mean_imputed[index_bool,]
  swKNN_imputed_SG<-swKNN_imputed[index_bool,]
  NIPALS_imputed_SG<-NIPALS_imputed[index_bool,]
  SVD_imputed_SG<-SVD_imputed[index_bool,]
  SVT_imputed_SG<-SVT_imputed[index_bool,]
  PPCA_imputed_SG<-PPCA_imputed[index_bool,]
  
  ABDS_SG<-c(rmse(data_grnd_full_SG,ABDS_imputed_SG),
             nrmse(data_grnd_full_SG,ABDS_imputed_SG))
  Half_min_SG<-c(rmse(data_grnd_full_SG,Half_min_imputed_SG),
                 nrmse(data_grnd_full_SG,Half_min_imputed_SG))
  Mean_SG<-c(rmse(data_grnd_full_SG,Mean_imputed_SG),
             nrmse(data_grnd_full_SG,Mean_imputed_SG))
  swKNN_SG<-c(rmse(data_grnd_full_SG,swKNN_imputed_SG),
              nrmse(data_grnd_full_SG,swKNN_imputed_SG))
  PPCA_SG<-c(rmse(data_grnd_full_SG,PPCA_imputed_SG),
             nrmse(data_grnd_full_SG,PPCA_imputed_SG))
  NIPALS_SG<-c(rmse(data_grnd_full_SG,NIPALS_imputed_SG),
               nrmse(data_grnd_full_SG,NIPALS_imputed_SG))
  SVD_SG<-c(rmse(data_grnd_full_SG,SVD_imputed_SG),
            nrmse(data_grnd_full_SG,SVD_imputed_SG))
  SVT_SG<-c(rmse(data_grnd_full_SG,SVT_imputed_SG),
            nrmse(data_grnd_full_SG,SVT_imputed_SG))
  
  ################# SG High Group ##################
  data_SG_high<-NULL
  Half_min_imputed_SG_high<-NULL
  Mean_imputed_SG_high<-NULL
  ABDS_imputed_SG_high<-NULL
  swKNN_imputed_SG_high<-NULL
  NIPALS_imputed_SG_high<-NULL
  SVD_imputed_SG_high<-NULL
  SVT_imputed_SG_high<-NULL
  PPCA_imputed_SG_high<-NULL
  
  
  length_SG<-lengths(SGlist$mg.list)
  start_x <- 1
  start_y <- 1
  for(i in 1:length(nRep)){
    data_SG_high <- rbind(data_SG_high,data_grnd_full_SG[start_x:(start_x + length_SG[i] - 1), 
                                                         start_y : (start_y + nRep[i] - 1)])
    ABDS_imputed_SG_high <- rbind(ABDS_imputed_SG_high,ABDS_imputed_SG[start_x:(start_x + length_SG[i] - 1), 
                                                                       start_y : (start_y + nRep[i] - 1)])
    PPCA_imputed_SG_high <- rbind(PPCA_imputed_SG_high,PPCA_imputed_SG[start_x:(start_x + length_SG[i] - 1), 
                                                                       start_y : (start_y + nRep[i] - 1)])
    Half_min_imputed_SG_high <- rbind(Half_min_imputed_SG_high,Half_min_imputed_SG[start_x:(start_x + length_SG[i] - 1), 
                                                                                   start_y : (start_y + nRep[i] - 1)])
    Mean_imputed_SG_high <- rbind(Mean_imputed_SG_high,Mean_imputed_SG[start_x:(start_x + length_SG[i] - 1), 
                                                                       start_y : (start_y + nRep[i] - 1)])
    swKNN_imputed_SG_high <- rbind(swKNN_imputed_SG_high,swKNN_imputed_SG[start_x:(start_x + length_SG[i] - 1), 
                                                                          start_y : (start_y + nRep[i] - 1)])
    NIPALS_imputed_SG_high <- rbind(NIPALS_imputed_SG_high,NIPALS_imputed_SG[start_x:(start_x + length_SG[i] - 1), 
                                                                             start_y : (start_y + nRep[i] - 1)])
    SVD_imputed_SG_high <- rbind(SVD_imputed_SG_high,SVD_imputed_SG[start_x:(start_x + length_SG[i] - 1), 
                                                                    start_y : (start_y + nRep[i] - 1)])
    SVT_imputed_SG_high <- rbind(SVT_imputed_SG_high,SVT_imputed_SG[start_x:(start_x + length_SG[i] - 1), 
                                                                    start_y : (start_y + nRep[i] - 1)])
    start_x <- start_x + length_SG[i]
    start_y <- start_y + nRep[i]
  }
  
  ABDS_SG_high<-c(rmse(data_SG_high,ABDS_imputed_SG_high),
                  nrmse(data_SG_high,ABDS_imputed_SG_high))
  PPCA_SG_high<-c(rmse(data_SG_high,PPCA_imputed_SG_high),
                  nrmse(data_SG_high,PPCA_imputed_SG_high))
  Half_min_SG_high<-c(rmse(data_SG_high,Half_min_imputed_SG_high),
                      nrmse(data_SG_high,Half_min_imputed_SG_high))
  Mean_SG_high<-c(rmse(data_SG_high,Mean_imputed_SG_high),
                  nrmse(data_SG_high,Mean_imputed_SG_high))
  swKNN_SG_high<-c(rmse(data_SG_high,swKNN_imputed_SG_high),
                   nrmse(data_SG_high,swKNN_imputed_SG_high))
  NIPALS_SG_high<-c(rmse(data_SG_high,NIPALS_imputed_SG_high),
                    nrmse(data_SG_high,NIPALS_imputed_SG_high))
  SVD_SG_high<-c(rmse(data_SG_high,SVD_imputed_SG_high),
                 nrmse(data_SG_high,SVD_imputed_SG_high))
  SVT_SG_high<-c(rmse(data_SG_high,SVT_imputed_SG_high),
                 nrmse(data_SG_high,SVT_imputed_SG_high))
  
  ################# Output #################
  ABDS_result <- c(rmse_ABDS_all,nrmse_ABDS_all,ABDS_SG,ABDS_SG_high)
  Mean_result <- c(rmse_Mean_all,nrmse_Mean_all,Mean_SG,Mean_SG_high)
  Half_min_result <- c(rmse_Half_min_all,nrmse_Half_min_all,Half_min_SG,Half_min_SG_high)
  swKNN_result <- c(rmse_swKNN_all,nrmse_swKNN_all,swKNN_SG,swKNN_SG_high)
  PPCA_result <- c(rmse_PPCA_all,nrmse_PPCA_all,PPCA_SG,PPCA_SG_high)
  NIPALS_result <- c(rmse_NIPALS_all,nrmse_NIPALS_all,NIPALS_SG,NIPALS_SG_high)
  SVD_result <- c(rmse_SVD_all,nrmse_SVD_all,SVD_SG,SVD_SG_high)
  SVT_result <- c(rmse_SVT_all,nrmse_SVT_all,SVT_SG,SVT_SG_high)
  
  return(rbind(ABDS_result,Mean_result,Half_min_result,swKNN_result,PPCA_result,NIPALS_result,SVD_result,SVT_result))
}


test1<-test_function(data = data ,rate = 0.2,Total_missing = 30000,pram = 10,sd_scaler=4)


result<-NULL
rate_range<-c(0.2,0.4,0.6,0.8)
Total_missing_range<-c(round(31099*16*0.2),round(31099*16*0.3),round(31099*16*0.4))
for (i in (1:length(rate_range))){
  for (j in (1:length(Total_missing_range))){
    temp<-test_function(data = data,rate=rate_range[i],Total_missing = Total_missing_range[j],pram = 12)
    named_vec<-c(temp[,1],temp[,2],temp[,3],temp[,4],temp[,5],temp[,6])
    result<-cbind(result,named_vec)
  }
}

write.csv(result,'result 19380.csv')
