#' @name MGpI
#' @title Mechanism-integrated group-wise pre-imputation
#' @description
#' A mechanism-integrated group-wise pre-imputation (MGpI) strategy that explicitly considers 
#' mixed missing mechanisms across different phenotypic groups
#' @param readin Read-in data, with features on the row and samples (grouped) on the columns.
#' Samples of the same group should be placed together.
#' @param nRep The number of samples in each group, e.g., c(5,4,3,2), corresponding to the readin data.
#' @param min_option The options for the minimum value, with 'global' (default) meaning to use the global minimum value across 
#' the whole data matrix, or 'local' meaning to use the feature-wise minimum value.
#' @param zero.ignore Whether to ignore the zero values in the input data. If set to TRUE, zero values in the input data will
#' not be seen as missing values thus will not be imputed. Default is FALSE (to see zero-values as missing and impute them).
#' @param missing_rate_threshold A value between 0 and 1. The upper threshold of sample-wise missing-rate allowed for features 
#' in any of each group, with 1 as default (all features are allowed).
#' @param sd_scaler For scaling up the standard deviation of the distribution obtained from the non-missing value. Default
#' is 1 (no scaling).
#' 
#' @return the normalized data.
#' 
#' @examples 
#' MGpIed <- MGpI(readin=data, nRep=c(5,4,3,2))
#' @export

MGpI<-function(readin,nRep,min_option='global',zero.ignore=FALSE,missing_rate_threshold=1, sd_scaler=1){
  
  before_remove<-readin
  rownames<-row.names(readin)
  
  message("MGpI: Pre-processing")
  
  # Pre-processing & Initialization
  
  for (i in 1:dim(before_remove)[2]){   
    temp<-as.numeric(before_remove[,i])
    before_remove[,i]<-temp
  }
  
  if (zero.ignore==TRUE){
    zero_position<-(before_remove==0) 
    # Store the position of zero values. If zero.ignore is set to TRUE, these zero values
    # will be imputed, but set back to 0 in the end
  }
  
  before_remove[is.na(before_remove)]<-0   
  after_remove<-NULL 
  index_name<-NULL
  
  
  # Select features with missing rate less than (a given threshold) in any of each group
  
  before_impute<-before_remove
  
  nRep_sum<-vector() ## suppose nRep is 21 39 34 2, nRep_sum is 21 60 94 96
  for (i in 1:length(nRep)){
    nRep_sum[i]<-sum(nRep[1:i])
  }
  nRep_plus<-c(0,nRep_sum)  # nRep_plus will be 0 21 60 94 96
  non_0<-list() # this stores the number of non-0 values in each group
  
  
  if (missing_rate_threshold!=1){
    message("MGpI: Selecting features with missing rate less than the threshold you choose in any of each group")
    
    for (j in 1:dim(before_remove)[1]){
      for (i in 1:length(nRep)){
        non_0[i]<-sum(before_remove[j,(nRep_plus[i]+1):(nRep_plus[i+1])] !=0)
        # the number of non-0 values in each group (iteratively, all features)
      }
      
      if (any(as.matrix(non_0)>as.matrix(nRep)*(1-missing_rate_threshold))){ 
        after_remove<-rbind(after_remove,before_remove[j,])
        index_name<-c(index_name,j)
        # remove unqualified features
      }
    }
    
    before_impute<-after_remove
    
    row.names(before_impute)<-rownames[index_name]
    
  }
  
  
  # Imputation
  
  ## Initializing
  
  after_impute<-NULL 
  # what's after imputation
  
  protein_segment<-list() 
  # Consider a protein is recorded in 10 samples where 5 samples are in one group and other 5 in the other one.
  # Data in the two groups are separately moved into the list as two elements, where each one is a 5-value vector.
  # They are treated as two 'protein_segment' here, independently for imputation.
  
  log_segment<-list()
  # log-transformed protein_segment
  
  overall_LLOD<-0
  # count of the LLOD value(s)
  
  ## 'local' refers to using feature-wise minimum value for imputation. 
  ## The other choice is 'global' that uses global minimum value.
  
  if (min_option == 'local'){
    message("MGpI: Using weighted feature-wise minimum value & group-mean for feature-wise pre-imputation")
    
    for (j in 1:dim(before_impute)[1]){
      
      this_protein<-before_impute[j,]
      
      if (length(this_protein)==length(this_protein[this_protein==0])){
        message("MGpI: Features with all value missing detected! Skip to the next feature. Consider removing.")
        after_impute<-rbind(after_impute,this_protein)
        overall_LLOD<-overall_LLOD+length(this_protein)
        
      } else {
        
        # log-transformed non-zero feature-wise minimum value
        feature_min<-min(this_protein[this_protein!=0])
        feature_min_log<-log(feature_min) 
        
        for (i in 1:length(nRep)){
          non_0[i]<-sum(before_impute[j,(nRep_plus[i]+1):(nRep_plus[i+1])] !=0) 
          # loading the number of non-zero values into temporary storage
          protein_segment[[i]]<-before_impute[j,(nRep_plus[i]+1):(nRep_plus[i+1])]
          # imputing the current 'protein_segment'

          if (non_0[i]==0){
            protein_segment[[i]][protein_segment[[i]]==0]<-feature_min/2
            # if the current 'protein_segment' contains only zero values, impute them with min/2
            missing_LLOD_rate<-1
            
          } else if (non_0[i]==1){
            protein_segment[[i]][protein_segment[[i]]==0]<-feature_min/2
            # if the current 'protein_segment' contains only one non-zero value, impute the rest of them with min/2
            missing_LLOD_rate<-(length(protein_segment[[i]])-1)/length(protein_segment[[i]])
            
          } else {
            log_segment[[i]]<-log(protein_segment[[i]][protein_segment[[i]]!=0])
            # if the current 'protein_segment' contains more than one non-zero values
            missing_LLOD_rate<-pnorm(feature_min_log,mean = mean(log_segment[[i]]), sd = sd_scaler*sd(log_segment[[i]]))
            # estimate the proportion of missingness of the current 'protein_segment' ...
            protein_segment[[i]][protein_segment[[i]]==0]<-
              missing_LLOD_rate * feature_min/2 + 
              (1 - missing_LLOD_rate) * mean(protein_segment[[i]][protein_segment[[i]]!=0])
          } # ... and impute them accordingly
          overall_LLOD<-overall_LLOD + missing_LLOD_rate*length(protein_segment[[i]])
        }
      
        after_impute<-rbind(after_impute,unlist(protein_segment))
      }
    }
  }
  
  
  
  if (min_option == 'global'){
    message("MGpI: Using weighted global minimum value & feature-wise group-mean for feature-wise pre-imputation")
    # same process, using the global minimum instead of feature minimum
    for (j in 1:dim(before_impute)[1]){
      
      this_protein<-before_impute[j,]
      
      if (length(this_protein)==length(this_protein[this_protein==0])){
        message("MGpI: Features with all value missing detected! Skip to the next feature. Consider removing.")
        after_impute<-rbind(after_impute,this_protein)
        overall_LLOD<-overall_LLOD+length(this_protein)
        
      } else {
        
        global_min<-min(before_impute[before_impute!=0])
        global_min_log<-log(global_min)
        
        for (i in 1:length(nRep)){
          non_0[i]<-sum(before_impute[j,(nRep_plus[i]+1):(nRep_plus[i+1])] !=0)
          protein_segment[[i]]<-before_impute[j,(nRep_plus[i]+1):(nRep_plus[i+1])]
        }
        
        for (i in 1:length(nRep)){
          if (non_0[i]==0){
            protein_segment[[i]][protein_segment[[i]]==0]<-global_min/2
            missing_LLOD_rate<-1
            
          } else if (non_0[i]==1){
            protein_segment[[i]][protein_segment[[i]]==0]<-global_min/2
            missing_LLOD_rate<-(length(protein_segment[[i]])-1)/length(protein_segment[[i]])
            
          } else {
            log_segment[[i]]<-log(protein_segment[[i]][protein_segment[[i]]!=0])
            missing_LLOD_rate<-pnorm(global_min_log,mean = mean(log_segment[[i]]), sd = sd_scaler*sd(log_segment[[i]]))
            protein_segment[[i]][protein_segment[[i]]==0]<-
              missing_LLOD_rate * global_min/2 + 
              (1 - missing_LLOD_rate) * mean(protein_segment[[i]][protein_segment[[i]]!=0])
          }
          overall_LLOD<-overall_LLOD + missing_LLOD_rate*length(protein_segment[[i]])
        }
        
        after_impute<-rbind(after_impute,unlist(protein_segment))
      }
    }
  }
  

  row.names(after_impute)<-rownames[index_name]

  message("MGpI: Done.")
  
  if (zero.ignore==TRUE){
    after_impute[zero_position]<-0
  }
  
  return(after_impute)
  
}