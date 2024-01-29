library(MASS); library(ggplot2); library(mvtnorm); library(gtools); library(openxlsx); library(writexl); library(magrittr); 
library(dplyr);library(pROC)
source('eCOT.R')
source('scatter_plot.R')

#-----------
#----Basic Parameter Settings----
#-----------

nGroup                  <- 3
psDEG                   <- c(1/3, 1/3, 1/3)    # Absolute percentage over all genes
nGene_general           <- 2400
nGene_sg                <- 60
nGene_dsg               <- 450
overlap_thres           <- 0.98
#----------- 
#----General Genes----
#-----------

alpha                   <- rep(1, nGroup)      # Parameters for Dirichlet distribution
general_genes1          <- rdirichlet(1200, alpha)
alpha                   <- rep(4, nGroup)      # Parameters for Dirichlet distribution
general_genes2          <- rdirichlet(1200, alpha)
general_genes           <- rbind(general_genes1,general_genes2)
scatter_plot(general_genes,mg.info = TRUE, mg = list(1:nGene_general), mg.col = c('blue'))


#-----------
#----Signature Genes----
#-----------

ref <- vector("list", length = nGroup)
for (i in 1:nGroup) {
            ref[[i]] <- rep(0, nGroup)
            ref[[i]][i] <- 1
}



data_sg <- NULL
for (i in 1:nGroup) {
            data_sg <-
                        rbind(data_sg, matrix(rep(ref[[i]], nGene_sg * psDEG[i]), ncol = nGroup, byrow = TRUE))
}



noise <- rnorm(nGroup * data_sg * sum(psDEG), mean = 0, sd = 0.05)
data_sg                 <- abs(data_sg + noise)
num_sg <- nrow(data_sg)
num_sg
scatter_plot(data_sg, mg.info = TRUE, mg = list(1:num_sg), mg.col = c('red'))


#-----------
#----Down-regulated Signature Genes----
#-----------



ref <- vector("list", length = nGroup)
for (i in 1:nGroup) {
  ref[[i]] <- rep(1, nGroup)
  ref[[i]][i] <- 0.1
}

data_dsg <- NULL
for (i in 1:nGroup) {
  data_dsg <-
    rbind(data_dsg, matrix(rep(ref[[i]], nGene_dsg * psDEG[i]), ncol = nGroup, byrow = TRUE))
}


noise <- rnorm(nGroup * data_dsg * sum(psDEG), mean = 0, sd = 0.1)
data_dsg                 <- abs(data_dsg + noise)
num_dsg <- nrow(data_dsg)
num_dsg
scatter_plot(data_dsg, mg.info = TRUE, mg = list(1:num_dsg), mg.col = c('green'))

general_overlapper<-cos_sDEG(general_genes)[,1]
general_genes_filtered<-general_genes[(general_overlapper<overlap_thres),]
scatter_plot(general_genes_filtered)

data_combined<-rbind(data_dsg,data_sg,general_genes_filtered)
scatter_plot(data_combined,mg.info = TRUE, mg = c(list(1:num_dsg),list((num_dsg+1):(num_sg+num_dsg))), mg.col = c('red','green'))


#### Three Methods ####

# OVR-FC by my understanding

ovr_fc_geomean<-function(input){
  a_bc<-input[,1]/exp(mean(log(input[,2:3])))      # R U Sure this is the geo mean?
  b_ac<-input[,2]/exp(mean(log(input[,c(1,3)])))
  c_ab<-input[,3]/exp(mean(log(input[,1:2])))
  temp<-cbind(a_bc,b_ac,c_ab)
  min_fc<-apply(temp,1,min)
  return(min_fc)
}


result_ove_fc<-ovr_fc_geomean(data_combined)

gt<-c(rep(0,nGene_dsg),rep(1,(dim(data_combined)[1]-nGene_dsg)))
predict_fc<-result_ove_fc
predict_fc<-predict_fc/max(predict_fc)
roc_curve_fc <- roc(gt, predict_fc)

scatter_plot(data_combined,mg.info = TRUE, mg = list(which(predict_fc<sort(predict_fc)[nGene_dsg])), mg.col = c('red'))


# eCOT
predict_ecot<-1-as.numeric(eCOT(data_combined)[,1])
predict_ecot<-predict_ecot/max(predict_ecot)
gt<-c(rep(0,nGene_dsg),rep(1,(dim(data_combined)[1]-nGene_dsg)))
roc_curve_ecot <- roc(gt, predict_ecot)

scatter_plot(data_combined,mg.info = TRUE, mg = list(which(predict_ecot<sort(predict_ecot)[nGene_dsg])), mg.col = c('red'))


# T-Test

ttest_dsg<-function(input){
  result<-NULL
  for (i in 1:dim(input)[1]){
    a_bc<-t.test(x=log(input[i,2:3]), mu=log(input[i,1]), var.equal = TRUE,alternative = "greater")$p.value
    b_ac<-t.test(x=log(input[i,c(1,3)]), mu=log(input[i,2]), var.equal = TRUE,alternative = "greater")$p.value
    c_ab<-t.test(x=log(input[i,1:2]), mu=log(input[i,3]), var.equal = TRUE,alternative = "greater")$p.value
    temp<-c(a_bc,b_ac,c_ab)
    result<-rbind(result,c(min(temp),which.min(temp)))
  }
  return(result)
}

result_ttest<-ttest_dsg(data_combined)

predict_ttest<-result_ttest[,1]
predict_ttest<-predict_ttest/max(predict_ttest)
gt<-c(rep(0,nGene_dsg),rep(1,(dim(data_combined)[1]-nGene_dsg)))
roc_curve_ttest <- roc(gt, predict_ttest)

scatter_plot(data_combined,mg.info = TRUE, mg = list(which(predict_ttest<sort(predict_ttest)[nGene_dsg])), mg.col = c('red'))

#### ROC & AUC ####

plot((1-roc_curve_fc$specificities), roc_curve_fc$sensitivities, 
     type = "l", xlab = "False Positive Rate", ylab = "True Positive Rate", col = "blue",lwd = 1.5)
lines(1-roc_curve_ecot$specificities, roc_curve_ecot$sensitivities, 
      type = "l", col = "red",lwd = 1.5)
lines(1-roc_curve_ttest$specificities, roc_curve_ttest$sensitivities, 
      type = "l", col = "green",lwd = 1.5)
legend("bottomright", legend=c("eCOT", "T-Test", "OVR-FC"),col=c("red","green", "blue"), lty=1:3, cex=0.8)

# pROC
thres<-0.05

ind_fc<-(1-roc_curve_fc$specificities)<thres
ind_ecot<-(1-roc_curve_ecot$specificities)<thres
ind_ttest<-(1-roc_curve_ttest$specificities)<thres

plot((1-roc_curve_fc$specificities)[ind_fc], roc_curve_fc$sensitivities[ind_fc], 
     type = "l", xlab = "False Positive Rate", ylab = "True Positive Rate",xlim = c(0,thres),ylim = c(0,1), col = "blue",lwd = 1.5)
legend("right", legend=c("eCOT", "T-Test", "OVR-FC"),col=c("red","green", "blue"), lty=1:3, cex=0.8)
lines(1-roc_curve_ecot$specificities[ind_ecot], roc_curve_ecot$sensitivities[ind_ecot], 
     type = "l",xlim = c(0,thres),ylim = c(0,1), col = "red",lwd = 1.5)
lines(1-roc_curve_ttest$specificities[ind_ttest], roc_curve_ttest$sensitivities[ind_ttest], 
     type = "l",xlim = c(0,thres),ylim = c(0,1), col = "green",lwd = 1.5)



auc(roc_curve_fc)
auc(roc_curve_ecot)
auc(roc_curve_ttest)

# pAUC
auc(roc_curve_fc, partial.auc=c(1,1-thres), partial.auc.focus="sp")/thres
auc(roc_curve_ecot, partial.auc=c(1,1-thres), partial.auc.focus="sp")/thres
auc(roc_curve_ttest, partial.auc=c(1,1-thres), partial.auc.focus="sp")/thres
