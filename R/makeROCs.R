#' Multi-ROC Plotting for TCGA Gene Expression Data
#'
#' This function generates a multi-Receiver Operating Characteristic (ROC) plot using RNA-seq data from The Cancer Genome Atlas (TCGA) database and a user-specified target variable.  TCGA data is imported using the TCGA2STAT package.  Feature selection is performed, and the remaining variables (in this case, genes) are processed by five different machine learning classifiers: LASSO-Logistic, K-Nearest Neighbor, Random Forest, a Radial-Kernal Support Vector Machine, and a Sigmoid-Kernal Support Vector Machine.  An ROC curve is generated from each classifier, and plotted onto one graph.  
#'
#' @param type TCGA-supported acronyms that designate the type of cancer.  makeROCs currently supports ACC, BLCA, KIRC, KIRP, LIHC, and THCA. Instead of specifying a type, “Random” can be used to generate a randomly-chosen type of cancer, as well as a randomly chosen target variable
#'
#' @param target the variable to be predicted.  makeROCs currently supports prediction of tumor “stage” (which attempts to distinguish stage I tumors from stage II, III, and IV tumors) and the participant’s “gender”.
#'
#' @details Cancer type acronyms:
#' ACC	  Adrenocortical Carcinoma,
#' BLCA	  Bladder Urothelial Carcinoma,
#' KIRC	  Kidney Renal Clear Cell Carcinoma,
#' KIRP	  Kidney Renal Papillary Cell Carcinoma,
#' LIHC	  Liver Heptocellular Carcinoma,
#' THCA	  Thyroid Carcinoma
#'
#' @author Jacob Blamer, \email{jwilliamblamer@gmail.com}
#'
#' @examples
#' makeROCs("THCA", "Stage")
#' makeROCs("BLCA", "Gender")
#' makeROCs("Random")
#'
#' @import 
#' caret
#' randomForest
#' e1071
#' glmnet
#' pROC
#' TCGA2STAT
#' MXM
#' ROCR
#'
#'
#'
#'
#'
#'
#' @export
makeROCs <- function(type, target) {

  stageROCs <- function(type) {
    df1 <- getTCGA(type, "RNASeq2", "RPKM", clinical = T, cvars="pathologicstage")
    
    ## Data Prep
    df1.md <- df1$merged.dat
    stage_i <- which(df1.md[,2] == "stage i")
    stage_i0 <- which(df1.md[,2] != "stage i")
    df1.md[stage_i,2] <- "si"
    df1.md[stage_i0,2] <- "si0"
    df1.md <- df1.md[-which(is.na(df1.md[,2])),]
    
    ## Feature Selection
    print("[TIMESTAMP]: Feature Selection Started")
    print(Sys.time())
    mxm <- MXM::MMPC(target  = as.factor(df1.md[,2]),            
                     dataset = df1.md[,3:20503],            
                     max_k = floor(nrow(df1.md)/10),          
                     threshold = 0.05,                                         
                     test = NULL)
    mxm_idxs <- mxm@selectedVars
    print("[TIMESTAMP]: Feature Selection Finished")
    print(Sys.time())
    
    ## Train/Test
    idx <- sample(nrow(df1.md), nrow(df1.md)*.75)
    f1_trainset <- df1.md[idx,]
    f1_testset <- df1.md[-idx,]
    
    ########################
    #### LASSO Logistic ####
    ########################
    print("[TIMESTAMP]: LASSO Logistic Started")
    print(Sys.time())
    
    f1_lasso <- cv.glmnet(data.matrix(f1_trainset[,mxm_idxs]), f1_trainset[,2], nfolds = 10, family="binomial", alpha=1, type.measure = "auc")
    
    f1_lasso.prd <- predict(f1_lasso, newx = data.matrix(f1_testset[,mxm_idxs]),
                            type="response", s=f1_lasso$lambda.min)[,1]
    
    lassoROC <- roc(f1_testset[,2], f1_lasso.prd)
    
    print("[TIMESTAMP]: LASSO Logistic Finished")
    print(Sys.time())
    
    #############
    #### KNN ####
    #############
    print("[TIMESTAMP]: KNN Started")
    print(Sys.time())
    
    ctrl <- trainControl(method="cv", number=10)
    f1_knn <- train(PATHOLOGICSTAGE~., 
                    data = f1_trainset[,c(2,mxm_idxs)], 
                    method = "knn", 
                    trControl = ctrl, 
                    tuneLength = 20)
    
    knn_test.prd <- predict(f1_knn, newdata = f1_testset, type = "prob")
    knnROC <- roc(f1_testset[,2], knn_test.prd[,"si"])
    
    print("[TIMESTAMP]: KNN Finished")
    print(Sys.time())
    
    ############
    #### RF ####
    ############
    print("[TIMESTAMP]: RF Started")
    print(Sys.time())
    
    trcontrol <- trainControl(method='repeatedcv', 
                              number=10,            
                              repeats=3,            
                              search='grid')
    sqrt_p <- sqrt(length(mxm_idxs))
    tunegrid <- expand.grid(.mtry = floor((sqrt_p-(sqrt_p/2)):(sqrt_p+(sqrt_p/2))))         
    rf_grid <- train(PATHOLOGICSTAGE~., 
                     data = f1_trainset[,c(2,mxm_idxs)], 
                     method = "rf",
                     tuneGrid = tunegrid)
    rf <- randomForest(as.factor(PATHOLOGICSTAGE)~., 
                       data = f1_trainset[,c(2,mxm_idxs)], 
                       proximity=TRUE, 
                       mtry=rf_grid$bestTune[[1]])
    
    rf_prds <- predict(rf, f1_testset, type="prob")
    rfROC <- roc(f1_testset[,2], rf_prds[,"si"])
    
    print("[TIMESTAMP]: RF Finished")
    print(Sys.time())
    
    #####################
    #### SVM (Radial) ###
    #####################
    print("[TIMESTAMP]: SVM Radial Started")
    print(Sys.time())
    
    tune.out <- tune(svm, as.factor(PATHOLOGICSTAGE)~., 
                     data = f1_trainset[,c(2,mxm_idxs)],
                     kernel ="radial", probability=TRUE,
                     ranges = list(cost=c(0.001 , 0.01, 0.1, 1,5,10,100)))
    
    f1_svm.radial <- tune.out$best.model
    temp_radial_prds <- predict(as.vector(f1_svm.radial), f1_testset, probability = T)
    
    svm.radial_prds <- attributes(temp_radial_prds)$probabilities
    radROC <- roc(f1_testset[,2], svm.radial_prds[,"si"])
    
    print("[TIMESTAMP]: SVM Radial Finished")
    print(Sys.time())
    
    #######################
    #### SVM (Sigmoid) ####
    #######################
    print("[TIMESTAMP]: SVM Sigmoid Started")
    print(Sys.time())
    
    tune.out <- tune(svm, as.factor(PATHOLOGICSTAGE)~., 
                     data = f1_trainset[,c(2,mxm_idxs)],
                     kernel ="sigmoid", probability=TRUE,
                     ranges = list(cost=c(0.001 , 0.01, 0.1, 1,5,10,100)))
    
    f1_svm.sig <- tune.out$best.model
    temp_sig_prds <- predict(as.vector(f1_svm.sig), f1_testset, probability = T)
    
    svm.sig_prds <- attributes(temp_sig_prds)$probabilities
    sigROC <- roc(f1_testset[,2], svm.sig_prds[,"si"])
    
    ## Plot ROCs
    plot(lassoROC, col = 2, main=paste("ROC Curves: Predicting Tumor Stage in ", type))
    plot(knnROC, col = 4, add=TRUE)
    plot(rfROC, col = 6, add=TRUE)
    plot(radROC, col = 11, add=TRUE)
    plot(sigROC, col = 13, add=TRUE)
    legend(x = .25, y = .7, cex = 1, border = "white", bty="n", y.intersp = .55,
           legend = c(paste("Lasso Logistic (AUC: ", round(lassoROC$auc,2), ")"), 
                      paste("KNN (AUC: ", round(knnROC$auc,2), ")"), 
                      paste("Random Forest (AUC: ", round(rfROC$auc,2), ")"), 
                      paste("SVM - Radial (AUC: ", round(radROC$auc,2), ")"), 
                      paste("SVM - Sigmoid (AUC: ", round(sigROC$auc,2), ")")),
           fill = c(2,4,6,11,13))
    
  }
  
  
  genderROCs <- function(type) {
    df1 <- getTCGA(type, "RNASeq2", "RPKM", clinical = T, cvars="gender")
    
    ## Data Prep
    df1.md <- df1$merged.dat
    
    ## Feature Selection
    print("[TIMESTAMP]: Feature Selection Started")
    Sys.time()
    mxm <- MXM::MMPC(target  = as.factor(df1.md[,2]),            
                     dataset = df1.md[,3:20503],            
                     max_k = floor(nrow(df1.md)/10),          
                     threshold = 0.05,                                         
                     test = NULL)
    mxm_idxs <- mxm@selectedVars
    print("[TIMESTAMP]: Feature Selection Finished")
    Sys.time()
    
    ## Train/Test
    idx <- sample(nrow(df1.md), nrow(df1.md)*.75)
    f1_trainset <- df1.md[idx,]
    f1_testset <- df1.md[-idx,]
    
    ########################
    #### LASSO Logistic ####
    ########################
    print("[TIMESTAMP]: LASSO Logistic Started")
    Sys.time()
    
    f1_lasso <- cv.glmnet(data.matrix(f1_trainset[,mxm_idxs]), f1_trainset[,2], nfolds = 10, family="binomial", alpha=1, type.measure = "auc")
    
    f1_lasso.prd <- predict(f1_lasso, newx = data.matrix(f1_testset[,mxm_idxs]),
                            type="response", s=f1_lasso$lambda.min)[,1]
    
    lassoROC <- roc(f1_testset[,2], f1_lasso.prd)
    
    print("[TIMESTAMP]: LASSO Logistic Finished")
    Sys.time()
    
    ############
    #### KNN ####
    #############
    print("[TIMESTAMP]: KNN Started")
    Sys.time()
    
    ctrl <- trainControl(method="cv", number=10)
    f1_knn <- train(GENDER~., 
                    data = f1_trainset[,c(2,mxm_idxs)], 
                    method = "knn", 
                    trControl = ctrl, 
                    tuneLength = 20)
    
    knn_test.prd <- predict(f1_knn, newdata = f1_testset, type = "prob")
    knnROC <- roc(f1_testset[,2], knn_test.prd[,"male"])
    
    print("[TIMESTAMP]: KNN Finished")
    Sys.time()
    
    ############
    #### RF ####
    ############
    print("[TIMESTAMP]: RF Started")
    Sys.time()
    
    trcontrol <- trainControl(method='repeatedcv', 
                              number=10,            
                              repeats=3,            
                              search='grid')
    sqrt_p <- sqrt(length(mxm_idxs))
    tunegrid <- expand.grid(.mtry = floor((sqrt_p-(sqrt_p/2)):(sqrt_p+(sqrt_p/2))))         
    rf_grid <- train(GENDER~., 
                     data = f1_trainset[,c(2,mxm_idxs)], 
                     method = "rf",
                     tuneGrid = tunegrid)
    rf <- randomForest(as.factor(GENDER)~., 
                       data = f1_trainset[,c(2,mxm_idxs)], 
                       proximity=TRUE, 
                       mtry=rf_grid$bestTune[[1]])
    
    rf_prds <- predict(rf, f1_testset, type="prob")
    rfROC <- roc(f1_testset[,2], rf_prds[,"male"])
    
    print("[TIMESTAMP]: RF Finished")
    Sys.time()
    
    #####################
    #### SVM (Radial) ###
    #####################
    print("[TIMESTAMP]: SVM Radial Started")
    Sys.time()
    
    tune.out <- tune(svm, as.factor(GENDER)~., 
                     data = f1_trainset[,c(2,mxm_idxs)],
                     kernel ="radial", probability=TRUE,
                     ranges = list(cost=c(0.001 , 0.01, 0.1, 1,5,10,100)))
    
    f1_svm.radial <- tune.out$best.model
    temp_radial_prds <- predict(as.vector(f1_svm.radial), f1_testset, probability = T)
    
    svm.radial_prds <- attributes(temp_radial_prds)$probabilities
    radROC <- roc(f1_testset[,2], svm.radial_prds[,"male"])
    
    print("[TIMESTAMP]: SVM Radial Finished")
    Sys.time()
    
    #######################
    #### SVM (Sigmoid) ####
    #######################
    print("[TIMESTAMP]: SVM Sigmoid Started")
    Sys.time()
    
    tune.out <- tune(svm, as.factor(GENDER)~., 
                     data = f1_trainset[,c(2,mxm_idxs)],
                     kernel ="sigmoid", probability=TRUE,
                     ranges = list(cost=c(0.001 , 0.01, 0.1, 1,5,10,100)))
    
    f1_svm.sig <- tune.out$best.model
    temp_sig_prds <- predict(as.vector(f1_svm.sig), f1_testset, probability = T)
    
    svm.sig_prds <- attributes(temp_sig_prds)$probabilities
    sigROC <- roc(f1_testset[,2], svm.sig_prds[,"male"])
    
    ## Plot ROCs
    plot(lassoROC, col=2, main=paste("ROC Curves: Predicting Gender in ", type))
    plot(knnROC, col = 4, add=TRUE)
    plot(rfROC, col = 6, add=TRUE)
    plot(radROC, col = 11, add=TRUE)
    plot(sigROC, col = 13, add=TRUE)
    legend(x = .25, y = .7, cex = 1, border = "white", bty="n", y.intersp = .55,
           legend = c(paste("Lasso Logistic ( AUC: ", round(lassoROC$auc,4), ")"), 
                      paste("KNN ( AUC: ", round(knnROC$auc,4), ")"), 
                      paste("Random Forest ( AUC: ", round(rfROC$auc,4), ")"), 
                      paste("SVM - Radial ( AUC: ", round(radROC$auc,4), ")"), 
                      paste("SVM - Sigmoid ( AUC: ", round(sigROC$auc,4), ")")),
           fill = c(2,4,6,11,13))
  }



  
  if (type == "Random") {
    gen_type <- c("BLCA", "KIRC", "KIRP", "LIHC", "ACC", "THCA")
    gen_target <- c("Stage", "Gender")
    
    type <- sample(gen_type, 1)
    target <- sample(gen_target, 1)
    paste("Randomly Generated Cancer Type: ", type)
    print("Randomly Generated Target Variable: ", target)
    
    if (gen_target == "Stage") {
      stageROCs(type)
    }
    
    else if (gen_target == "Gender") {
      genderROCs(type)
    }
    
  }
  
  else if (type != "Random") {
    if (target == "Stage") {
      stageROCs(type)
    }
    
    else if (target == "Gender") {
      genderROCs(type)
    }
  }
  
  
  else {
    print("ERROR! Invalid Imput!")
  }
}


