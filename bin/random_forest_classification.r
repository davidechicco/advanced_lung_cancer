setwd(".")
options(stringsAsFactors = FALSE)
cat("\014")
set.seed(12)

EXP_ARG_NUM <- 2

TRAIN_SET_OVERSAMPLING_SYNTHETIC <- TRUE

threshold <- 0.5

fileNameData <- "../data/pone0210951_s006_dataset_EDITED.csv"
targetName <- "Sepsis"


MISSING_DATA_IMPUTATION <- TRUE

list.of.packages <- c("pacman")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
library("pacman")
pacman::p_load("randomForest", "ggplot2", "dplyr", "pastecs",  "mice", "ROSE")

source("./confusion_matrix_rates.r")
source("./utils.r")

NUM_METRICS <- 9
confMatDataFrame <- matrix(ncol=NUM_METRICS, nrow=1)
colnames(confMatDataFrame) <- c("MCC", "F1 score", "accuracy", "TP_rate", "TN_rate", "PPV", "NPV", "PR_AUC", "ROC_AUC")


FEATURE_RANKING_PLOT_DEPICTION <-  FALSE
TWO_FEATURES_PLOT <- FALSE

patients_data <- read.csv(fileNameData, header = TRUE, sep =",");
patients_data_original <- patients_data
cat("Read data from file ", fileNameData, "\n", sep="")

patients_data$"Histology_coding" <- NULL
patients_data$"Cause_Direct" <- NULL	
patients_data$"Cause_analysis" <- NULL
patients_data$"MICU_Cause_explain" <- NULL
patients_data$"culture_result" <- NULL
patients_data$"Deathcauseexp" <- NULL
patients_data$"ADM_cause_detailed" <- NULL
patients_data$"treatment_detailed" <- NULL
patients_data$"ICU_ADM_Y" <- NULL
patients_data$"RRTetc" <-NULL
patients_data$"serial_no" <-NULL
patients_data$"immunotherapy" <-NULL
patients_data$"admission_via" <-NULL
patients_data$"Underlying_etc" <-NULL




patients_data$"StageNum" <- -1
patients_data[which(patients_data$"Stage" == "ED"),]$"StageNum" <- 5
patients_data[which(patients_data$"Stage" == "IIIB"),]$"StageNum" <- 3
patients_data[which(patients_data$"Stage" == "IV"),]$"StageNum" <- 4
patients_data$"Stage" <- NULL

patients_data$"histology" %>% unique() %>% sort() %>% cat("\n", sep="\n")
patients_data$"histologyNum" <- -1
patients_data[which(patients_data$"histology" == "adenoca"),]$"histologyNum" <- 0
patients_data[which(patients_data$"histology" == "large_cell_NEC"),]$"histologyNum" <- 1
patients_data[which(patients_data$"histology" == "mixed"),]$"histologyNum" <- 2
patients_data[which(patients_data$"histology" == "neuroendocrine_tumor"),]$"histologyNum" <- 3
patients_data[which(patients_data$"histology" == "nonsmall_cell"),]$"histologyNum" <- 4
patients_data[which(patients_data$"histology" == "others"),]$"histologyNum" <- 5
patients_data[which(patients_data$"histology" == "poorly_differentiated"),]$"histologyNum" <- 6
patients_data[which(patients_data$"histology" == "sarcomatoid"),]$"histologyNum" <- 7
patients_data[which(patients_data$"histology" == "small_cell"),]$"histologyNum" <- 8
patients_data[which(patients_data$"histology" == "squamous"),]$"histologyNum" <- 9
patients_data$"histology" <- NULL

patients_data <- patients_data %>% select(order(colnames(patients_data)))


names(patients_data)[names(patients_data) == targetName] <- "target"
cat("The target feature is ", targetName, "\n", sep="")

cat("application of dplyr::select()\n")
patients_data <- patients_data%>%dplyr::select(-target,target)
target_index <- ncol(patients_data)


num_to_return <- 1
upper_num_limit <- 100000
exe_num <- sample(1:upper_num_limit, num_to_return)


# patients_data$"HCV.RNATaqman.Log.IU.ml." <- as.numeric(patients_data$"HCV.RNATaqman.Log.IU.ml.")

# patients_data$"category_0healthy_1hepatitis_2fibrosis_3cirrhorsis" <- NULL
if(MISSING_DATA_IMPUTATION==TRUE){

    # missing data imputation
    NUM_DATASETS <- 1
    imputed_data <- mice(patients_data, m=NUM_DATASETS, maxit = 50, method = 'pmm', seed = 500)
    patients_data <- complete(imputed_data, NUM_DATASETS)
}

patients_data[,target_index] <- as.factor(patients_data[,target_index])

# formula
allFeaturesFormula <- as.formula(paste(as.factor(colnames(patients_data)[target_index]), '.', sep=' ~ ' ))

# cycle of executions

execution_number <- 100
cat("Number of executions = ", execution_number, "\n", sep="")
for(exe_i in 1:execution_number)
{
    cat("[Execlution number ", exe_i, " out of ", execution_number, "]\n", sep="" )
    cat("[Randomizing the rows]\n")
    patients_data <- patients_data[sample(nrow(patients_data)),] # shuffle the rows

    totalElements <- dim(patients_data)[1]

    subsets_size <- totalElements

    target_label <- colnames(patients_data[target_index])
    cat("target_label = ", target_label, "\n", sep="")

    if (subsets_size != totalElements) {
        cat("ATTENTION: We are running the method on a subset of the original dataset, \n", sep="")
        cat(" containing only ", subsets_size, " elements \n", sep="")
        cat(" instead of ", totalElements, " elements \n", sep="")
    }

    patients_data <- patients_data[1:subsets_size, ]

    dataset_dim_retriever(patients_data)
    imbalance_retriever(patients_data[,target_index])

    training_set_perc <- 80
    INPUT_PERC_POS <- 50
    cat("[training set = ", training_set_perc,"%]\n", sep="")
    cat("[test set = ", (100-training_set_perc),"%]\n", sep="")

    artificialBalance <- FALSE
    balancedFlag <- FALSE # flag that sets everything to 50% 50% ratio

    if (artificialBalance == TRUE) {


        train_data_balancer_output <- train_data_balancer(patients_data, target_index, training_set_perc, INPUT_PERC_POS, balancedFlag)

        patients_data_train <- train_data_balancer_output[[1]]
        patients_data_test <- train_data_balancer_output[[2]]
        
        # Creating the subsets for the targets
        patients_data_train_labels <- patients_data_train[, target_index] # NEW
        patients_data_test_labels <- patients_data_test[, target_index]   # NEW

    } else {


        # the training set is the first 60% of the whole dataset
        training_set_first_index <- 1 # NEW
        training_set_last_index <- round(dim(patients_data)[1]*training_set_perc/100) # NEW

        # the test set is the last 40% of the whole dataset
        test_set_first_index <- training_set_last_index+1 # NEW
        test_set_last_index <- dim(patients_data)[1] # NEW

        cat("[Creating the training set and test set for the values]\n")
        patients_data_train <- patients_data[training_set_first_index:training_set_last_index, 1:(target_index)] # NEW
        patients_data_test <- patients_data[test_set_first_index:test_set_last_index, 1:(target_index)] # NEW
        
        # train_balanced_over <- ovun.sample(allFeaturesFormula, data = patients_data_train, method = "both",  p=0.5, N = (nrow(patients_data_train)*2))$data
        # patients_data_train <- train_balanced_over
         
         # https://www.analyticsvidhya.com/blog/2016/03/practical-guide-deal-imbalanced-classification-problems/
         
         cat("\ncheck\n")
         
         if(TRAIN_SET_OVERSAMPLING_SYNTHETIC == TRUE)
         {
            thisP <- 0.5
         
            data.rose <- ROSE(allFeaturesFormula, data = patients_data_train, p=thisP, seed = 1)$data
            patients_data_train <- data.rose
         }
        
        cat("[training set dimensions: ", dim(patients_data_train)[1], " patients]\n")

        cat("[test set dimensions: ", dim(patients_data_test)[1], " patients]\n")

        cat("[Creating the training set and test set for the labels \"1\"-\"0\"]\n")
        patients_data_train_labels <- patients_data_train[, target_index] # NEW
        patients_data_test_labels <- patients_data[test_set_first_index:test_set_last_index, target_index]   # NEW

    }


    dataset_dim_retriever(patients_data_train)
    imbalance_retriever(patients_data_train[, "target"])    

    cat("\n[Training the random forest classifier on the training set]\n")

    rf_new <- NULL
    rf_new <- randomForest(allFeaturesFormula, data=patients_data_train, importance=FALSE, proximity=TRUE)
    
    cat("\n[Applying the trained random forest classifier on the test set]\n")
    patients_data_test_PRED <- predict(rf_new, patients_data_test, type="prob")[,"1"]

    thisConfMat <- confusion_matrix_rates(patients_data_test_labels, patients_data_test_PRED, "@@@ Test set @@@")
    
    if (exe_i == 1)  confMatDataFrame <-  thisConfMat
    else  confMatDataFrame <- rbind(confMatDataFrame, thisConfMat)
    
 }
 
cat("\n\n\n=== final results ===\n")
cat("Number of executions = ", execution_number, "\n", sep="")

# statistics on the dataframe of confusion matrices
statDescConfMatr <- stat.desc(confMatDataFrame)
meanSigmaRowResults <- (statDescConfMatr)[c("mean","std.dev"),]
print(dec_three(statDescConfMatr))
cat("\n\n=== === === ===\n")
print(dec_three(meanSigmaRowResults))
cat("\n\n=== === === ===\n\n\n")

print(sort(colnames(patients_data)))

# printResultsLatex("Random forests", meanSigmaRowResults)

#   cat("\t", colnames(meanSigmaRowResults), "\\\\ \n", sep=" & ")
#     cat("mean ", as.character(dec_three((meanSigmaRowResults)["mean",])), sep=" & ")
#     cat("$\\sigma$", as.character(dec_three((meanSigmaRowResults)["std.dev",])), "\\\\ \n", sep=" & ")

computeExecutionTime()
