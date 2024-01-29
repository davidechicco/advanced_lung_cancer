setwd(".")
options(stringsAsFactors = FALSE)
cat("\014")
set.seed(12)
options(repos = list(CRAN="http://cran.rstudio.com/"))


# agregateTwoSortedRankings
agregateTwoSortedRankings <- function(dd, firstColumnName, secondColumnName) {

    cat("[function agregateTwoSortedRankings()]\n")

    # dd_sorted_MSE <- dd[order(-dd$firstColumnName), ]
    dd_sorted_firstColumn <- dd[order(-dd[[firstColumnName]]), ]
    # print(dd_sorted_firstColumn)
    
    dd_sorted_secondColumn <- dd[order(-dd[[secondColumnName]]), ]
    # print(dd_sorted_IncNodePurity);


    # varImpPlot(rf_output)
    dd_sorted_firstColumn_only <- dd_sorted_firstColumn
    dd_sorted_firstColumn_only[[secondColumnName]] <- NULL # we do not need the other values
    dd_sorted_firstColumn_only$firstColPos <- c(1:dim(dd_sorted_firstColumn_only)[1])
    
    dd_sorted_secondColumn_only <- dd_sorted_secondColumn
    dd_sorted_secondColumn_only[[firstColumnName]] <- NULL # we do not need the other values
    dd_sorted_secondColumn_only$secondColPos <- c(1:dim(dd_sorted_secondColumn_only)[1])

    dd_sorted_firstColumn_only$features <- rownames(dd_sorted_firstColumn_only)
    dd_sorted_secondColumn_only$features <- rownames(dd_sorted_secondColumn_only)

    # let's sort alphabetically
    dd_sorted_firstColumn_only <- dd_sorted_firstColumn_only[order(dd_sorted_firstColumn_only$"features"), ]
    dd_sorted_secondColumn_only <- dd_sorted_secondColumn_only[order(dd_sorted_secondColumn_only$"features"), ]
    
    mergedRanking <- cbind(dd_sorted_firstColumn_only, dd_sorted_secondColumn_only)

    mergedRankingAlphaBeta <- mergedRanking[order(mergedRanking$"features"), ]
    mergedRankingAlphaBeta$posSum <- mergedRankingAlphaBeta$firstColPos + mergedRankingAlphaBeta$secondColPos

    mergedRankingGeneralRank <- mergedRankingAlphaBeta[order(mergedRankingAlphaBeta$"posSum"), ]
    mergedRankingGeneralRank$finalPos <- c(1:dim(mergedRankingGeneralRank)[1])
    
    # remove duplicate columns
    temp <- mergedRankingGeneralRank[, !duplicated(colnames(mergedRankingGeneralRank))]
    mergedRankingGeneralRank <- temp

    # print(mergedRankingGeneralRank)
    
    return (mergedRankingGeneralRank)

}


EXP_ARG_NUM <- 2

args = commandArgs(trailingOnly=TRUE)
if (length(args)<EXP_ARG_NUM) {
  stop("At least one argument must be supplied (input files)", call.=FALSE)
} else {
  # default output file
  TOP_FEATURES_NUMBER <- as.numeric(args[1])
  NUM_ITERATIONS <- as.numeric(args[2])
}

# TOP_FEATURES_NUMBER <- 3

cat("TOP_FEATURES_NUMBER = ", TOP_FEATURES_NUMBER, "\n", sep="")
cat("NUM_ITERATIONS = ", NUM_ITERATIONS, "\n", sep="")

execution_FS_number <- NUM_ITERATIONS
execution_classification_number <- NUM_ITERATIONS # 100


threshold <- 0.5
pROSE <- 0.55

# fileNameData <- "../data/pone0210951_s006_dataset_EDITED.csv"
fileNameData <- "../data/pone0210951_s006_dataset_EDITED_IMPUTED.csv"
targetName <- "Sepsis"


MISSING_DATA_IMPUTATION <- FALSE
TRAIN_SET_OVERSAMPLING_SYNTHETIC_ROSE <- FALSE

if(TRAIN_SET_OVERSAMPLING_SYNTHETIC_ROSE == FALSE) cat("TRAIN_SET_OVERSAMPLING_SYNTHETIC_ROSE == FALSE\n")

TRAIN_SET_OVERSAMPLING_SYNTHETIC_SMOTE <- TRUE

list.of.packages <- c("pacman")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library("pacman")
pacman::p_load("randomForest", "ggplot2", "dplyr", "pastecs",  "mice", "ROSE", "DMwR")

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
patients_data$"Cause_detail_YUJ" <- NULL
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


# removed features because of absent documentation
patients_data$"cancer_status" <- NULL
patients_data$"Cause_Ca_related" <- NULL
patients_data$"Cause_Direct_coding_YUJ" <- NULL
patients_data$"Cr" <- NULL
patients_data$"Death" <- NULL
patients_data$"Dx" <- NULL
patients_data$"GW_Death" <- NULL
patients_data$"hemodynamic" <- NULL
patients_data$"highCr" <- NULL
patients_data$"Ix" <- NULL
patients_data$"K" <- NULL
patients_data$"line" <- NULL
patients_data$"Na" <- NULL
patients_data$"Ninf" <- NULL
patients_data$"NIV" <- NULL
patients_data$"NIVMV" <- NULL
patients_data$"NOVA" <- NULL
patients_data$"pallRT" <- NULL
patients_data$"Planned" <- NULL
patients_data$"post_ICU_HD" <- NULL
patients_data$"prev_op" <- NULL
patients_data$"Procedure" <- NULL
patients_data$"procedure_related" <- NULL
patients_data$"Steroid_ICU" <- NULL
patients_data$"Txoff" <- NULL




patients_data$"StageNum" <- -1
patients_data[patients_data$"Stage" == "ED",]$"StageNum" <- 5
patients_data[patients_data$"Stage" == "IIIB",]$"StageNum" <- 3
patients_data[patients_data$"Stage" == "IV",]$"StageNum" <- 4
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

imbal_number <- 126

# patients_data$"HCV.RNATaqman.Log.IU.ml." <- as.numeric(patients_data$"HCV.RNATaqman.Log.IU.ml.")

# patients_data$"category_0healthy_1hepatitis_2fibrosis_3cirrhorsis" <- NULL
if(MISSING_DATA_IMPUTATION==TRUE){

    # missing data imputation
    NUM_DATASETS <- 1
    imputed_data <- mice(patients_data, m=NUM_DATASETS, maxit = 50, method = 'pmm', seed = 500)
    patients_data <- complete(imputed_data, NUM_DATASETS)
}

aggregateRankings <- NULL


TRAINING_SET_RATIO <- 0.9
TEST_SET_RATIO <- 1 - TRAINING_SET_RATIO

cat("TRAINING_SET_RATIO = ", TRAINING_SET_RATIO, "\n", sep="")

allExecutionsFinalRanking <- data.frame(Doubles=double(),
                 Ints=integer(),
                 Factors=factor(),
                 Logicals=logical(),
                 Characters=character(),
                 stringsAsFactors=FALSE)
                 


cat("Number of classification executions = ", execution_classification_number, "\t", sep="")
for(exe_class_i in 1:execution_classification_number)
{      

    patients_data <- patients_data[sample(nrow(patients_data)),] # shuffle the rows

    patients_training_set_index_start <- 1
    patients_training_set_index_end <- round(nrow(patients_data) * TRAINING_SET_RATIO)
    patients_test_set_index_start <- patients_training_set_index_end + 1
    patients_test_set_index_end <- nrow(patients_data)

    patients_training_set <- patients_data[patients_training_set_index_start:patients_training_set_index_end,]
    patients_test_set <- patients_data[patients_test_set_index_start:patients_test_set_index_end,]

    cat("Number of feature selection executions = ", execution_FS_number, "\n", sep="")
    for(exe_fs_i in 1:execution_FS_number)
    {

        cat("\n\n\n Execution number feature selection = ", exe_fs_i, " & classif =  ", exe_class_i,"\n", sep="")
        cat("[Randomizing the rows]\n")
        patients_training_set <- patients_training_set[sample(nrow(patients_training_set)),] # shuffle the rows

        this_target_index <- patients_training_set %>% ncol()

        
        # formula
        allFeaturesFormula <- as.formula(paste(as.factor(colnames(patients_training_set)[this_target_index]), '.', sep=' ~ ' ))
        if(TRAIN_SET_OVERSAMPLING_SYNTHETIC_ROSE == TRUE)
            {
                cat("ROSE oversampling\n")
                data.rose <- ROSE(allFeaturesFormula, data = patients_training_set, p=pROSE, seed = 1)$data
                patients_training_set <- data.rose
            }

        if(TRAIN_SET_OVERSAMPLING_SYNTHETIC_SMOTE == TRUE)
            {

                this_patients_training_set <- patients_training_set
                this_patients_training_set$"target" <- as.factor(this_patients_training_set$"target")

                # cat("this_patients_training_set %>% str()\n")
                # this_patients_training_set %>% str() %>% print()

                cat("SMOTE oversampling\n")
                newData <- SMOTE(allFeaturesFormula, data = this_patients_training_set, perc.over = imbal_number, perc.under = imbal_number)
                patients_training_set <- newData
            }


        cat("application of randomForest()\n")
        rf_output <- randomForest(as.factor(patients_training_set$target) ~ ., data=patients_training_set, importance=TRUE, proximity=TRUE)
            

        dd <- as.data.frame(rf_output$importance);
        mergedRankingGeneralRank <- agregateTwoSortedRankings(dd, "MeanDecreaseAccuracy", "MeanDecreaseGini")
        
        rownames(mergedRankingGeneralRank) <- (removeDot(removeUnderscore(rownames(mergedRankingGeneralRank))))
        mergedRankingGeneralRank$features <- removeDot(removeUnderscore(mergedRankingGeneralRank$features))
        
        # print(mergedRankingGeneralRank[, c("finalPos", "MeanDecreaseAccuracy", "MeanDecreaseGini"), drop=FALSE])

        finalRankingOneExecution <- mergedRankingGeneralRank[, c("features", "finalPos", "MeanDecreaseAccuracy", "MeanDecreaseGini"), drop=FALSE]
        finalRankingOneExecutionAlphaBeta <- finalRankingOneExecution[order(finalRankingOneExecution$"features"), , drop=FALSE]

        if (exe_fs_i == 1) {
            allExecutionsFinalRanking <- finalRankingOneExecutionAlphaBeta
        } else {
            
            allExecutionsFinalRanking$MeanDecreaseAccuracy <- allExecutionsFinalRanking$MeanDecreaseAccuracy + finalRankingOneExecutionAlphaBeta$MeanDecreaseAccuracy
            allExecutionsFinalRanking$MeanDecreaseGini <- allExecutionsFinalRanking$MeanDecreaseGini + finalRankingOneExecutionAlphaBeta$MeanDecreaseGini
            allExecutionsFinalRanking$finalPos <- allExecutionsFinalRanking$finalPos + finalRankingOneExecutionAlphaBeta$finalPos
        }
    }



    allExecutionsFinalRanking$MeanDecreaseAccuracy <- allExecutionsFinalRanking$MeanDecreaseAccuracy / execution_FS_number
    allExecutionsFinalRanking$MeanDecreaseGini <- allExecutionsFinalRanking$MeanDecreaseGini / execution_FS_number
    allExecutionsFinalRanking$finalPos <- allExecutionsFinalRanking$finalPos / execution_FS_number

    # # let's eliminate the target index from the rank
    # targetRow <-  which(allExecutionsFinalRanking==targetName)
    # allExecutionsFinalRanking <- allExecutionsFinalRanking[-c( which(allExecutionsFinalRanking==targetName)), ]

    cat("\n\n\n\n== final ranking after ", execution_FS_number, " executions == \n", sep="")

    allExecutionsFinalRanking_mse_Gini <-  allExecutionsFinalRanking[, c("MeanDecreaseAccuracy", "MeanDecreaseGini")]
    aggregateRankings <- agregateTwoSortedRankings(allExecutionsFinalRanking_mse_Gini, "MeanDecreaseAccuracy", "MeanDecreaseGini")

    aggregateRankings$"features" <- gsub(" ", "_", aggregateRankings$"features")
    rownames(aggregateRankings) <- gsub(" ", "_", rownames(aggregateRankings))
    rownames(allExecutionsFinalRanking_mse_Gini) <- gsub(" ", "_", rownames(allExecutionsFinalRanking_mse_Gini))

    aggregateRankings$"meanDecreaseAccuracyPosition" <- -1
    aggregateRankings$"meanDecreaseGiniPosition" <- -1
    aggregateRankings[order(-aggregateRankings$"MeanDecreaseAccuracy"),]$"meanDecreaseAccuracyPosition" <- seq(1:nrow(aggregateRankings))
    aggregateRankings[order(-aggregateRankings$"MeanDecreaseGini"),]$"meanDecreaseGiniPosition" <- seq(1:nrow(aggregateRankings))

    aggregateRankings$"aggregatedPosition" <- -1
    aggregateRankings$"aggregatedPosition" <- aggregateRankings$"meanDecreaseAccuracyPosition" + aggregateRankings$"meanDecreaseGiniPosition"
    
    
    # print(aggregateRankings[, c("finalPos", "MeanDecreaseAccuracy", "MeanDecreaseGini")])

    print(aggregateRankings[order(aggregateRankings$"aggregatedPosition"), ])
    # print(allExecutionsFinalRanking_mse_Gini[order(-allExecutionsFinalRanking_mse_Gini["MeanDecreaseAccuracy"]), ])


    SAVE_RANKING_TO_CSV <- TRUE
    if(SAVE_RANKING_TO_CSV) {
        outputFileName <- paste0("../results/feature_ranking_SMOTE_rand", exe_num, ".csv")
        write.csv(aggregateRankings, file=outputFileName, row.names=FALSE)
        cat("saved file ", outputFileName, "\n")
    }

    top_features_num <- TOP_FEATURES_NUMBER
    # selectedFeaturesNames <- rownames((allExecutionsFinalRanking_mse_Gini[order(-allExecutionsFinalRanking_mse_Gini["MeanDecreaseAccuracy"]), ])[1:top_features_num,])
    # selectedFeaturesNames <- rownames((aggregateRankings[order(aggregateRankings["aggregatedPosition"]), ])[1:top_features_num,])
    selectedFeaturesNames <- rownames((aggregateRankings[order(aggregateRankings$"aggregatedPosition"), ])[1:top_features_num,])


    cat("number of selected top features: ", top_features_num, "\n", sep="")
    cat("selected top features: \n")
    print(selectedFeaturesNames)


    patients_training_set_reduced_features <- patients_training_set[, c(selectedFeaturesNames, "target")]
    patients_test_set_reduced_features <- patients_test_set[, c(selectedFeaturesNames, "target")]

    # Allocation of the training set and of the test set
    training_set_with_target <- patients_training_set_reduced_features
    this_target_index <- training_set_with_target %>% ncol()

    # formula
    allFeaturesFormula <- as.formula(paste(as.factor(colnames(patients_training_set_reduced_features)[this_target_index]), '.', sep=' ~ ' ))
	  if(TRAIN_SET_OVERSAMPLING_SYNTHETIC_ROSE == TRUE)
	      {
		  thisP <- pROSE

          cat("ROSE oversampling\n")
		  data.rose <- ROSE(allFeaturesFormula, data = training_set_with_target, p=thisP, seed = 1)$data
          patients_training_set_reduced_features <- data.rose
	      }

    cat("\n[Training Random Forests classifier on the training set with only the top ", top_features_num ," features]\n")
    rf_new <- NULL
    allFeaturesFormula <- as.formula(paste(as.factor("target"), '.', sep=' ~ ' ))
    patients_training_set_reduced_features$"target" <- as.factor(patients_training_set_reduced_features$"target")
    rf_new <- randomForest(allFeaturesFormula, data=patients_training_set_reduced_features, importance=FALSE, proximity=TRUE, type="classification")
	
    cat("\n[Applying the trained Random Forests classifier on the test set with only the top ", top_features_num ," features]\n")
    patients_data_test_PRED <- predict(rf_new, patients_test_set_reduced_features, type="prob")[,"1"]
    patients_data_test_labels <- patients_test_set_reduced_features$"target"
    thisConfMat <- confusion_matrix_rates(patients_data_test_labels, patients_data_test_PRED, "@@@ Test set @@@")

    if (exe_class_i == 1)  { 
        confMatDataFrame <-  thisConfMat 
    } else { 
        confMatDataFrame <- rbind(confMatDataFrame, thisConfMat) 
    }
    
 }
  
 
cat("\n\n\n=== final results ===\n")
cat("Number of executions = ", execution_classification_number, "\n", sep="")

# statistics on the dataframe of confusion matrices
statDescConfMatr <- stat.desc(confMatDataFrame)
meanSigmaRowResults <- (statDescConfMatr)[c("mean","std.dev"),]
print(dec_three(statDescConfMatr))
cat("\n\n=== === === ===\n")
print(dec_three(meanSigmaRowResults))
cat("\n\n=== === === ===\n\n\n")


computeExecutionTime()

if(TRAIN_SET_OVERSAMPLING_SYNTHETIC_ROSE == TRUE) cat("oversampled with pROSE = ", pROSE, "\n", sep="")
if(TRAIN_SET_OVERSAMPLING_SYNTHETIC_SMOTE == TRUE) cat("oversampled with SMOTE with perc.over & perc.under = ", imbal_number, "\n", sep="")


FEATURE_RANKING_PLOT_DEPICTION <- TRUE
if (FEATURE_RANKING_PLOT_DEPICTION == TRUE) {
    
        # print(colnames(dd_sorted_IncNodePurity_only))


        numberOfFeaturesToPlot <- 20
        cat("We'll plot only the top ", numberOfFeaturesToPlot, " features\n", sep="")
        theseAggregateRankings <- aggregateRankings[(1:numberOfFeaturesToPlot),]

        folder <- "../results/"
        mkdirResultsCommand <- paste0("mkdir -p ", folder)
        system(mkdirResultsCommand)
        cat("applied command: ", mkdirResultsCommand, "\n", sep="")
        x_upper_lim <- -1
          
         barPlotOfRanking(theseAggregateRankings, theseAggregateRankings$MeanDecreaseAccuracy, theseAggregateRankings$features, theseAggregateRankings$firstColPos, exe_num, "features", "MeanDecreaseAccuracy", x_upper_lim, folder)
         
         barPlotOfRanking(theseAggregateRankings, theseAggregateRankings$MeanDecreaseGini, theseAggregateRankings$features, theseAggregateRankings$secondColPos, exe_num, "features", "MeanDecreaseGini", x_upper_lim, folder)
            
}
