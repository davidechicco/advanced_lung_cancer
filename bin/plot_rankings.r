setwd(".")
options(stringsAsFactors = FALSE)
cat("\014")
# set.seed(11)

rankingFile <- "../results/feature_ranking_RF_results00.csv"
aggregateRankings <- read.csv(rankingFile, header = TRUE, sep =",");
cat("Read data from file ", rankingFile, "\n", sep="")

rownames(aggregateRankings) <- gsub(" ", "_", rownames(aggregateRankings))

aggregateRankings$"meanDecreaseAccuracyPosition" <- -1
aggregateRankings$"meanDecreaseGiniPosition" <- -1
aggregateRankings[order(-aggregateRankings$meanDecreaseAccuracy),]$"meanDecreaseAccuracyPosition" <- seq(1:nrow(aggregateRankings))
aggregateRankings[order(-aggregateRankings$meanDecreaseGiniScore),]$"meanDecreaseGiniPosition" <- seq(1:nrow(aggregateRankings))

aggregateRankings$"aggregatedPosition" <- -1
aggregateRankings$"aggregatedPosition" <- aggregateRankings$"meanDecreaseAccuracyPosition" + aggregateRankings$"meanDecreaseGiniPosition"

generalAggregateRankings <- aggregateRankings[order(aggregateRankings$aggregatedPosition),]

x_upper_lim <- -1

resultsDirPath <- "../results/"

ratioTopFeaturesToShow <- 0.2
topFeaturesNum <- round(generalAggregateRankings %>% nrow() * ratioTopFeaturesToShow, 0)

selectedSubDataFrame <- generalAggregateRankings[c(1:topFeaturesNum),]

barPlotOfRanking(selectedSubDataFrame, selectedSubDataFrame$"meanDecreaseAccuracy", selectedSubDataFrame$"feature", selectedSubDataFrame$"meanDecreaseAccuracyPosition", exe_num, "feature", "meanDecreaseAccuracy", x_upper_lim, resultsDirPath)
barPlotOfRanking(selectedSubDataFrame, selectedSubDataFrame$"meanDecreaseGiniScore", selectedSubDataFrame$"feature", selectedSubDataFrame$"meanDecreaseGiniPosition", exe_num, "feature", "meanDecreaseGiniScore", x_upper_lim, resultsDirPath)
