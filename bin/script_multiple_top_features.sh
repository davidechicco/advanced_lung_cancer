#!/bin/bash
#
#$ -cwd
#$ -S /bin/bash
#
set -o nounset -o pipefail -o errexit
set -o xtrace

random_num=$((1 + $RANDOM % 10000))

# topFeatNum=4
# numIterations=10
# fileName="../results/output_"$random_num"_top"$topFeatNum"features_SMOTE_"$numIterations"iter.txt";
# Rscript random_forest_feature_selection_and_classification_afterwards.r $topFeatNum $numIterations > $fileName 2> $fileName; tail $fileName

grep_lim=1

start=2
end=10
numIterations=100
for((i=start;i<=end;++i));
do
    topFeatNum=$i
    fileName="../results/output_"$random_num"_top"$topFeatNum"features_SMOTE_"$numIterations"iter.txt";
    Rscript random_forest_feature_selection_and_classification_afterwards.r $topFeatNum $numIterations > $fileName 2> $fileName; tail $fileName

#     echo -e "\n"
#     tail -20 "../results/output_5655_top"$i"features.txt" | grep -$grep_lim "MCC"
done


#
# start=11
# end=20
# for((i=start;i<=end;++i));
# do
#     topFeatNum=$i
#     fileName="../results/output_"$random_num"_top"$topFeatNum"features.txt";
#     Rscript random_forest_feature_selection_and_classification_afterwards.r $topFeatNum > $fileName 2> $fileName; tail $fileName
#
#     echo -e "\n"
#     tail -20 "../results/output_2372_top"$i"features.txt" | grep -$grep_lim "MCC"
# done
#
#
# start=21
# end=30
# for((i=start;i<=end;++i));
# do
#     topFeatNum=$i
#     fileName="../results/output_"$random_num"_top"$topFeatNum"features.txt";
#     Rscript random_forest_feature_selection_and_classification_afterwards.r $topFeatNum > $fileName 2> $fileName; tail $fileName
#
#     echo -e "\n"
#     tail -20 "../results/output_625_top"$i"features.txt" | grep -$grep_lim "MCC"
# done



echo "The end"
