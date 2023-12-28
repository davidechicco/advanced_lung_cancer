setwd(".")
options(stringsAsFactors = FALSE)
cat("\014")
set.seed(142)

# R tutorial
# 
# https://datavizpyr.com/how-to-make-umap-plot-in-r/

MISSING_DATA_IMPUTATION <- FALSE

library("pacman")
p_load("tidyverse")
p_load("palmerpenguins")
p_load("umap")
p_load("dplyr")
p_load("mice")
p_load("ROSE")

SAVE_SINGLE_PLOT <- TRUE
if(SAVE_SINGLE_PLOT) cat("The script will save the single plot into a file\n")
SAVE_AGE_SPLIT_PLOT <- TRUE
if(SAVE_AGE_SPLIT_PLOT) cat("The script will save the age split into a file\n")

# PDF_OR_PNG <- "png"
PDF_OR_PNG <- "pdf"

setFontSize <- 20

# Load the patients' data
fileNameData <- "../data/pone0210951_s006_dataset_EDITED_IMPUTED.csv"
patients_data <- read.csv(fileNameData, header = TRUE, sep =",");
patients_data_original <- patients_data

cat("File read: ", fileNameData, "\n", sep="")


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
# patients_data$"Stage" <- NULL

# Small cell lung cancer (SCLC) is usually classified into a two-stage system, limited (LD) and extensive disease (ED
# https://pubmed.ncbi.nlm.nih.gov/12234695/

patients_data[which(patients_data$"histology" == "adenoca"),]$"histology" <- "adenocarcinoma"

patients_data$"histology" %>% unique() %>% sort() %>% cat("\n", sep="\n")
patients_data$"histologyNum" <- -1
patients_data[which(patients_data$"histology" == "adenocarcinoma"),]$"histologyNum" <- 0
patients_data[which(patients_data$"histology" == "large_cell_NEC"),]$"histologyNum" <- 1
patients_data[which(patients_data$"histology" == "mixed"),]$"histologyNum" <- 2
patients_data[which(patients_data$"histology" == "neuroendocrine_tumor"),]$"histologyNum" <- 3
patients_data[which(patients_data$"histology" == "nonsmall_cell"),]$"histologyNum" <- 4
patients_data[which(patients_data$"histology" == "others"),]$"histologyNum" <- 5
patients_data[which(patients_data$"histology" == "poorly_differentiated"),]$"histologyNum" <- 6
patients_data[which(patients_data$"histology" == "sarcomatoid"),]$"histologyNum" <- 7
patients_data[which(patients_data$"histology" == "small_cell"),]$"histologyNum" <- 8
patients_data[which(patients_data$"histology" == "squamous"),]$"histologyNum" <- 9
# patients_data$"histology" <- NULL

patients_data <- patients_data %>% select(order(colnames(patients_data)))


if(MISSING_DATA_IMPUTATION==TRUE){

    # missing data imputation
    NUM_DATASETS <- 1
    imputed_data <- mice(patients_data, m=NUM_DATASETS, maxit = 50, method = 'pmm', seed = 500)
    patients_data <- complete(imputed_data, NUM_DATASETS)
}


random_numbers <- 1
lower_limit_for_random <- 1
upper_limit_for_random <- 10000
random_number <- sample(lower_limit_for_random:upper_limit_for_random, random_numbers)

# Let's do the  analysis for each visit's subset
patients_data_by_visit <- patients_data

cat("Table information: ", nrow(patients_data_by_visit), " rows and ", ncol(patients_data_by_visit), " columns\n", sep="")

# To perform UMAP using Palmer Penguin’s dataset, we will use numerical columns and ignore non-numerical columns as meta data (like we did it for doing tSNE analysis in R). First, let us remove any missing data and add unique row ID.
# penguins <- penguins %>% drop_na() %>% mutate(ID=row_number()) 

# identify columns with only identical values
list_columns_with_identical_values <-  apply(patients_data_by_visit, 2, function(a) length(unique(a))==1) %>% as.data.frame()
colnames(list_columns_with_identical_values)[1] <- "identical_values" 
list_columns_with_identical_values[which(list_columns_with_identical_values$"identical_values" == TRUE),]
list_columns_with_identical_values$"feature" <- rownames(list_columns_with_identical_values)
columns_with_identical_values <- list_columns_with_identical_values[which(list_columns_with_identical_values$identical_values == TRUE),]$"feature"

if(columns_with_identical_values %>% length() > 0) {
    cat("The columns ", sep="") 
    cat(columns_with_identical_values, sep=" ") 
    cat(" have only identical values, so they will be removed from the UMAP application\n", sep="" )
}

patients_data_temp <- patients_data_by_visit %>% drop_na() %>% select(-all_of(columns_with_identical_values))

patients_data_processed <- patients_data_temp %>% drop_na() %>% mutate(ID=row_number()) 

# Let us create a dataframe with all categorical variables with the unique row ID.
# penguins_meta <- penguins %>%  select(ID, species, island, sex)

# Let us select numerical columns using is.numeric() function with select(), standardise the data using scale() function before applying umap() function to perform UMAP.
# umap_fit <- penguins %>% select(where(is.numeric)) %>%  column_to_rownames("ID") %>% scale() %>%  umap()


thisNeighborsNumber <- 15
if(patients_data_processed %>% nrow() <= thisNeighborsNumber) thisNeighborsNumber <- round(patients_data_processed %>% nrow() /2)


umap_fit_patients <- patients_data_processed %>% select(where(is.numeric)) %>%  column_to_rownames("ID") %>% scale() %>%  umap(n_neighbors = thisNeighborsNumber)

# The umap result object is a list object and the layout variable in the list contains two umap components that we are interested in. We can extract the components and save it in a dataframe. Also, we merge the UMAP components with the meta data associated with the data.  
# umap_df <- umap_fit$"layout" %>% 
# as.data.frame() %>% 
# rename(UMAP1="V1",UMAP2="V2") %>% 
# mutate(ID=row_number())%>%
# inner_join(penguins_meta, by="ID")

umap_df_patients <- umap_fit_patients$"layout" %>% 
    as.data.frame() %>% 
    rename(UMAP1="V1",UMAP2="V2") %>% 
    mutate(ID=row_number())%>%
    inner_join(patients_data_processed, by="ID") 

# We can make UMAP plot, a scatter plot with the two UMAP components colored by variables of interest that are part of the data. In this example, we have added color by species variable and shape by sex variable.
# umap_df %>%
#   ggplot(aes(x = UMAP1, 
#  y = UMAP2, 
#  color = species,
#  shape = sex))+
#   geom_point()+
#   labs(x = "UMAP1",
#y = "UMAP2",
#   subtitle = "UMAP plot")

pointSize <- 10
umap_df_patients$"patient_sex" <- ""
umap_df_patients[which(umap_df_patients$"sex"==0),]$"patient_sex" <- "woman"
umap_df_patients[which(umap_df_patients$"sex"==1),]$"patient_sex" <- "man"

umap_df_patients$"sepsis_status" <- ""
umap_df_patients[which(umap_df_patients$"Sepsis"==0),]$"sepsis_status" <- "no_sepsis"
umap_df_patients[which(umap_df_patients$"Sepsis"==1),]$"sepsis_status" <- "yes_sepsis"

plot_title <- paste0("UMAP plot for the advanced lung cancer and sepsis dataset")

first_chosen_variable <- "sepsis_status"
second_chosen_variable <- "patient_sex"

umap_plot_single <- umap_df_patients %>%
    ggplot(aes(x = UMAP1, 
    y = UMAP2, 
    color = get(first_chosen_variable),
    shape = get(second_chosen_variable))) +
    geom_point(size=pointSize, alpha=0.5)+
    labs(x = "UMAP1",
    y = "UMAP2",
    subtitle = plot_title) + 
    theme(text=element_text(size=setFontSize))

plotWidth <- 20
plotHeight <- 10
if(SAVE_SINGLE_PLOT == TRUE) {  


plotFileName <- paste0("../results/UMAP_plot_", first_chosen_variable, "_and_", second_chosen_variable, "_single_plot_rand", random_number, ".", PDF_OR_PNG)
ggsave(umap_plot_single, file=plotFileName, width = plotWidth, height = plotHeight)
cat("saved ", plotFileName, " file\n", sep="")
}

# In the second example of UMAP plot, we have used the same UMAP components, but this time we have added facetting based on island variable to see the relationship between species and island more clearly.

# umap_df %>%
#   ggplot(aes(x = UMAP1, 
#  y = UMAP2,
#  color = species)) +
#   geom_point(size=3, alpha=0.5)+
#   facet_wrap(~island)+
#   labs(x = "UMAP1",
#y = "UMAP2",
#subtitle="UMAP plot")+
#   theme(legend.position="bottom")


thisMin <- summary(umap_df_patients$age)[[1]]  %>% round()
thisFirstQuarter <-  summary(umap_df_patients$age)[[2]] %>% round()
thisSecondQuarter  <- summary(umap_df_patients$age)[[4]]   %>% round()
thisThirdQuarter <-  summary(umap_df_patients$age)[[5]] %>% round()
thisMax <- summary(umap_df_patients$age)[[6]]   %>% round()

umap_df_patients$"age_cate" <- 0
if(umap_df_patients[which(umap_df_patients$"age" >= thisMin & umap_df_patients$"age"<  thisFirstQuarter),] %>% nrow() >= 1)  umap_df_patients[which(umap_df_patients$"age" >= thisMin & umap_df_patients$"age"<  thisFirstQuarter),]$"age_cate" <- paste0("age", thisMin, "_", toString(thisFirstQuarter-1))
if(umap_df_patients[which(umap_df_patients$"age" >= thisFirstQuarter & umap_df_patients$"age"<  thisSecondQuarter),] %>% nrow() >= 1)  umap_df_patients[which(umap_df_patients$"age" >= thisFirstQuarter & umap_df_patients$"age"<  thisSecondQuarter),]$"age_cate" <- paste0("age", thisFirstQuarter, "_", toString(thisSecondQuarter-1))
if(umap_df_patients[which(umap_df_patients$"age" >= thisSecondQuarter & umap_df_patients$"age"<  thisThirdQuarter),] %>% nrow() >= 1)  umap_df_patients[which(umap_df_patients$"age" >= thisSecondQuarter & umap_df_patients$"age"<  thisThirdQuarter),]$"age_cate" <- paste0("age", thisSecondQuarter, "_", toString(thisThirdQuarter-1))
if(umap_df_patients[which(umap_df_patients$"age" >= thisThirdQuarter & umap_df_patients$"age"<=  thisMax),] %>% nrow() >= 1) umap_df_patients[which(umap_df_patients$"age" >= thisThirdQuarter & umap_df_patients$"age"<=  thisMax),]$"age_cate" <- paste0("age", thisThirdQuarter, "_", thisMax)

stage_colors <- c("ED" = "red", "IIIB" = "forestgreen", "IV" = "blue")

# theme_set(theme_bw(18))


first_chosen_variable <- "Stage"
second_chosen_variable <- "sepsis_status"

umap_plot_multi <- umap_df_patients %>%
    ggplot(aes(x = UMAP1, 
    y = UMAP2,
    color = get(first_chosen_variable),
    shape = get(second_chosen_variable))) +
    geom_point(size=pointSize, alpha=0.5)+
    facet_wrap(~age_cate)+
    labs(x = "UMAP1",
    y = "UMAP2",
    subtitle=plot_title)+
    theme(legend.position="left") + geom_jitter() + scale_color_manual(values = stage_colors) + theme(text=element_text(size=setFontSize))


if(SAVE_AGE_SPLIT_PLOT == TRUE) {
    plotFileName <- paste0("../results/UMAP_plot_", first_chosen_variable, "_and_", second_chosen_variable, "_multiple_plot_rand", random_number, ".", PDF_OR_PNG)
    ggsave(umap_plot_multi, file=plotFileName, width = plotWidth, height = plotHeight)
    cat("saved ", plotFileName, " file\n", sep="")
}



first_chosen_variable <- "histology"
second_chosen_variable <- "sepsis_status"
histology_colors <- c("adenocarcinoma" = "red", "large_cell_NEC" = "forestgreen", "mixed" = "blue", "neuroendocrine_tumor" = "yellow", "nonsmall_cell" = "purple", "others" = "black", "poorly_differentiated" = "orange", "sarcomatoid" = "white",  "small_cell" = "pink", "squamous" = "green" )


umap_plot_multiB <- umap_df_patients %>%
    ggplot(aes(x = UMAP1, 
    y = UMAP2,
    color = get(first_chosen_variable),
    shape = get(second_chosen_variable))) +
    geom_point(size=pointSize, alpha=0.5)+
    facet_wrap(~age_cate)+
    labs(x = "UMAP1",
    y = "UMAP2",
    subtitle=plot_title)+
    theme(legend.position="left") + geom_jitter() + scale_color_manual(values = histology_colors) + theme(text=element_text(size=setFontSize))


if(SAVE_AGE_SPLIT_PLOT == TRUE) {
    plotFileName <- paste0("../results/UMAP_plot_", first_chosen_variable, "_and_", second_chosen_variable, "_multiple_plot_rand", random_number, ".", PDF_OR_PNG)
    ggsave(umap_plot_multiB, file=plotFileName, width = plotWidth, height = plotHeight)
    cat("saved ", plotFileName, " file\n", sep="")
}


cat(" : : : The end : : :\n")
