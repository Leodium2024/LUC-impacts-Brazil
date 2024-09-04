library(readxl)
library(writexl)
library(dismo)
library(ggplot2)
library(caret)
library(dplyr)
library(usdm)
library(pROC)
library(flexsdm)
library(sf)
library(spdep)
library(pROC) 
library(ROCR)
library(MASS)
library(geosphere)
library(spatstat)
library(modEvA)
library(tidyr)


########################################################################################
########################################################################################
################      Functions to extract data of a species     #######################
########################################################################################
########################################################################################

Get_Species_Ott_id <-function (Absence_data, Presence_data){
  
  Presence_ottid <- unique(Presence_data$ott_id)
  Presence_Species<- data.frame(Presence_ottid)
  Absence_ottid <- unique(Absence_data$ott_id)
  Absence_Species<- data.frame(Absence_ottid)
  #We keep the species only if we have info on its ranges
  Merged <- merge(Presence_Species, Absence_Species, by.x = "Presence_ottid", by.y = "Absence_ottid", all = FALSE)
  
  return (Merged)
  
}
Combine_Absence_Presence_per_ott_id <- function(Absence_data, Presence_data, ott_id, NB_absence) {
  
  #Combine in one dataframe the absence and presence data and Assoicate the LU value based on the observation year
  
  Variable_Presence <- c("ott_id","decimalLatitude","decimalLongitude","year", 
                         "Bio_1","Bio_2","Bio_3", "Bio_4","Bio_5","Bio_6","Bio_7", "Bio_8","Bio_9","Bio_10","Bio_11", "Bio_12","Bio_13","Bio_14","Bio_15", "Bio_16","Bio_17","Bio_18","Bio_19", 
                         "Elevation","Slope",
                         "Forest1990","F1990","Grass1990","Past1990","Other1990","Agri1990","Forest1991","F1991","Grass1991","Past1991","Other1991","Agri1991","Forest1992","F1992","Grass1992","Past1992","Other1992","Agri1992","Forest1993","F1993","Grass1993","Past1993","Other1993","Agri1993","Forest1994","F1994","Grass1994","Past1994","Other1994","Agri1994","Forest1995","F1995","Grass1995","Past1995","Other1995","Agri1995","Forest1996","F1996","Grass1996","Past1996","Other1996","Agri1996","Forest1997","F1997","Grass1997","Past1997","Other1997","Agri1997","Forest1998","F1998","Grass1998","Past1998","Other1998","Agri1998","Forest1999","F1999","Grass1999","Past1999","Other1999","Agri1999",
                         "Forest2000","F2000","Grass2000","Past2000","Other2000","Agri2000","Forest2001","F2001","Grass2001","Past2001","Other2001","Agri2001","Forest2002","F2002","Grass2002","Past2002","Other2002","Agri2002","Forest2003","F2003","Grass2003","Past2003","Other2003","Agri2003","Forest2004","F2004","Grass2004","Past2004","Other2004","Agri2004","Forest2005","F2005","Grass2005","Past2005","Other2005","Agri2005","Forest2006","F2006","Grass2006","Past2006","Other2006","Agri2006","Forest2007","F2007","Grass2007","Past2007","Other2007","Agri2007","Forest2008","F2008","Grass2008","Past2008","Other2008","Agri2008","Forest2009","F2009","Grass2009","Past2009","Other2009","Agri2009",
                         "Forest2010","F2010","Grass2010","Past2010","Other2010","Agri2010","Forest2011","F2011","Grass2011","Past2011","Other2011","Agri2011","Forest2012","F2012","Grass2012","Past2012","Other2012","Agri2012","Forest2013","F2013","Grass2013","Past2013","Other2013","Agri2013","Forest2014","F2014","Grass2014","Past2014","Other2014","Agri2014","Forest2015","F2015","Grass2015","Past2015","Other2015","Agri2015","Forest2016","F2016","Grass2016","Past2016","Other2016","Agri2016","Forest2017","F2017","Grass2017","Past2017","Other2017","Agri2017","Forest2018","F2018","Grass2018","Past2018","Other2018","Agri2018","Forest2019","F2019","Grass2019","Past2019","Other2019","Agri2019",
                         "Forest2020","F2020","Grass2020","Past2020","Other2020","Agri2020","Forest2021","F2021","Grass2021","Past2021","Other2021","Agri2021","Forest2022","F2022","Grass2022","Past2022","Other2022","Agri2022"
  )
  
  Variable_Absence <- c("ott_id","decimalLatitude","decimalLongitude","year", 
                        "Bio_1","Bio_2","Bio_3", "Bio_4","Bio_5","Bio_6","Bio_7", "Bio_8","Bio_9","Bio_10","Bio_11", "Bio_12","Bio_13","Bio_14","Bio_15", "Bio_16","Bio_17","Bio_18","Bio_19", 
                        "Elevation","Slope",
                        "Forest2021","F2021","Grass2021","Past2021","Other2021","Agri2021"
  )
  #1) Filter to only keep the ott_id researched
  
  # Filter Absence_data
  Absence_filtered <- Absence_data[
    Absence_data$ott_id == ott_id & 
      complete.cases(Absence_data[, Variable_Absence]), 
    ]
  
  # Filter Presence_data
  Presence_filtered <- Presence_data[
    Presence_data$ott_id == ott_id & 
      complete.cases(Presence_data[, Variable_Presence]), 
    ]
  
  #2) Stop the function if we don't have absence data (no info on species range)
  if (nrow(Absence_filtered) == 0) {
    stop("No absence data")
  }
  
  #3) Randomly select the as much absence data than presence data (if not enough absence data, we take all the absence data available)
  num_rows_presence <- nrow(Presence_filtered)
  num_rows_absence <- nrow(Absence_filtered)
  
  if (num_rows_presence <= NB_absence && num_rows_absence >= NB_absence) {
    Absence_filtered <- Absence_filtered[sample(nrow(Absence_filtered), NB_absence), ]
  } else {
    Absence_filtered <- Absence_filtered
  }
  
  #4) creation of a binary absence_Presence variable (1= presence; 0 = absence)
  Absence_filtered$Presence_Absence <- 0
  Presence_filtered$Presence_Absence <- 1
  
  #5) Selection of the variable of interest and merge of absence presence
  
  Absence_final <- Absence_filtered[, c(Variable_Absence, "Presence_Absence")]
  Presence_final <- Presence_filtered[, c(Variable_Presence, "Presence_Absence")]
  Final_data <- bind_rows(Absence_final, Presence_final)
  
  #6) Selection of the land use year selection
  
  Final_data<- Get_LU_year(Final_data)
  
  #7) Weight 
  
  Final_data<- Get_weight(Final_data)
  
  return(as.data.frame(Final_data))
}
Get_LU_year <- function (data) {
  
  # This function associate create the LU variable depending on the year of observation
  # This function is called in the "Combine_Absence_Presence_per_ott_id" funtion
  
  # 1) Add new empty columns
  data$Forest <- NA
  data$Forestry <- NA
  data$Agriculture <- NA
  data$Pasture <- NA
  data$Grassland <- NA
  data$Other <- NA
  
  # 2) Assign values based on the year
  # Remove rows where the year is under 1990 
  data <- subset(data, year >= 1990)
  data <- subset(data, year <= 2022)
  
  # Assign values for years between 1990 and 2021
  for (year in 1990:2022) {
    data$Forest <- ifelse(data$year == year, data[[paste0("Forest", year)]], data$Forest)
    data$Forestry <- ifelse(data$year == year, data[[paste0("F", year)]], data$Forestry)
    data$Agriculture <- ifelse(data$year == year, data[[paste0("Agri", year)]], data$Agriculture)
    data$Pasture <- ifelse(data$year == year, data[[paste0("Past", year)]], data$Pasture)
    data$Grassland <- ifelse(data$year == year, data[[paste0("Grass", year)]], data$Grassland)
    data$Other <- ifelse(data$year == year, data[[paste0("Other", year)]], data$Other)
  }
  
  
  
  
  # 3 Remove the original columns
  data <- data[, !(names(data) %in% c("Forest1990","F1990","Grass1990","Past1990","Other1990","Agri1990","Forest1991","F1991","Grass1991","Past1991","Other1991","Agri1991","Forest1992","F1992","Grass1992","Past1992","Other1992","Agri1992","Forest1993","F1993","Grass1993","Past1993","Other1993","Agri1993","Forest1994","F1994","Grass1994","Past1994","Other1994","Agri1994","Forest1995","F1995","Grass1995","Past1995","Other1995","Agri1995","Forest1996","F1996","Grass1996","Past1996","Other1996","Agri1996","Forest1997","F1997","Grass1997","Past1997","Other1997","Agri1997","Forest1998","F1998","Grass1998","Past1998","Other1998","Agri1998","Forest1999","F1999","Grass1999","Past1999","Other1999","Agri1999",
                                      "Forest2000","F2000","Grass2000","Past2000","Other2000","Agri2000","Forest2001","F2001","Grass2001","Past2001","Other2001","Agri2001","Forest2002","F2002","Grass2002","Past2002","Other2002","Agri2002","Forest2003","F2003","Grass2003","Past2003","Other2003","Agri2003","Forest2004","F2004","Grass2004","Past2004","Other2004","Agri2004","Forest2005","F2005","Grass2005","Past2005","Other2005","Agri2005","Forest2006","F2006","Grass2006","Past2006","Other2006","Agri2006","Forest2007","F2007","Grass2007","Past2007","Other2007","Agri2007","Forest2008","F2008","Grass2008","Past2008","Other2008","Agri2008","Forest2009","F2009","Grass2009","Past2009","Other2009","Agri2009",
                                      "Forest2010","F2010","Grass2010","Past2010","Other2010","Agri2010","Forest2011","F2011","Grass2011","Past2011","Other2011","Agri2011","Forest2012","F2012","Grass2012","Past2012","Other2012","Agri2012","Forest2013","F2013","Grass2013","Past2013","Other2013","Agri2013","Forest2014","F2014","Grass2014","Past2014","Other2014","Agri2014","Forest2015","F2015","Grass2015","Past2015","Other2015","Agri2015","Forest2016","F2016","Grass2016","Past2016","Other2016","Agri2016","Forest2017","F2017","Grass2017","Past2017","Other2017","Agri2017","Forest2018","F2018","Grass2018","Past2018","Other2018","Agri2018","Forest2019","F2019","Grass2019","Past2019","Other2019","Agri2019",
                                      "Forest2020","F2020","Grass2020","Past2020","Other2020","Agri2020","Forest2021","F2021","Grass2021","Past2021","Other2021","Agri2021","Forest2022","F2022","Grass2022","Past2022","Other2022","Agri2022"
  ))]
  
  # 4) data to keep only the rows where all values in the columns specified by cols_to_check are non-missing and non-negative.
  cols_to_check <- c('Forest', 'Pasture', 'Grassland', 'Agriculture', 'Forestry',"Bio_1","Bio_4", "Bio_12","Bio_15", "Elevation", "Slope")
  data <- data[complete.cases(data[, cols_to_check]) & !rowSums(data[, cols_to_check] < 0), ]
  
  
  return (data)
}
Get_weight <- function(data) {
  # Calculate the number of presence and absence
  N_Pres <- sum(data$Presence_Absence == 1)
  N_Absence <- sum(data$Presence_Absence == 0)
  
  # Calculate the weights
  Weight_Presence <- (N_Pres+N_Absence)/N_Pres 
  Weight_Absence <-  (N_Pres+N_Absence)/N_Absence
  
  # Assign weights based on Presence_Absence values
  data$weights <- ifelse(data$Presence_Absence == 1, Weight_Presence, Weight_Absence)
  
  return(data)
}
check_not_all_zeros <- function(df, col_names) {
  # Initialize an empty list to store column names that have at least one non-zero value
  cols_with_non_zeros <- character()
  
  # Iterate over each column name provided in the list
  for (col_name in col_names) {
    # Check if the column exists in the dataframe
    if (col_name %in% names(df)) {
      # Check if there is at least one non-zero value in the column
      if (any(df[[col_name]] != 0)) {
        # If there is at least one non-zero value, add the column name to the list
        cols_with_non_zeros <- c(cols_with_non_zeros, col_name)
      }
    }
  }
  
  return(cols_with_non_zeros)
}


########################################################################################
########################################################################################
##########################      Calibration Functions     ##############################
########################################################################################
########################################################################################

check_multicollinearity <- function(data, predictor, threshold = 10) {
  
  # Use vifstep to check multicollinearity
  Variable_to_test <- data[, predictor]
  result <- usdm::vifstep(Variable_to_test, th = threshold)
  
  return(result)
}
Variable_Selection <- function (data, predictor, response){
  
  Predictor_Response <- c(predictor, response)
  Variable_to_test <- data[, Predictor_Response]
  multicollinearity_result<- check_multicollinearity(data, predictor)
  Remaining = exclude(Variable_to_test, multicollinearity_result)
  
  return(Remaining)
  
}
Extract_coefficients <- function(model, Variable_names, ott_id) {
  
  # Add Intercept to the list of Variable_names
  #Variable_names <- c("Intercept", Variable_names)
  
  # Extract coefficients from the model
  coefficients <- coef(model)
  Coefficients <- data.frame(coefficients, Variable_names)
  
  
  # Extract coefficients for each column
  Forest <- ifelse("Forest" %in% Coefficients$Variable_names, Coefficients$coefficients[Coefficients$Variable_names == "Forest"], 0)
  Pasture <- ifelse("Pasture" %in% Coefficients$Variable_names, Coefficients$coefficients[Coefficients$Variable_names == "Pasture"], 0)
  Grassland <- ifelse("Grassland" %in% Coefficients$Variable_names, Coefficients$coefficients[Coefficients$Variable_names == "Grassland"], 0)
  Agriculture <- ifelse("Agriculture" %in% Coefficients$Variable_names, Coefficients$coefficients[Coefficients$Variable_names == "Agriculture"], 0)
  Forestry <- ifelse("Forestry" %in% Coefficients$Variable_names, Coefficients$coefficients[Coefficients$Variable_names == "Forestry"], 0)
  Bio_1 <- ifelse("Bio_1" %in% Coefficients$Variable_names, Coefficients$coefficients[Coefficients$Variable_names == "Bio_1"], 0)
  Bio_4 <- ifelse("Bio_4" %in% Coefficients$Variable_names, Coefficients$coefficients[Coefficients$Variable_names == "Bio_4"], 0)
  Bio_12 <- ifelse("Bio_12" %in% Coefficients$Variable_names, Coefficients$coefficients[Coefficients$Variable_names == "Bio_12"], 0)
  Bio_15 <- ifelse("Bio_15" %in% Coefficients$Variable_names, Coefficients$coefficients[Coefficients$Variable_names == "Bio_15"], 0)
  Elevation <- ifelse("Elevation" %in% Coefficients$Variable_names, Coefficients$coefficients[Coefficients$Variable_names == "Elevation"], 0)
  Slope <- ifelse("Slope" %in% Coefficients$Variable_names, Coefficients$coefficients[Coefficients$Variable_names == "Slope"], 0)
  Intercept <- ifelse("(Intercept)" %in% Coefficients$Variable_names, Coefficients$coefficients[Coefficients$Variable_names == "(Intercept)"], 0)
  #Intercept <- ifelse("Intercept" %in% Coefficients$Variable_names, Coefficients$coefficients[Coefficients$Variable_names == "Intercept"], 0)
  
  # Check if all coefficients (excluding Intercept) are zero
  if (Forest == 0 && Pasture == 0 && Grassland == 0 && Agriculture == 0 &&
      Forestry == 0 && Bio_1 == 0 && Bio_4 == 0 && Bio_12 == 0 &&
      Bio_15 == 0 && Elevation == 0 && Slope == 0) {
    return(NULL)
  }
  
  # Construct the new row
  new_row <- c(ott_id, Forest, Pasture, Grassland, Agriculture, Forestry, 
               Bio_1,Bio_4,Bio_12,Bio_15, Elevation,Slope, Intercept)
  
  return(new_row)
}
Stepwise_Selection <- function (data,response, predictors){
  
  formula <- paste(response, "~", paste(predictors, collapse = " + "))
  
  base.model <- glm(Presence_Absence ~ 1, data = data, family = binomial, weights = data$weights)
  scope.model <- glm(formula, data = data, family = binomial, weights = data$weights)
  step.model <- stepAIC(base.model, direction = "both", scope = list(lower = base.model, upper = scope.model), trace = FALSE)
  
  
  return(step.model)
  
}

########################################################################################
########################################################################################
############################      Validation Functions   ###############################
########################################################################################
########################################################################################


find_optimal_threshold <- function(model, data) {
  # Predict probabilities
  probabilities <- predict(model, newdata = data, type = "response")
  
  # Extract actual values
  actual <- data$Presence_Absence
  
  # Iterate over a sequence of thresholds and calculate TSS for each
  thresholds <- seq(0, 1, by = 0.001)
  tss_values <- sapply(thresholds, function(t) calculate_tss(probabilities, actual, t))
  
  # Find the threshold with the maximum TSS
  optimal_threshold <- thresholds[which.max(tss_values)]
  
  return(optimal_threshold)
}
calculate_tss <- function(probabilities, actual, threshold) {
  predicted <- ifelse(probabilities >= threshold, 1, 0)
  confusion_matrix <- table(predicted, actual)
  
  # Ensure all entries in confusion matrix exist
  TP <- ifelse("1" %in% rownames(confusion_matrix) && "1" %in% colnames(confusion_matrix), confusion_matrix["1", "1"], 0)
  TN <- ifelse("0" %in% rownames(confusion_matrix) && "0" %in% colnames(confusion_matrix), confusion_matrix["0", "0"], 0)
  FP <- ifelse("1" %in% rownames(confusion_matrix) && "0" %in% colnames(confusion_matrix), confusion_matrix["1", "0"], 0)
  FN <- ifelse("0" %in% rownames(confusion_matrix) && "1" %in% colnames(confusion_matrix), confusion_matrix["0", "1"], 0)
  
  sensitivity <- ifelse(TP + FN > 0, TP / (TP + FN), 0)
  specificity <- ifelse(TN + FP > 0, TN / (TN + FP), 0)
  
  TSS <- sensitivity + specificity - 1
  return(TSS)
}

cross_validate_logistic <- function(data, response, predictors, k = 5, th) {
  # Ensure the response variable is a factor
  data[[response]] <- as.factor(data[[response]])
  
  # Create k equally sized folds
  folds <- createFolds(data[[response]], k = k, list = TRUE, returnTrain = FALSE)
  
  # Initialize metric accumulators
  TPR <- 0
  TNR <- 0
  SORENSEN <- 0
  JACCARD <- 0
  F_meas <- 0
  OR <- 0
  TSS <- 0
  KAPPA <- 0
  AUC <- 0
  CBI <- 0
  IMAE <- 0
  
  # Cross-validation loop
  for (i in 1:k) {
    # Split data into training and test sets
    test_indices <- folds[[i]]
    test_data <- data[test_indices, ]
    train_data <- data[-test_indices, ]
    
    # Apply weights if needed
    test_data <- Get_weight(test_data)
    train_data <- Get_weight(train_data)
    
    # Fit the logistic model on the training data
    model <- Logistic_model(train_data, response, predictors)
    
    # Store the evaluation result
    eval_result <- calculate_metrics(model, test_data, th)
    TPR <- TPR + eval_result[["TPR"]]
    TNR <- TNR + eval_result[["TNR"]]
    SORENSEN <- SORENSEN + eval_result[["SORENSEN"]]
    JACCARD <- JACCARD + eval_result[["JACCARD"]]
    F_meas <- F_meas + eval_result[["F_meas"]]
    OR <- OR + eval_result[["OR"]]
    TSS <- TSS + eval_result[["TSS"]]
    KAPPA <- KAPPA + eval_result[["KAPPA"]]
    AUC <- AUC + eval_result[["AUC"]]
    CBI <- CBI + eval_result[["CBI"]]
    IMAE <- IMAE + eval_result[["IMAE"]]
  }
  
  # Calculate average metrics
  metric <- data.frame(
    TPR = TPR / k,
    TNR = TNR / k,
    SORENSEN = SORENSEN / k,
    JACCARD = JACCARD / k,
    F_meas = F_meas / k,
    OR = OR / k,
    TSS = TSS / k,
    KAPPA = KAPPA / k,
    AUC = AUC / k,
    CBI = CBI / k,
    IMAE = IMAE / k
  )
  
  return(metric)
}
calculate_metrics <- function(model, data, threshold) {
  
  # Predict probabilities
  prob <- predict(model, newdata = data, type = "response")
  predicted <- ifelse(prob> threshold, 1, 0)
  
  # Extract actual values
  actual <- data$Presence_Absence
  # Ensure actual and predicted are factors with the same levels
  actual <- factor(actual, levels = c(0, 1))
  predicted <- factor(predicted, levels = c(0, 1))
  
  # Confusion matrix
  cm <- table(Predicted = predicted, Actual = actual)
  
  TP <- cm["1", "1"]
  TN <- cm["0", "0"]
  FP <- cm["1", "0"]
  FN <- cm["0", "1"]
  
  # True Positive Rate (TPR) or Sensitivity
  TPR <- TP / (TP + FN)
  
  # True Negative Rate (TNR) or Specificity
  TNR <- TN / (TN + FP)
  
  # Sørensen-Dice Coefficient
  SORENSEN <- (2 * TP) / (2 * TP + FP + FN)
  
  # Jaccard Index
  JACCARD <- TP / (TP + FP + FN)
  
  # F-measure on presence-background
  F_meas <- 2 * TP / (2 * TP + FP + FN)
  
  # Omission Rate (OR)
  OR <- FN / (TP + FN)
  
  # True Skill Statistic (TSS)
  TSS <- TPR + TNR - 1
  
  # Kappa statistic
  total <- sum(cm)
  expected <- (sum(cm["1",]) * sum(cm[, "1"]) + sum(cm["0",]) * sum(cm[, "0"])) / total^2
  KAPPA <- (TP + TN - expected) / (total - expected)
  
  # Area Under Curve (AUC)
  ROC<-AUC(obs = actual, pred = prob)
  AUC<-ROC[["AUC"]]
  
  # Boyce Index
  bins <- seq(0, 1, length.out =  101)
  bin_midpoints <- (bins[-1] + bins[-length(bins)]) / 2
  pred_bins <- cut(prob, bins, include.lowest = TRUE)
  observed_freq <- table(pred_bins, actual)[, 2]
  expected_freq <- table(pred_bins)
  ratio <- observed_freq / expected_freq
  boyce_index <- cor(bin_midpoints, ratio, method = "spearman", use = "complete.obs")
  
  # Inverse Mean Absolute Error (IMAE)
  if (!is.null(prob)) {
    mae <- mean(abs(prob - as.numeric(actual)))
    IMAE <- 1 / mae
  } else {
    IMAE <- NA
  }
  
  # Return a list of all metrics
  
  metric <- list(
    TPR = TPR,
    TNR = TNR,
    SORENSEN = SORENSEN,
    JACCARD = JACCARD,
    F_meas = F_meas,
    OR = OR,
    TSS = TSS,
    KAPPA = KAPPA,
    AUC = AUC,
    CBI = boyce_index,
    IMAE = IMAE
  )
  
  return(metric)
} 
calculate_nni <- function(data) {
  data <- data %>% filter(Presence_Absence == 1)
  data <- data %>%
    st_as_sf(coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)
  coords <- st_coordinates(data)
  bbox <- st_bbox(data)
  window <- owin(xrange = c(bbox["xmin"], bbox["xmax"]), yrange = c(bbox["ymin"], bbox["ymax"]))
  points_ppp <- ppp(coords[, 1], coords[, 2], window = window)
  
  nn_distances <- nndist(points_ppp)
  D_observed <- mean(nn_distances)
  lambda <- npoints(points_ppp) / area.owin(window)
  D_expected <- 1 / (2 * sqrt(lambda))
  NNI <- D_observed / D_expected
  return(NNI)
}

Extract_variable_importance<- function(model,ott_id){
  
  importance <- varImp(model_stepwise, normalize=TRUE)
  
  Forest <- ifelse("Forest" %in% names(importance), importance["Forest"], 0)
  Pasture <- ifelse("Pasture" %in% names(importance), importance["Pasture"], 0)
  Grassland <- ifelse("Grassland" %in% names(importance), importance["Grassland"], 0)
  Agriculture <- ifelse("Agriculture" %in% names(importance), importance["Agriculture"], 0)
  Forestry <- ifelse("Forestry" %in% names(importance), importance["Forestry"], 0)
  Bio_1 <- ifelse("Bio_1" %in% names(importance), importance["Bio_1"], 0)
  Bio_4 <- ifelse("Bio_4" %in% names(importance), importance["Bio_4"], 0)
  Bio_12 <- ifelse("Bio_12" %in% names(importance), importance["Bio_12"], 0)
  Bio_15 <- ifelse("Bio_15" %in% names(importance), importance["Bio_15"], 0)
  Elevation <- ifelse("Elevation" %in% names(importance), importance["Elevation"], 0)
  Slope <- ifelse("Slope" %in% names(importance), importance["Slope"], 0)
  
  # Construct the new row
  new_row <- c(ott_id, Forest, Pasture, Grassland, Agriculture, Forestry, 
               Bio_1,Bio_4,Bio_12,Bio_15, Elevation,Slope)
  
  return(new_row)
  
  
}

########################################################################################
########################################################################################
###################################       Command      #################################
########################################################################################
########################################################################################


#Input Data
Absence_data <- read_xlsx("Absence.xlsx")
Presence_data <- read_xlsx("Presence.xlsx")


#Output Data
Calibrating_data <- data.frame()
Coef_df_AIC <- data.frame(matrix(ncol = (length(Predictors_All)+2), nrow = 0))
Variable_importance_df<- data.frame(matrix(ncol = (length(Predictors_All)+2), nrow = 0))
Performance_AIC <- data.frame()

#Predictors
Predictors_All<- c('Forest', 'Pasture', 'Grassland', 'Agriculture', 'Forestry',"Bio_1","Bio_4", "Bio_12","Bio_15", "Elevation", "Slope")
Response <- c('Presence_Absence')

########################################################################################
########################################################################################
###########################      Process the code      #################################
########################################################################################
########################################################################################

#Species included 
Species<- Get_Species_Ott_id(Absence_data, Presence_data) #ott id of the species
num_Species <- (nrow(Species)) # Number of species included

#One iteration per species in which the logistic regression is calibrated and validated
for (i in 1:num_Species) { 
  #ott id of the species
  Ott_id = Species[i,1] 
  print(Ott_id)
  
  #Extract the data for that species
  Abs_Pres = Combine_Absence_Presence_per_ott_id(Absence_data, Presence_data, Ott_id,1000)#Number of absence point has to be determined
  Calibrating_data <- rbind(Calibrating_data, Abs_Pres) #To store the data used for the calibration
  Predictors<-check_not_all_zeros(Abs_Pres,Predictors_All) #Check if the variable have values
  
  #Exclusion  of the variable based on the colinearity
  Abs_Pres_Sorted <- Variable_Selection(Abs_Pres, Predictors, Response)
  Variable_names <- names(Abs_Pres_Sorted)
  Variable_names <- subset(Variable_names, Variable_names != "Presence_Absence")
  
  
  #Calibration including a backward and forward selection based on AIC
  model_stepwise<- Stepwise_Selection(Abs_Pres,Response, Variable_names)
  if (is.null(model_stepwise)) {next}
  
  #Extraction of the coeficients
  summary_model <- summary(model_stepwise)
  coefficients<-coef(summary_model)
  Variable_names_AIC <- rownames(coefficients)
  Coef <- Extract_coefficients(model_stepwise, Variable_names_AIC, Ott_id)
  
  
  #Performances and Varibale importance
  if (is.null(Coef)) {cat("All predictors were eliminated, no model returned.\n")
    next} 
  else {
    cat("Model was successfully refined.\n")
    Coef_df_AIC <- rbind(Coef_df_AIC, Coef) #Store the coeficient of this species with the others
    
    #Find optimal Threshold and validation
    Final_var <- setdiff(Variable_names_AIC, "(Intercept)")
    optimal_threshold <- find_optimal_threshold(model_stepwise, Abs_Pres)
    Evaluation<-cross_validate_logistic(Abs_Pres, Response, Final_var, k = 5, optimal_threshold)
    Evaluation$Ott_id<-Ott_id
    Evaluation$thr_value<-optimal_threshold
    Evaluation$N_Abs<-Abs_Pres %>% filter(Presence_Absence == 0) %>% nrow()
    Evaluation$N_Pre<-Abs_Pres %>% filter(Presence_Absence == 1) %>% nrow()
    Evaluation$AIC<-AIC(model_stepwise)
    Evaluation$BIC<-BIC(model_stepwise)
    Evaluation$NNI<-calculate_nni(Abs_Pres)
    Performance_AIC <- rbind(Performance_AIC,Evaluation)#Store in the performance output varibale
    
    #Calculate variable importance
    Variable_importance<-Extract_variable_importance(model_stepwise,Ott_id)
    Variable_importance_df <- rbind(Variable_importance_df, Variable_importance)
    
  }
  
}

colnames(Coef_df_AIC) <- c('Ott_id', Predictors_All, 'Intercept')
colnames(Variable_importance_df) <- c('Ott_id', Predictors_All)
write_xlsx(Coef_df_AIC,"Species_coef_LR_AIC_1000.xlsx")
write_xlsx(Variable_importance_df,"Variable_importance_LR_AIC_1000.xlsx")
write_xlsx(Performance_AIC,"Performance_LR_AIC_1000.xlsx")
write_xlsx(Calibrating_data,"Calibrating_data_LR_1000.xlsx")
