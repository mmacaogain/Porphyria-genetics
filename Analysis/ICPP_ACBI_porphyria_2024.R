setwd("C:/Users/mmaca/Code/git/Porphyria-genetics/Analysis/")
library(reshape2)
library(knitr)
library(randomForest)
library(caret)
library(pROC)
library(dplyr)
library(tidyverse

data<-read.csv("../Scores.csv")

# Define the thresholds for each predictor
thresholds <- c(
  PP2 = 0.908,  # Example threshold for PolyPhen-2
  SIFT = 0.05,  #Example threshold for SIFT
  Consurf = 5,  # Example threshold (assuming 5 as a placeholder)
  MA = 3.5,  # Example threshold for MutationAssessor
  PROVEAN = -2.5,  # Example threshold for PROVEAN
  FATHMM = -6.0,  # Example threshold for FATHMM
  REVEL = 0.7,  # Example threshold for REVEL
  CADD = 25,  # Example threshold for CADD (PHRED score)
  AlphaMissense_score = 0.5,  # Example threshold for AlphaMissense
  DiscovEHR = 1,  # Treat as 1 if not NA, 0 if NA
  BayesDel_addAF_score_0.07 = 0.0  # Example threshold for BayesDel_addAF
)


# Function to calculate the consensus score for each row
calculate_consensus <- function(row, thresholds) {
  count_above_threshold <- 0
  total_predictors <- 0
  
  for (predictor in names(thresholds)) {
    value <- row[[predictor]]
    
    # Handle non-NA values
    if (!is.na(value)) {
      total_predictors <- total_predictors + 1
      
      # Special handling for DiscovEHR
      if (predictor == "DiscovEHR") {
        if (value >= thresholds[predictor]) {
          count_above_threshold <- count_above_threshold + 1
        }
      } else if (value > thresholds[predictor]) {
        count_above_threshold <- count_above_threshold + 1
      }
    }
    
    # Handle NA values for DiscovEHR (count as 0)
    if (is.na(value) && predictor == "DiscovEHR") {
      total_predictors <- total_predictors + 1
    }
  }
  
  # Calculate the percentage consensus score
  percentage_consensus <- (count_above_threshold / total_predictors) * 100
  
  return(percentage_consensus)
}

# Apply the function to each row in the dataset
data$Consensus_Score <- apply(data, 1, calculate_consensus, thresholds = thresholds)

# View the resulting data frame
head(data)

# Save the updated data to a new CSV file if needed
# write.csv(data, "output_file_with_consensus_score.csv", row.names = FALSE)

#set cut off values:
cutz <- data.frame(predictornames = 
               c("PP2", "SIFT", "Consurf","MA", "PROVEAN","FATHMM", "CONDEL", "REVEL","CADD","AlphaMissense_score", "DiscovEHR", "BayesDel_addAF_score_0.07"))

cutz$cutoffs <- c(0.908, 0.05,  5,  3.5,  -2.5, -6, 0.469,  0.7,  25,  0.5, 0.07)

PP2 = 0.908
SIFT= 0.05
Consurf= 5
MA= 3.5
PROVEAN= -2.5
FATHMM= -6
CONDEL= 0.469
REVEL= 0.7
CADD= 25
AlphaMissense_score= 0.5
DiscovEHR= NA
BayesDel_addAF_score_0.07= 0.07


#Original cuts: "PP2=0.908", "SIFT=0.05", "Consurf=5","MA=3.5", "PROVEAN=-2.5","CONDEL=0.5", "REVEL=0.5","CADD=25","Consensus=7.9".

#graph prediction scores
meltadata<-melt(PredScores)
meltadata$Gene <- factor(meltadata$Gene, c("HMBS", "PPOX","CPOX"))
meltadata$cutoffs <- cutz$cutoffs[match(meltadata$variable, cutz$predictornames )]
submeltadata <- subset(meltadata, variable %in% c("Consensus"))
subsubmeltadata<-subset(submeltadata, Functional_status != "undetermined")
ggplot(subsubmeltadata, aes(x=Functional_status, y=value, group=Functional_status, colour = Functional_status ))+
  scale_colour_manual(name = "Functional status", values = c("purple", "red" , "dark grey"))+
  scale_shape_manual(values = c(19, 19, 1))+
  geom_violin(colour = "black")+
  geom_jitter(width = 0.2, height = 0.1, shape = 1)+
  #facet_wrap(~submeltadata$variable, scales = "free", ncol = 3)+
  theme_classic()+
  theme(legend.position="bottom")+
  labs(x = "Functional status", y = "SJH in silico prediction score")+
  geom_hline(data = subsubmeltadata, aes(yintercept = cutoffs), linetype = 'dotted')

#c("Prot","Gene","PP2","SIFT","Consurf","MA","PROVEAN","CONDEL","REVEL","CADD","Consensus")

PredScoresPred<-subset(PredScores, Functional_status != "undetermined")
PredScoresPred<-subset(PredScoresPred, MA != "NA")

PredScoresPred$y<-ifelse(PredScoresPred$Functional_status == "pathogenic",1,0)
PredScoresPred$cons<-ifelse(PredScoresPred$Consensus >= 8,1,0)

#logistic regression model
PredScoreModel = glm(data = PredScoresPred, y~cons,family=binomial(link=logit))
#minimal
#PredScoreModel = glm(data = PredScoresPred, y~SIFT.05+MA.3.5+PROVEAN..2.5,family=binomial(link=logit))
summary(PredScoreModel) 

# Calculate sensitivity and false positive measures for logit model "Consensus"
fity_ypos <- PredScoreModel$fitted[PredScoresPred$y == 1]
fity_yneg <- PredScoreModel$fitted[PredScoresPred$y == 0]
sort_fity <- sort(PredScoreModel$fitted.values)

sens <- 0
spec_c <- 0

for (i in length(sort_fity):1){
  sens <- c(sens, mean(fity_ypos >= sort_fity[i]))
  spec_c <- c(spec_c, mean(fity_yneg >= sort_fity[i]))
  
} 

# plot ROC curves
plot(spec_c, sens, xlim = c(0, 1), ylim = c(0, 1), type = "l", 
     xlab = "false positive rate", ylab = "true positive rate", col = 'red', lwd = 2)
abline(0, 1, col= "black")

#Confusion matrix: Log.Reg.
threshold=0.5
predicted_values<-ifelse(predict(PredScoreModel,type="response")>threshold,1,0)
actual_values<-PredScoreModel$y
conf_matrix<-table(predicted_values,actual_values)
confusionMatrix(conf_matrix,positive = "1")

#Confusion matrix: RF
threshold=0.5
predicted_values<-ifelse(predict(rf,type="response")>threshold,1,0)
actual_values<-rf$y
conf_matrix<-table(predicted_values,actual_values)
confusionMatrix(conf_matrix,positive = "1")

###confusion matrix: individual variables
predicted_values<-PredScoresPred$cons
actual_values<-PredScoresPred$y
conf_matrix<-table(predicted_values,actual_values)
confusionMatrix(conf_matrix, positive = "1")




# Plot the ROC curve for the best model
plot(varImp(rf_cv_model))  # Variable importance plot

# To visualize the ROC curve using the test set predictions (if you have a test set)
# Assuming you split the data into train/test sets as before:
rf_probs_test <- predict(rf_cv_model, test_data, type = "prob")[, 2]  # Probability for class 1 (Pathogenic)
rf_roc_test <- roc(test_data$Pathogenic_Classification_Binary, rf_probs_test)

# Plot the ROC curve for the test set
plot(rf_roc_test, col = "green", lwd = 2, main = "Random Forest Test ROC Curve")

