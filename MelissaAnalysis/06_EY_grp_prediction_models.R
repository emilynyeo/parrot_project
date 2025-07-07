# Making a good classification model 

#### Load libraries ####
pacman::p_load(tidyverse, vegan, MASS, phyloseq, tibble, ANCOMBC, ggplot2, coin,
               MicrobiotaProcess, patchwork, reshape2, ggnewscale, VennDiagram,
               UpSetR, gridExtra, grid, WGCNA, gt, caret, randomForest, pROC, e1071)
library(e1071)      # For confusionMatrix
#BiocManager::install(c("GO.db", "impute", "preprocessCore"))

#### Load data ####
# Also converting to matrices with rownames for phyloseq
setwd("/Users/emily/projects/research/parrot_project/MelissaAnalysis")
load("01_process_and_clean_data/phyloseq16n.rds")
load("01_process_and_clean_data/phyloseq16ren.rds")

load("01_process_and_clean_data/phyloseq18n.rds")
load("01_process_and_clean_data/phyloseq18ren_nohost.rds")
load("01_process_and_clean_data/phyloseq18m.rds")

# Split into classes phylo
load("02_split_18_classes/phyloseq18nano_nohost.rds")
load("02_split_18_classes/phyloseq18nano_plants.rds")
load("02_split_18_classes/phyloseq18nano_microeuk.rds")

load("02_split_18_classes/phyloseq18miseq_nohost.rds")
load("02_split_18_classes/phyloseq18miseq_plants.rds")
load("02_split_18_classes/phyloseq18miseq_microeuk.rds")

# ** Merge phyloseq objects **
phyloseq16all <- merge_phyloseq(phyloseq16n, phyloseq16ren)
phyloseq18all <- phyloseq18nano_nohost
phyloseq18all <- prune_samples(sample_sums(phyloseq18all) >0, phyloseq18all)
phyloseq18miseq_nohost <- prune_samples(sample_sums(phyloseq18miseq_nohost) >0, phyloseq18miseq_nohost)

# ** Combine Crate Data into one **
sample_df <- as.data.frame(sample_data(phyloseq16all))

# Create new variable with recoded values
sample_df$cap_wild_crate <- as.character(sample_df$Captive.Wild)
sample_df$cap_wild_crate[sample_df$Captive.Wild %in% c("Crate 1: Unknown origin", 
                                                       "Crate 2: Unknown origin")] <- "New crates"
# factor with levels
sample_df$cap_wild_crate <- factor(sample_df$cap_wild_crate,
                                   levels = c("Captive", 
                                              "Wild, free ranging", 
                                              "Wild, seized from traffickers", 
                                              "New crates"))
# Replace with updated data
sample_data(phyloseq16all) <- sample_data(sample_df)
levels(sample_data(phyloseq16all)$cap_wild_crate)

# Create new phyloseq object without "New crates"
phyloseq_no_new_crates <- subset_samples(phyloseq16all, cap_wild_crate != "New crates")
levels(sample_data(phyloseq_no_new_crates)$cap_wild_crate)


### MANUAL MODEL ### 
ps <- phyloseq16all

# Transform to relative abundance
ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))
otu <- as.data.frame(otu_table(ps_rel))
meta <- as.data.frame(sample_data(ps_rel))

# Combine OTU and metadata
otu_t <- as.data.frame(t(otu))
data <- cbind(otu_t, Group = meta$cap_wild_crate)
data$Group <- as.factor(data$Group)
levels(data$Group) <- c("Captive", "Wild_free", "Wild_seized", "New_crates")

#  Create Train/Test Split
set.seed(123)
trainIndex <- createDataPartition(data$Group, p = 0.8, list = FALSE)
trainData <- data[trainIndex, ]
testData <- data[-trainIndex, ]

# Caret k-fold cross-validation (10-fold) & model training (Random Forest here):
ctrl <- trainControl(method = "cv", number = 10,
                     classProbs = TRUE,
                     summaryFunction = multiClassSummary,
                     savePredictions = TRUE)

model <- train(Group ~ ., data = trainData,
               method = "rf",
               trControl = ctrl,
               metric = "Accuracy")

# Evaluate Model on Test Set
pred <- predict(model, newdata = testData)
confusion <- confusionMatrix(pred, testData$Group)
print(confusion)

# per-class sensitivity and specificity:
confusion$byClass[, c("Sensitivity", "Specificity")]

# macro/micro averages:
mean(confusion$byClass[, "Sensitivity"], na.rm = TRUE)
mean(confusion$byClass[, "Specificity"], na.rm = TRUE)

cm_df <- as.data.frame(confusion$table)
colnames(cm_df) <- c("Prediction", "Reference", "Freq")

confusion_plot <- ggplot(cm_df, aes(x = Reference, y = Prediction, fill = Freq)) +
  geom_tile() +
  geom_text(aes(label = Freq), color = "white", size = 5) +
  scale_fill_gradient(low = "steelblue", high = "red") +
  labs(title = "Confusion Matrix Heatmap", x = "Actual", y = "Predicted") +
  theme_minimal()

importance <- varImp(model)$importance
importance$Feature <- rownames(importance)

importance_plot <- ggplot(importance, aes(x = reorder(Feature, Overall), y = Overall)) +
  geom_col(fill = "darkgreen") +
  coord_flip() +
  labs(title = "Variable Importance", x = "Feature", y = "Importance") +
  theme_minimal()

#### MANUAL MICROBIOME PREDICTION MODEL ####

# function for easy tuning 

run_microbiome_model <- function(physeq_obj,
                                 group_var,
                                 method = "rf",
                                 p = 0.8,
                                 cv_folds = 10,
                                 seed = 123,
                                 plot_results = TRUE) {
  
  set.seed(seed)
  # Transform to relative abundance
  ps_rel <- transform_sample_counts(physeq_obj, function(x) x / sum(x))
  # Extract OTU and metadata
  otu <- as.data.frame(otu_table(ps_rel))
  meta <- as.data.frame(sample_data(ps_rel))
  
  # Ensure OTUs are samples in rows
  if (!taxa_are_rows(ps_rel)) {
    otu_t <- as.data.frame(otu)
  } else {
    otu_t <- as.data.frame(t(otu))
  }
  
  # Combine with metadata
  data <- cbind(otu_t, Group = meta[[group_var]])
  data$Group <- as.factor(data$Group)
  data$Group <- as.factor(make.names(data$Group))
  
  # Create Train/Test Split
  trainIndex <- createDataPartition(data$Group, p = p, list = FALSE)
  trainData <- data[trainIndex, ]
  testData <- data[-trainIndex, ]
  
  # Train Control
  ctrl <- trainControl(method = "cv",
                       number = cv_folds,
                       classProbs = TRUE,
                       summaryFunction = multiClassSummary,
                       savePredictions = TRUE)
  
  # Model training
  model <- train(Group ~ ., data = trainData,
                 method = method,
                 trControl = ctrl,
                 metric = "Accuracy")
  
  # Test set prediction
  pred <- predict(model, newdata = testData)
  
  # Confusion matrix
  confusion <- confusionMatrix(pred, testData$Group)
  
  # Sensitivity and specificity (per class and averages)
  sens_spec <- confusion$byClass[, c("Sensitivity", "Specificity")]
  mean_sens <- mean(confusion$byClass[, "Sensitivity"], na.rm = TRUE)
  mean_spec <- mean(confusion$byClass[, "Specificity"], na.rm = TRUE)
  
  confusion_plot <- NULL
  importance_plot <- NULL
  
  if (plot_results) {
    # Confusion Matrix Heatmap
    cm_df <- as.data.frame(confusion$table)
    colnames(cm_df) <- c("Prediction", "Reference", "Freq")
    
    confusion_plot <- ggplot(cm_df, aes(x = Reference, y = Prediction, fill = Freq)) +
      geom_tile() +
      geom_text(aes(label = Freq), color = "white", size = 5) +
      scale_fill_gradient(low = "steelblue", high = "red") +
      labs(title = "Confusion Matrix Heatmap", x = "Actual", y = "Predicted") +
      theme_minimal()
    
    # Variable Importance Plot (if supported)
    if ("varImp" %in% methods(class = class(model))) {
      importance <- varImp(model)$importance
      importance$Feature <- rownames(importance)
      importance <- importance %>% filter(Overall > 0) %>%
                    slice_max(order_by = Overall, n = 100)
      
      importance_plot <- ggplot(importance, aes(x = reorder(Feature, Overall), y = Overall)) +
        geom_col(fill = "darkgreen") +
        coord_flip() +
        labs(title = "Variable Importance", x = "Feature", y = "Importance") +
        theme_minimal()
    }
  }
  
  # Return results and plots
  list(
    model = model,
    predictions = pred,
    confusion_matrix = confusion,
    sensitivity_specificity = sens_spec,
    mean_sensitivity = mean_sens,
    mean_specificity = mean_spec,
    confusion_plot = confusion_plot,
    importance_plot = importance_plot
  )
}

results_rf <- run_microbiome_model(
  physeq_obj = phyloseq16all,
  group_var = "cap_wild_crate",
  method = "rf",
  p = 0.8,
  cv_folds = 10,
  plot_results = TRUE)
