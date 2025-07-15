# Making a good classification model 

#### Load libraries ####
pacman::p_load(tidyverse, vegan, MASS, phyloseq, tibble, ANCOMBC, ggplot2, coin,
               MicrobiotaProcess, patchwork, reshape2, ggnewscale, VennDiagram,
               UpSetR, gridExtra, grid, WGCNA, gt, caret, randomForest, pROC, e1071)
library(calecopal)
bigsur = c("#E4DECE", "#ECBD95", "#9BB1BB", "#79ACBD", "#346575", "#0B4221")

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
ps <- phyloseq_no_new_crates

#### MANUAL MICROBIOME PREDICTION MODEL ####

# largely based on this: https://microbiome.netlify.app/classification-using-microbiome.html

run_microbiome_model <- function(physeq_obj,
                                   group_var,
                                   method = "rf",
                                   p = 0.75,
                                   cv_folds = 10,
                                   seed = 123,
                                   plot_results = TRUE) {
    library(phyloseq)
    library(caret)
    library(ggplot2)
    library(dplyr)
    
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
    
    # Combine OTU data and group variable
    data <- cbind(otu_t, Group = meta[[group_var]])
    data$Group <- as.factor(make.names(data$Group))  # ensure valid R names
    
    # Train/test split
    trainIndex <- createDataPartition(data$Group, p = p, list = FALSE)
    trainData <- data[trainIndex, ]
    testData  <- data[-trainIndex, ]
    
    # Cross-validation control
    ctrl <- trainControl(method = "cv",
                         number = cv_folds,
                         classProbs = TRUE,
                         summaryFunction = multiClassSummary,
                         savePredictions = TRUE)
    
    # Train model
    model <- train(Group ~ ., data = trainData,
                   method = method,
                   trControl = ctrl,
                   metric = "Accuracy")
    
    # Predict on test data
    pred <- predict(model, newdata = testData)
    
    # Confusion matrix
    confusion <- confusionMatrix(pred, as.factor(testData$Group))
    
    # Sensitivity and specificity
    sens_spec  <- confusion$byClass[, c("Sensitivity", "Specificity")]
    mean_sens  <- mean(confusion$byClass[, "Sensitivity"], na.rm = TRUE)
    mean_spec  <- mean(confusion$byClass[, "Specificity"], na.rm = TRUE)
    
    # Initialize plots
    confusion_plot  <- NULL
    importance_plot <- NULL
    
    if (plot_results) {
      ## Confusion matrix heatmap
      cm_df <- as.data.frame(confusion$table)
      colnames(cm_df) <- c("Prediction", "Reference", "Freq")
      
      confusion_plot <- ggplot(cm_df, aes(x = Reference, y = Prediction, fill = Freq)) +
        geom_tile() +
        geom_text(aes(label = Freq), color = "white", size = 5) +
        scale_fill_gradient(low = "steelblue", high = "#FF4D4D") +
        labs(title = "Confusion Matrix Heatmap", x = "Actual", y = "Predicted") +
        theme_minimal()
      
      ## Variable importance plot
      importance_result <- try(varImp(model), silent = TRUE)
      
      if (!inherits(importance_result, "try-error") && "importance" %in% names(importance_result)) {
        importance <- importance_result$importance
        
        if ("Overall" %in% colnames(importance)) {
          importance$Feature <- rownames(importance)
          importance <- importance %>%
            filter(Overall > 0) %>%
            slice_max(order_by = Overall, n = 60)
          
          importance_plot <- ggplot(importance, aes(x = reorder(Feature, Overall), y = Overall)) +
            geom_col(fill = "darkgreen") +
            coord_flip() +
            labs(title = "Top 60 Variable Importance",
                 x = "Feature", y = "Importance") +
            theme_minimal()
        }
      }
    }
    
    # Return everything as a named list
    return(list(
      model                = model,
      predictions          = pred,
      confusion_matrix     = confusion,
      sensitivity_specificity = sens_spec,
      mean_sensitivity     = mean_sens,
      mean_specificity     = mean_spec,
      confusion_plot       = confusion_plot,
      importance_plot      = importance_plot
    ))
  }
  
# Options for model inputs : 
model_options <- modelLookup() %>% filter(forClass == TRUE)
  
results_rf <- run_microbiome_model(
  physeq_obj = ps,
  group_var = "cap_wild_crate",
  method = "rf",
  p = 0.75, seed = 241, cv_folds = 10,
  plot_results = TRUE)


results_pcaNNet <- run_microbiome_model(
  physeq_obj = ps,
  group_var = "cap_wild_crate",
  method = "pcaNNet",
  p = 0.75, seed = 241, cv_folds = 10,
  plot_results = TRUE)

results_xgbTree <- run_microbiome_model(
  physeq_obj = ps,
  group_var = "cap_wild_crate",
  method = "xgbTree",
  p = 0.75, seed = 241, cv_folds = 10,
  plot_results = TRUE)

save(results_rf, 
     results_pcaNNet,
     results_xgbTree,
     file = "06_grp_prediction_models/rf_pcannet_xgbtre_no_crate.RData")

ggsave(filename = "06_grp_prediction_models/rf_confusion_plot_no_crates.png",
       plot = results_rf$confusion_plot, width = 8, height = 6, dpi = 300)

ggsave(filename = "06_grp_prediction_models/rf_importance_plot_no_crates.png",
       plot = results_rf$importance_plot, width = 8, height = 6, dpi = 300)

# Save pcaNNet plots
ggsave(filename = "06_grp_prediction_models/pcannet_confusion_plot_no_crates.png",
       plot = results_pcaNNet$confusion_plot, width = 8, height = 6, dpi = 300)

ggsave(filename = "06_grp_prediction_models/pcannet_importance_plot_no_crates.png",
       plot = results_pcaNNet$importance_plot, width = 8, height = 6, dpi = 300)

# Save xgbTree plots
ggsave(filename = "06_grp_prediction_models/xgbtree_confusion_plot_no_crates.png",
       plot = results_xgbTree$confusion_plot, width = 8, height = 6, dpi = 300)

ggsave(filename =  "06_grp_prediction_models/xgbtree_importance_plot_no_crates.png",
       plot = results_xgbTree$importance_plot, width = 8, height = 6, dpi = 300)
  
### Plot and compare models 

# Step 1: Extract and format the metrics from each model
extract_metrics <- function(model_result, model_name) {
  metrics <- model_result$confusion_matrix$byClass[, c("Sensitivity", "Specificity", "Balanced Accuracy")]
  metrics_df <- as.data.frame(metrics)
  metrics_df$Class <- rownames(metrics_df)
  metrics_df$Model <- model_name
  return(metrics_df)
}

# Extract metrics from all models
rf_df       <- extract_metrics(results_rf, "Random Forest")
pcaNNet_df  <- extract_metrics(results_pcaNNet, "pcaNNet")
xgbTree_df  <- extract_metrics(results_xgbTree, "xgbTree")

# Combine all into df and make long
all_metrics <- bind_rows(rf_df, pcaNNet_df, xgbTree_df)
metrics_long <- all_metrics %>%
  pivot_longer(cols = c(Sensitivity, Specificity, `Balanced Accuracy`),
               names_to = "Metric", values_to = "Value")

# Clean up class labels
metrics_long$Class <- gsub("\\.\\.", " ", gsub("\\.", " ", metrics_long$Class))

mod_results_compare <- ggplot(metrics_long, aes(x = Class, y = Value, fill = Model)) +
  geom_col(position = position_dodge()) +
  facet_wrap(~Metric, scales = "free_y") +
  theme_minimal(base_size = 14, base_family = "") +  # Increase base font size
  labs(title = "Model Performances by Parrot Group",
    y = "Metric Value", x = "Group") +
  scale_fill_manual(values = cal_palette("superbloom3")) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),  # White background
    panel.background = element_rect(fill = "white", color = NA), 
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(face = "bold", size = 16),      
    axis.text = element_text(face = "bold", size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_text(face = "bold", size = 16),
    legend.title = element_text(face = "bold", size = 14),
    legend.text = element_text(face = "bold", size = 12),
    plot.title = element_text(face = "bold", size = 18, hjust = 0.5))

ggsave(filename = "06_grp_prediction_models/compare_models_no_crates.png",mod_results_compare, height=15, width=18)


### Testing out network embeddings 

data(enterotype)
enterotype = subset_samples(ps, !is.na(ps))

plot_net(enterotype, maxdist = 0.4, point_label = "cap_wild_crate")
plot_net(enterotype, maxdist = 0.3, color = "cap_wild_crate", shape="Country")
ig <- make_network(enterotype, max.dist=0.7)
plot_network(ig, enterotype, color="cap_wild_crate", 
             shape="Country", line_weight=0.4, label=NULL)

ig <- make_network(enterotype, dist.fun="bray", max.dist=0.7)
plot_network(ig, enterotype, color="Site.name", 
             shape="cap_wild_crate", line_weight=0.5, label=NULL)

# Trying another

library(rSETSe)

g <- generate_peels_network(type = "E") %>%
     prepare_edges(k = 500, distance = 1) %>%
     prepare_categorical_force(., node_names = "name", force_var = "class") 

#Embedds using the bi-connected auto-parametrization algorithm.
#This method is strongly reccomended, it tends to be much faster and almost always converges
embeddings <- setse_bicomp(g,
                           force = "class_A",
                           tol = sum(abs(vertex_attr(g, "class_A")))/1000,
                           hyper_tol = 0.1,
                           hyper_iters = 3000,
                           verbose = T)



