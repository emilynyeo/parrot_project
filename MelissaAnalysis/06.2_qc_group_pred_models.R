# Making a good classification model 

#### Load libraries ####
pacman::p_load(tidyverse, vegan, MASS, phyloseq, tibble, ANCOMBC, ggplot2, coin,
               MicrobiotaProcess, patchwork, reshape2, ggnewscale, VennDiagram,
               UpSetR, gridExtra, grid, WGCNA, gt, caret, randomForest, pROC, 
               e1071, calecopal, purrr, stringr, dplyr, ggtree)

bigsur = c("#E4DECE", "#ECBD95", "#9BB1BB", "#79ACBD", "#346575", "#0B4221")

#BiocManager::install(c("GO.db", "impute", "preprocessCore"))

#### Load data ####
run_microbiome_model <- function(physeq_obj,
                                 group_var,
                                 method = "rf",
                                 p = 0.75,
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


# Step 1: Load all phyloseq RDS files
phyloseq_files <- list.files("01.1_qc_checks", pattern = "^phyloseq_.*\\.rds$", full.names = TRUE)

file_df <- data.frame(
  path = phyloseq_files,
  file = basename(phyloseq_files),
  stringsAsFactors = FALSE) %>%
  mutate(
    clean_id = sub("^phyloseq_", "", tools::file_path_sans_ext(file)),
    base = str_extract(clean_id, "^otu\\d{2}[a-z]+(_nohost)?_mat_filt"),
    suffix_thresh = str_replace(clean_id, "^otu\\d{2}[a-z]+(_nohost)?_mat_filt_?", ""),
    clean_id = str_remove(clean_id, "_mat_filt"),
    base = str_remove(base, "_mat_filt"))

# Step 3: Define grouping for merge (e.g., 16n + 16ren)

merge_groups <- list(
  "phyloseq16n" = c("otu16n", "otu16ren"),
  "phyloseq18n" = c("otu18n_nohost", "otu18ren_nohost"))

non_merge_groups <- list(
  "phyloseq16m" = c("otu16m"),
  "phyloseq18m" = c("otu18m_nohost"))

all_groups <- list(
  merge = merge_groups,
  no_merge = non_merge_groups)

# Step 4: Loop over suffix/threshold combos to merge and model
#suffixes_thresholds <- unique(file_df$suffix_thresh)
# suffixes_thresholds <- unique(file_df$suffix_thresh)
for (sfx in suffixes_thresholds) {
  for (merge_status in names(all_groups)) {
    group_list <- all_groups[[merge_status]]
    
    for (group_name in names(group_list)) {
      base_parts <- group_list[[group_name]]
      
      # Get matching files
      matching_paths <- file_df %>%
        filter(base %in% base_parts, suffix_thresh == sfx) %>%
        arrange(base)
      
      # Skip if necessary parts not found
      if (merge_status == "merge" && nrow(matching_paths) < 2) {
        message("Skipping ", group_name, "_", sfx, " (not all parts found)")
        next
      }
      if (merge_status == "no_merge" && nrow(matching_paths) < 1) {
        message("Skipping ", group_name, "_", sfx, " (file not found)")
        next
      }
      
      # Load phyloseq objects
      phylo_objs <- lapply(matching_paths$path, function(path) {
        env <- new.env()
        obj_name <- load(path, envir = env)
        env[[obj_name]]
      })
      
      # Merge if needed, or use single object
      ps_merged <- if (merge_status == "merge") {
        do.call(merge_phyloseq, phylo_objs)
      } else {
        phylo_objs[[1]]
      }
      
      clean_id <- paste0(group_name, "_", sfx)
      
#for (sfx in suffixes_thresholds) {
#  for (group_name in names(merge_groups)) {
#    base_parts <- merge_groups[[group_name]]
    
    # Get matching files
#    matching_paths <- file_df %>%
#      filter(base %in% base_parts, suffix_thresh == sfx) %>%
#      arrange(base)  # consistent order
    
#    if (nrow(matching_paths) < 2) {
#      message("Skipping ", group_name, "_", sfx, " (not all parts found)")
#      next
#    }
    
    # Load and merge
#    phylo_objs <- lapply(matching_paths$path, function(path) {
#      env <- new.env()
#      obj_name <- load(path, envir = env)
#      env[[obj_name]]
#    })
    
#    ps_merged <- do.call(merge_phyloseq, phylo_objs)
    
    # Clean up group label
#    clean_id <- paste0(group_name, "_", sfx)
    
    # ** Combine Crate Data into one **
    sample_df <- as.data.frame(sample_data(ps_merged))
    
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
    sample_data(ps_merged) <- sample_data(sample_df)
    levels(sample_data(ps_merged)$cap_wild_crate)
    
    # Check for required variable
    if (!"cap_wild_crate" %in% colnames(sample_data(ps_merged))) {
      message("Skipping ", clean_id, " - missing 'cap_wild_crate'")
      next
    }
    
    # Create new phyloseq object without "New crates"
    phyloseq_no_new_crates <- subset_samples(ps_merged, cap_wild_crate != "New crates")
    levels(sample_data(phyloseq_no_new_crates)$cap_wild_crate)
    
    # Run models
    cat("Running models for:", clean_id, "\n")
    
    # Options for model inputs : 
    model_options <- modelLookup() %>% filter(forClass == TRUE)
    
    results_rf <- run_microbiome_model(
      physeq_obj = phyloseq_no_new_crates,
      group_var = "cap_wild_crate",
      method = "rf",
      p = 0.75, seed = 241, cv_folds = 4,
      plot_results = TRUE)
    
    #results_pcaNNet <- run_microbiome_model(
    #  physeq_obj = phyloseq_no_new_crates,
    #  group_var = "cap_wild_crate",
    #  method = "pcaNNet",
    #  p = 0.75, seed = 241, cv_folds = 4,
    #  plot_results = TRUE)
    
    results_xgbTree <- run_microbiome_model(
      physeq_obj = phyloseq_no_new_crates,
      group_var = "cap_wild_crate",
      method = "xgbTree",
      p = 0.75, seed = 241, cv_folds = 4,
      plot_results = TRUE)
    
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
    metrics_long$Class <- gsub("\\.\\.", " ", 
                               gsub("\\.", " ", metrics_long$Class)) # clean label
    
    write.csv(metrics_long,
              file = paste0("06.2_qc_group_pred_models/compare_mods_no_crates_metrics", clean_id,".csv"),
              row.names = FALSE)
    
    mod_results_compare <- ggplot(metrics_long, aes(x = Class, y = Value, fill = Model)) +
      geom_col(position = position_dodge()) +
      facet_wrap(~Metric, scales = "free_y") +
      theme_minimal(base_size = 14, base_family = "") +  # Increase base font size
      labs(title = paste0("Model Performances by Parrot Group for QC: ", clean_id),
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
    
    ggsave(filename = paste0("06.2_qc_group_pred_models/compare_models_no_crates", clean_id,".png"),
           mod_results_compare, height=15, width=18)
    
    
    # Save models
    save(results_rf, results_pcaNNet, results_xgbTree,
         file = paste0("06.2_qc_group_pred_models/models_rf_pcannet_xgbtree", clean_id, ".RData"))
    
    # Save RF Plots 
    ggsave(filename = paste0("06.2_qc_group_pred_models/rf_confusion_plot_", clean_id, ".png"),
           plot = results_rf$confusion_plot, width = 8, height = 6, dpi = 300)
    
    ggsave(filename = paste0("06.2_qc_group_pred_models/rf_importance_plot_", clean_id, ".png"),
           plot = results_rf$importance_plot, width = 8, height = 6, dpi = 300)
    
    # Save pcaNNet plots
    ggsave(filename = paste0("06.2_qc_group_pred_models/pcannet_confusion_plot_", clean_id, ".png"),
           plot = results_pcaNNet$confusion_plot, width = 8, height = 6, dpi = 300)
    
    ggsave(filename = paste0("06.2_qc_group_pred_models/pcannet_importance_plot_", clean_id, ".png"),
           plot = results_pcaNNet$importance_plot, width = 8, height = 6, dpi = 300)
    
    # Save xgbTree plots
    ggsave(filename = paste0("06.2_qc_group_pred_models/xgbtree_confusion_plot_", clean_id, ".png"),
           plot = results_xgbTree$confusion_plot, width = 8, height = 6, dpi = 300)
    
    ggsave(filename = paste0("06.2_qc_group_pred_models/xgbtree_importance_plot_", clean_id, ".png"),
           plot = results_xgbTree$importance_plot, width = 8, height = 6, dpi = 300)
    
    # Save plots
    plots <- list(
      rf_conf = results_rf$confusion_plot,
      rf_imp = results_rf$importance_plot,
      pca_conf = results_pcaNNet$confusion_plot,
      pca_imp = results_pcaNNet$importance_plot,
      xgb_conf = results_xgbTree$confusion_plot,
      xgb_imp = results_xgbTree$importance_plot)
    
    for (plot_type in names(plots)) {
      p <- plots[[plot_type]]
      if (!is.null(p)) {
        ggsave(
          filename = paste0("06.2_qc_group_pred_models/", plot_type, "_", clean_id, ".png"),
          plot = p, width = 8, height = 6, dpi = 300
        )
      }
    }
  }
}}


## Compare balanced accuracy of all models 

model_files <- list.files("06.2_qc_group_pred_models", 
                          pattern = "^models_rf_pcannet_xgbtree.*\\.RData$", 
                          full.names = TRUE)

extract_metrics <- function(model_result, model_name, id) {
  metrics <- model_result$confusion_matrix$byClass[, c("Sensitivity", "Specificity", "Balanced Accuracy")]
  metrics_df <- as.data.frame(metrics)
  metrics_df$Class <- rownames(metrics_df)
  metrics_df$Model <- model_name
  metrics_df$ID <- id
  return(metrics_df)
}
  
all_metrics_list <- list()
for (file in model_files) {
    # Extract clean_id from filename
    clean_id <- sub("models_rf_pcannet_xgbtree", "", 
                    tools::file_path_sans_ext(basename(file)))
    
    # Load models into environment
    env <- new.env()
    load(file, envir = env)
    
    # Check if all models are present
    models_present <- c("results_rf", "results_pcaNNet", "results_xgbTree") %in% ls(env)
    
    if (!all(models_present)) {
      warning(paste("Missing model in file:", file))
      next
    }
  
    # Extract metrics
    rf_df <- extract_metrics(env$results_rf, "Random Forest", clean_id)
    xgb_df <- extract_metrics(env$results_xgbTree, "xgbTree", clean_id)
    pca_df <- extract_metrics(env$results_pcaNNet, "pcaNNet", clean_id)
    
    # Combine and store
    all_metrics_list[[clean_id]] <- bind_rows(rf_df, xgb_df, pca_df)
}

all_metrics_df <- bind_rows(all_metrics_list)

all_metrics_df <- all_metrics_df %>% 
  mutate(ID = str_replace_all(ID, "_", " "),              #underscores to spaces
         #ID = str_replace(ID, "(.* )", "r\\1"),           #add "r" before last group
         ID = str_replace(ID, "thr", "p")) 

balanced_df <- all_metrics_df %>%
  select(ID, Model, Class, `Balanced Accuracy`) %>% 
  filter(Model != "pcaNNet")

all_plot <- ggplot(balanced_df, aes(x = ID, y = `Balanced Accuracy`, fill = Class)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = cal_palette("superbloom3")) +
  facet_wrap(~ Model, scales = "free_x") +
  theme_minimal(base_size = 14) +
  labs(title = "Balanced Accuracy Comparison (from saved models)",
       y = "Balanced Accuracy", x = "Group") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(filename = paste0("06.2_qc_group_pred_models/compare_RF_XGBTREE_models_no_crates.png"),
       all_plot, height=15, width=18)

