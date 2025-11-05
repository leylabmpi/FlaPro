# ============================================================================
# Functions for the model_notebook.ipynb notebook
# ============================================================================
# 
# 1. data handling
# 2. splitting
# 3. model training function
# 4. some vizualization functions
#
# ============================================================================


# ============================================================================
# CORE  FUNCTIONS
# ============================================================================


merge_fasta_with_df <- function(fasta_file, df, id_column = "Flagellin_ID") {
  # Read the FASTA file
  sequences <- readAAStringSet(fasta_file)
  
  # Extract sequence names from FASTA headers
  seq_names <- names(sequences)
  
  # Some FASTA headers might have additional information after a space or pipe
  # Extract just the IDs by splitting at first space or pipe
  seq_ids <- sapply(strsplit(seq_names, "[ |]"), `[`, 1)
  
  # Convert sequences to character strings
  seq_strings <- as.character(sequences)
  
  # Create a dataframe of sequence info
  seq_df <- data.frame(
    ID = seq_ids,
    sequence = seq_strings,
    stringsAsFactors = FALSE
  )
  
  # Check which IDs from the dataframe are in the FASTA file
  matching_ids <- df[[id_column]] %in% seq_df$ID
  
  if (sum(matching_ids) == 0) {
    warning("No matching IDs found between dataframe and FASTA file")
    return(df)
  }
  
  # Merge the sequences into the original dataframe
  # Create a match index
  match_idx <- match(df[[id_column]], seq_df$ID)
  
  # Add the sequences to the dataframe
  df$sequence <- NA
  df$sequence[!is.na(match_idx)] <- seq_df$sequence[match_idx[!is.na(match_idx)]]
  
  # Return the dataframe with sequences
  return(df)
}
# df_with_sequences <- merge_fasta_with_df("path.fasta", df)



split_train_test <- function(df, split = 123) {
    
    #### Split data into training and testing sets. ####
    #df : includes parameters and phenotype
    # split : random seed (def = 123) 
    #####################################################
    
    set.seed(split)
    ind <- partition(df$Exp.phenotype, p = c(train = 0.7, test = 0.3)) 
    train <- df[ind$train, ]  # training set
    test <- df[ind$test, ]  # testing set
    return(list(train, test))
}

split_train_test_cluster <- function(df, homology_clusters, split = 123) {
  set.seed(split)
  
  unique_clusters <- unique(homology_clusters)
  
  # Sample 70% of clusters for training
  train_clusters <- sample(unique_clusters, round(0.7 * length(unique_clusters)))
  
  # Create train and test sets based on clusters
  train <- df[homology_clusters %in% train_clusters, ]
  test <- df[!(homology_clusters %in% train_clusters), ]
  
  return(list(train, test))
}


create_distance_matrix_robust <- function(sequences, type = "protein", method = "dist") {
  # Remove any NA or empty sequences
  valid_idx <- which(!is.na(sequences) & sequences != "")
  sequences <- sequences[valid_idx]
  
  if (length(sequences) == 0) {
    stop("No valid sequences found")
  }
  
  # Print diagnostics
  cat("Creating distance matrix for", length(sequences), "sequences\n")
  
  # Clean sequences: remove non-standard characters
  if (type == "protein") {
    # Standard amino acid codes
    valid_chars <- c("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V")
    
    # Clean sequences
    clean_sequences <- sapply(sequences, function(seq) {
      chars <- strsplit(seq, "")[[1]]
      chars[!toupper(chars) %in% valid_chars] <- "X"  # Replace non-standard with X
      paste0(chars, collapse="")
    })
    
    sequences <- clean_sequences
  }
  
  # Convert sequences to appropriate format
  if (type == "dna") {
    seq_set <- DNAStringSet(sequences)
  } else if (type == "protein") {
    seq_set <- AAStringSet(sequences)
  }
  
  # Perform alignment
  library(msa)
  cat("Performing sequence alignment...\n")
  alignment <- msa(seq_set)
  
  if (type == "dna") {
    aligned_seqs <- DNAStringSet(as(alignment, "DNAStringSet"))
    dna_bin <- as.DNAbin(aligned_seqs)
    dist_matrix <- as.matrix(dist.dna(dna_bin, model = "raw"))
  } else {
    aligned_seqs <- AAStringSet(as(alignment, "AAStringSet"))
    
    cat("Calculating protein distance...\n")
    
    if (method == "dist") {
      # Simple Hamming distance
      aa_matrix <- as.matrix(aligned_seqs)
      n <- nrow(aa_matrix)
      dist_matrix <- matrix(0, n, n)
      
      for (i in 1:(n-1)) {
        for (j in (i+1):n) {
          pos_diff <- sum(aa_matrix[i,] != aa_matrix[j,], na.rm = TRUE)
          dist_matrix[i,j] <- dist_matrix[j,i] <- pos_diff
        }
      }
    } else if (method == "blosum") {
      library(Biostrings)
      data(BLOSUM62)
      
      # Define dist_matrix at function scope level
      dist_matrix <- NULL
      
      # Try BLOSUM62 with error handling
      success <- tryCatch({
        dist_matrix <- as.matrix(stringDist(aligned_seqs, method = "substitutionMatrix", 
                                           substitutionMatrix = BLOSUM62))
        TRUE
      }, error = function(e) {
        cat("BLOSUM calculation failed, using manual implementation\n")
        FALSE
      })
      
      # If automatic method failed, use manual implementation
      if (!success) {
        aa_matrix <- as.matrix(aligned_seqs)
        n <- nrow(aa_matrix)
        dist_matrix <- matrix(0, n, n)
        
        # Convert BLOSUM62 (similarity) to distance
        blosum_dist <- max(BLOSUM62) - BLOSUM62
        
        for (i in 1:(n-1)) {
          for (j in (i+1):n) {
            # Calculate BLOSUM-based distance
            dist_sum <- 0
            valid_pos <- 0
            
            for (pos in 1:ncol(aa_matrix)) {
              aa1 <- aa_matrix[i, pos]
              aa2 <- aa_matrix[j, pos]
              
              # Skip gaps or unknown
              if (aa1 == "-" || aa2 == "-" || aa1 == "X" || aa2 == "X") {
                next
              }
              
              # Add BLOSUM distance if both amino acids are valid
              if (aa1 %in% rownames(blosum_dist) && aa2 %in% colnames(blosum_dist)) {
                dist_sum <- dist_sum + blosum_dist[aa1, aa2]
                valid_pos <- valid_pos + 1
              }
            }
            
            # Normalize by valid positions
            if (valid_pos > 0) {
              dist_matrix[i,j] <- dist_matrix[j,i] <- dist_sum / valid_pos
            } else {
              dist_matrix[i,j] <- dist_matrix[j,i] <- max(blosum_dist)
            }
          }
        }
      }
    }
  }
  return(dist_matrix)
  #return(list(
  #  dist_matrix = dist_matrix#,
    #valid_indices = valid_idx
  #))
}


split_train_test_distance_direct <- function(df, dist_matrix, split = 123, train_ratio = 0.7) {
  set.seed(split)
  
  # Hierarchical clustering
  hc <- hclust(as.dist(dist_matrix))
  
  # Determine number of clusters
  k <- 7 #max(2, min(7, ceiling(nrow(df)/10)))
  clusters <- cutree(hc, k = k)
  
  # Assign names to clusters using rownames from df
  names(clusters) <- rownames(df)
  
  # Stratified sampling
  train_indices <- c()
  
  if ("Exp.phenotype" %in% colnames(df)) {
    for (cluster_id in unique(clusters)) {
      cluster_indices <- which(clusters == cluster_id)
      
      for (pheno in unique(df$Exp.phenotype[cluster_indices])) {
        pheno_indices <- cluster_indices[df$Exp.phenotype[cluster_indices] == pheno]
        n_sample <- max(1, round(train_ratio * length(pheno_indices)))
        sampled <- sample(pheno_indices, min(n_sample, length(pheno_indices)))
        train_indices <- c(train_indices, sampled)
      }
    }
  } else {
    for (cluster_id in unique(clusters)) {
      cluster_indices <- which(clusters == cluster_id)
      n_sample <- max(1, round(train_ratio * length(cluster_indices)))
      train_indices <- c(train_indices, sample(cluster_indices, n_sample))
    }
  }
  
  # Create train and test sets
  train <- df[train_indices, ]
  test <- df[-train_indices, ]
  
  # Add cluster assignments to the dataframe
  df$cluster <- clusters
  
  return(list(
    train = train,
    test = test,
    clusters = clusters,
    train_indices = train_indices,
    df_with_clusters = df
  ))
}

filter_low_variance <- function(df, threshold) {
                cat('* N of features initially:', ncol(df), '\n')
                # Calculate variance as before
                param_cols <- names(df)[!(names(df) %in% c("Exp.phenotype", "Cluster_CDHit", "sequence", 'pident_FlaB', 'length_FlaB', 'mismatch_FlaB', 'pident_FliC', 'length_FliC', 'mismatch_FliC', 'alntmscore_foldseek', 'rmsd_foldseek'))]
                var_values <- sapply(df[param_cols], function(x) {
                x_clean <- suppressWarnings(as.numeric(as.character(x)))
                if(all(is.na(x_clean))) return(NA)
                return(var(x_clean, na.rm = TRUE))
                })    

                var_df <- data.frame(
                Parameter = names(var_values),
                Variance = var_values
                )
                # Get parameter order based on variance
                param_order <- var_df[order(-var_df$Variance), "Parameter"]
                low_var_params <- var_df$Parameter[var_df$Variance <= threshold]

                cat('* N of features with variance <=', threshold, ":", length(low_var_params), '\n')
                cat('* Namely:', low_var_params, '\n')

                df <- df %>% select(-all_of(low_var_params))

                cat('* N of features after filtration:', ncol(df))
    return(df)
}


train_model_nestedcv <- function(train, 
                                 outer_folds_k = 10, 
                                 inner_folds_k = 10, 
                                 method = "rf", 
                                 cv_cores = 1, 
                                 seed = 123) {
                                  
  #### Train a model using nested cross-validation. ####
  # train : train data
  # outer_folds_k : outer folds. For Leave-One-Out, set outer_folds_k = length(y). (def = 10)
  # inner_folds_k : nested (inner) folds (def = 10)
  # method : method ("rf", "ranger", "svm", etc.)  (deaf = "rf")
  # cv_cores : mumber of cores for parallelization (def = 1)
  # seed : def = 123
  #####################################################

  # outcome variable is a factor
  train$Exp.phenotype <- as.factor(train$Exp.phenotype)
  
  # x = features, y = response (last column is the phenotype)
  y <- train$Exp.phenotype
  x <- train[, -ncol(train)]
  
  # STEP 1: Create OUTER folds (for the outer loop)
  out_folds <- caret::createFolds(y, k = outer_folds_k)
  
  # STEP 2: Create INNER folds (for the inner loop)
  # For each outer split, create new folds on the training subset
  in_folds <- lapply(out_folds, function(idx) {
    # The training portion in the outer loop
    train_y <- y[-idx]
    caret::createFolds(train_y, k = inner_folds_k)
  })
  
  # STEP 3: Run NESTED CV (using the nestedcv package)
  # Here, pass the folds to nestcv.train

  set.seed(seed)
  
  res <- nestcv.train(y, x, 
                      method      = method,    
                      cv.cores    = cv_cores,  
                      inner_folds = in_folds,  
                      outer_folds = out_folds) 
  
  return(res)
}


evaluate_v3_robust <- function(df, iterations = 3) {
  library(caret)
  library(dplyr)
  
  # Initialize results data frame
  results <- data.frame(
    Iteration = integer(),
    Method = character(),
    TestAccuracy = numeric(),
    TrainAccuracy = numeric(),
    Sensitivity = numeric(),
    Specificity = numeric(),
    PPV = numeric(),
    NPV = numeric(),
    F1 = numeric(),
    Kappa = numeric(),
    BalancedAccuracy = numeric(),
    TrainSize = integer(),
    TestSize = integer(),
    TP = integer(),
    FP = integer(),
    TN = integer(),
    FN = integer(),
    positive_class = character(),
    negative_class = character()
  )
  
  # Process iterations sequentially for reliability
  for (i in 1:iterations) {
    cat("Processing iteration", i, "of", iterations, "\n")
    
    tryCatch({
      # 2. Random split
      set.seed(i)
      random_indices <- sample(1:nrow(df), round(0.7 * nrow(df)))
      random_train <- df[random_indices,]
      random_test <- df[-random_indices,]
      
      random_train_shaffle <- random_train
      random_train_shaffle[["Exp.phenotype"]] <- sample(random_train_shaffle[["Exp.phenotype"]])
      
      # Define datasets
      datasets <- list(
        random_split = list(train = random_train, test = random_test),
        random_cluster_shaffle = list(train = random_train_shaffle, test = random_test)
      )
      
      # Process each dataset
      for (method_name in names(datasets)) {
        cat("  Method:", method_name, "\n")
        train_data <- datasets[[method_name]]$train
        test_data <- datasets[[method_name]]$test
        
        # Remove columns consistently across all split types
        cols_to_remove <- c("Cluster_CDHit", "cluster", "sequence")
        for (col in cols_to_remove) {
          if (col %in% colnames(train_data)) {
            train_data <- train_data %>% select(-all_of(col))
            if (col %in% colnames(test_data)) {
              test_data <- test_data %>% select(-all_of(col))
            }
          }
        }
        
        # Ensure consistent variables
        common_vars <- intersect(names(train_data), names(test_data))
        train_data <- train_data[, common_vars]
        test_data <- test_data[, common_vars]
        
        # Ensure factor levels match
        for (col in names(train_data)) {
          if (is.factor(train_data[[col]])) {
            all_levels <- union(levels(train_data[[col]]), levels(test_data[[col]]))
            train_data[[col]] <- factor(train_data[[col]], levels = all_levels)
            test_data[[col]] <- factor(test_data[[col]], levels = all_levels)
          }
        }
        
        # Use memory-efficient garbage collection
        gc()
        
        # Train model with error handling
        model <- tryCatch({
          train_model_nestedcv(train_data, seed = i)
        }, error = function(e) {
          cat("  Error training model:", e$message, "\n")
          NULL
        })
        
        if (is.null(model)) next
        
        # Predict on test set
        test_preds <- predict(model, newdata = test_data)
        test_accuracy <- sum(test_preds == test_data$Exp.phenotype) / nrow(test_data)
        
        # Predict on training set
        train_preds <- predict(model, newdata = train_data)
        train_accuracy <- sum(train_preds == train_data$Exp.phenotype) / nrow(train_data)
        
        # Calculate confusion matrix with error handling
        tryCatch({
          # Ensure factors have same levels for confusion matrix
          levels_to_use <- union(levels(test_data$Exp.phenotype), unique(test_preds))
          test_preds_factor <- factor(test_preds, levels = levels_to_use)
          test_phenotype_factor <- factor(test_data$Exp.phenotype, levels = levels_to_use)
          
          cm <- confusionMatrix(test_preds_factor, test_phenotype_factor)
          
          # Get the positive/negative classes
          positive_class <- levels(test_preds_factor)[2]
          negative_class <- levels(test_preds_factor)[1]
          
          # Extract statistics
          sensitivity <- cm$byClass["Sensitivity"]
          specificity <- cm$byClass["Specificity"]
          ppv <- cm$byClass["Pos Pred Value"]
          npv <- cm$byClass["Neg Pred Value"]
          f1 <- cm$byClass["F1"]
          kappa <- cm$overall["Kappa"]
          balanced_acc <- cm$byClass["Balanced Accuracy"]
          
          # Get confusion matrix values
          cm_table <- cm$table
          if(nrow(cm_table) == 2 && ncol(cm_table) == 2) {
            TP <- cm_table[2, 2]  # True positive
            FP <- cm_table[2, 1]  # False positive
            TN <- cm_table[1, 1]  # True negative
            FN <- cm_table[1, 2]  # False negative
          } else {
            TP <- FP <- TN <- FN <- NA
          }
        }, error = function(e) {
          cat("  Error in confusion matrix calculation:", e$message, "\n")
          sensitivity <- specificity <- ppv <- npv <- f1 <- kappa <- balanced_acc <- NA
          TP <- FP <- TN <- FN <- NA
          positive_class <- negative_class <- NA_character_
        })
        
        # Add results for this method
        results <- rbind(results, 
                        data.frame(
                          Iteration = i, 
                          Method = method_name, 
                          TestAccuracy = test_accuracy,
                          TrainAccuracy = train_accuracy,
                          Sensitivity = sensitivity,
                          Specificity = specificity,
                          PPV = ppv,
                          NPV = npv,
                          F1 = f1,
                          Kappa = kappa,
                          BalancedAccuracy = balanced_acc,
                          TrainSize = nrow(train_data),
                          TestSize = nrow(test_data),
                          TP = TP,
                          FP = FP,
                          TN = TN,
                          FN = FN,
                          positive_class = positive_class,
                          negative_class = negative_class
                        ))
        
        # Save intermediate results after each method
        #saveRDS(results, file = paste0("results_intermediate_", i, "_", method_name, ".rds"))
      }
    }, error = function(e) {
      cat("Error in iteration", i, ":", e$message, "\n")
    })
    
    # Save results after each iteration
    #saveRDS(results, file = paste0("results_iteration_", i, ".rds"))
  }
  
  return(results)
}

# plot evaluation results

plot_evaluation_results <- function(results, class_labels = c("Negative", "Positive")) {
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(ggforce)
  library(gridExtra)
  #library(ggpattern)  # Add ggpattern for pattern fills
  
  # Extract method prefix for coloring and order methods
  results <- results %>%
    mutate(
      Method_Prefix = gsub("_.*", "", Method),  # Extract prefix (text before first underscore)
      Has_Shaffle = grepl("_shaffle", Method),   # Flag if method contains "_shaffle"
      # Create a method base for ordering (remove "_shaffle" if present)
      Method_Base = gsub("_shaffle", "", Method)
    )
    
  # Special handling for ordering - create a custom sort function
  # This ensures consistent ordering across all method types
  custom_sort <- function(df) {
    # Group and sort the methods
    groups <- split(df, df$Method_Prefix)
    sorted_df <- data.frame()
    
    for (prefix in names(groups)) {
      group_data <- groups[[prefix]]
      # Get unique base methods
      base_methods <- unique(group_data$Method_Base)
      
      for (base in base_methods) {
        # Get rows for this base method
        base_rows <- group_data[group_data$Method_Base == base,]
        # Sort by shuffle status (non-shuffled first)
        base_rows <- base_rows[order(base_rows$Has_Shaffle),]
        # Add to result
        sorted_df <- rbind(sorted_df, base_rows)
      }
    }
    return(sorted_df)
  }
  
  # Apply custom sorting
  sorted_methods <- custom_sort(distinct(results, Method, Method_Prefix, Method_Base, Has_Shaffle))
  
  # Apply the factor with the sorted levels
  results$Method_Factor <- factor(results$Method, levels = sorted_methods$Method)
  
  # Summary statistics
  summary_stats <- results %>%
    group_by(Method, Method_Prefix, Has_Shaffle, Method_Factor) %>%
    summarize(across(c(TestAccuracy, TrainAccuracy, Sensitivity, 
                     Specificity, PPV, NPV, F1, Kappa, BalancedAccuracy),
                    list(mean = mean, sd = sd), na.rm = TRUE),
              .groups = "drop")
  
  # Reshape for visualization
  metrics_long <- results %>%
    select(Method, Method_Prefix, Has_Shaffle, Method_Factor, Iteration, TestAccuracy, Sensitivity, Specificity, 
           PPV, NPV, F1, BalancedAccuracy) %>%
    pivot_longer(cols = c(TestAccuracy, Sensitivity, Specificity, 
                         PPV, NPV, F1, BalancedAccuracy),
                names_to = "Metric", values_to = "Value")
  
  # Create train vs test accuracy plot
  acc_data <- results %>%
    select(Method, Method_Prefix, Has_Shaffle, Method_Factor, Iteration, TestAccuracy, TrainAccuracy) %>%
    pivot_longer(cols = c(TestAccuracy, TrainAccuracy),
                names_to = "AccuracyType", values_to = "Value")
  
  
  # Check if confusion matrix data is available
  has_cm_data <- all(c("TP", "FP", "TN", "FN") %in% colnames(results))
  
  if (has_cm_data) {
    # Calculate mean and sd for each confusion matrix cell by method
    cm_summary <- results %>%
      group_by(Method, Method_Prefix, Has_Shaffle, Method_Factor) %>%
      summarize(
        TP_mean = mean(TP, na.rm = TRUE),
        TP_sd = sd(TP, na.rm = TRUE),
        FP_mean = mean(FP, na.rm = TRUE),
        FP_sd = sd(FP, na.rm = TRUE),
        TN_mean = mean(TN, na.rm = TRUE),
        TN_sd = sd(TN, na.rm = TRUE),
        FN_mean = mean(FN, na.rm = TRUE),
        FN_sd = sd(FN, na.rm = TRUE),
        .groups = "drop"
      )
    
    # Check if we have class label information in the results
    has_class_labels <- all(c("positive_class", "negative_class") %in% colnames(results))
    
    # Get the class labels
    neg_label <- class_labels[1]
    pos_label <- class_labels[2]
    
    # If we have class labels in the results, use the most common ones
    if (has_class_labels) {
      # Filter out NA values
      valid_labels <- results %>% 
        filter(!is.na(positive_class) & !is.na(negative_class))
      
      if (nrow(valid_labels) > 0) {
        # Get the most frequent positive and negative class labels
        pos_counts <- table(valid_labels$positive_class)
        neg_counts <- table(valid_labels$negative_class)
        
        if (length(pos_counts) > 0 && length(neg_counts) > 0) {
          pos_label <- names(pos_counts)[which.max(pos_counts)]
          neg_label <- names(neg_counts)[which.max(neg_counts)]
        }
      }
    }
    
    # Function to create a confusion matrix plot for a specific method
    create_cm_plot <- function(method_data, method_name) {
      # Create data for the plot
      plot_data <- data.frame(
        Position = c("TP", "FP", "FN", "TN"),
        Mean = c(method_data$TP_mean, method_data$FP_mean, 
                 method_data$FN_mean, method_data$TN_mean),
        SD = c(method_data$TP_sd, method_data$FP_sd, 
               method_data$FN_sd, method_data$TN_sd)
      )
      
      # Create the text to display in each cell
      plot_data$Label <- sprintf("%.1f ± %.1f", plot_data$Mean, plot_data$SD)
      
      # Position the cells in a confusion matrix layout
      # Standard confusion matrix layout:
      # [TN, FP]
      # [FN, TP]
      plot_data$row <- c(2, 1, 2, 1)  # Row positions (TP is bottom right)
      plot_data$col <- c(2, 2, 1, 1)  # Column positions (TN is top left)
      
      # Add value column for color intensity based on mean value
      plot_data$Value <- plot_data$Mean
      
      # Create the plot with blue gradient
      p <- ggplot(plot_data, aes(x = factor(col), y = factor(row))) +
        geom_tile(aes(fill = Value), color = "black", size = 0.5) +
        scale_fill_gradient(low = "#e6f3ff", high = "#084b82") +  
        geom_text(aes(label = Label), size = 4, color = "black") +
        labs(title = paste("Confusion Matrix for", method_name),
             x = "Predicted Class", y = "Actual Class") +
        theme_minimal() +
        theme(legend.position = "none",
              axis.text.x = element_text(size = 12),
              axis.text.y = element_text(size = 12),
              plot.title = element_text(hjust = 0.5, size = 14)) +
        scale_x_discrete(labels = c(neg_label, pos_label)) +  # Predicted class
        scale_y_discrete(labels = c(neg_label, pos_label))    # Actual class
      
      return(p)
    }
    
    # Create confusion matrix plots for each method
    methods <- unique(results$Method)
    cm_plots <- list()
    
    for (method in methods) {
      method_data <- cm_summary %>% filter(Method == method)
      if (nrow(method_data) > 0) {
        cm_plots[[method]] <- create_cm_plot(method_data, method)
      }
    }
    
    # Calculate how many plots per row based on the number of methods
    plots_per_row <- min(2, length(methods))
    
    # Combine all confusion matrix plots
    combined_cm_plot <- do.call(grid.arrange, c(cm_plots, ncol = plots_per_row))
    
    # Create a summary table of confusion matrix values
    cm_table <- cm_summary %>%
      mutate(
        TP = sprintf("%.1f ± %.1f", TP_mean, TP_sd),
        FP = sprintf("%.1f ± %.1f", FP_mean, FP_sd),
        TN = sprintf("%.1f ± %.1f", TN_mean, TN_sd),
        FN = sprintf("%.1f ± %.1f", FN_mean, FN_sd)
      ) %>%
      select(Method, TP, FP, TN, FN)
    
    # Return all results including confusion matrix
    return(list(
      #metrics_plot = p1, 
      #accuracy_plot = p2, 
      summary = summary_stats,
      confusion_matrix_plots = cm_plots,
      combined_confusion_matrix = combined_cm_plot,
      confusion_matrix_summary = cm_table
    ))
  } else {
    # Return original results if no confusion matrix data
    return(list(metrics_plot = p1, accuracy_plot = p2, summary = summary_stats))
  }
}
