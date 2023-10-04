options(scipen = 999)

library(GEOquery)
library(preprocessCore)
library(sva)
library(randomForest)  # Add the randomForest package
library(pheatmap)
library(ggplot2)
library(limma)
library(AnnotationDbi)
library(clusterProfiler)
library(DOSE)
library(glmnet)
library(caret)
library(dplyr)
library(caret)

# Get DataSet from GEO Database
gse <- getGEO("GSE66494", GSEMatrix = TRUE)

# Get the Labels (possible future use)
gseData <- pData(gse[[1]])

# Import Dataset and Save
Expr_GSM1623299 <- exprs(gse[[1]])

# Create the labels:
statuses <- c(
  rep(1, 47),  # Original CKD patients
  rep(0, 5),   # Original controls
  rep(1, 5),   # Following 5 CKD patients
  rep(0, 3)    # Last 3 original controls
)

# Number of iterations
num_iterations <- 2

#Model accuracy
accuracy_list <- numeric(num_iterations)

# Create a data frame to store the importance of each gene for each seed
gene_importance_df <- data.frame(Gene=rownames(important_genes_sorted))


for(i in 1:num_iterations) {
  
  # Set the seed for shuffling the data
  set.seed(127 + i)  # Each iteration will have a different seed
  
  control_indices <- which(statuses == 0)
  control_shuffled <- sample(control_indices)
  
  # Shuffle the CKD-positive patients:
  ckd_indices <- which(statuses == 1)
  ckd_shuffled <- sample(ckd_indices)
  
  # Update the indices:
  first_40_ckd <- ckd_shuffled[1:40]
  next_4_controls <- control_shuffled[1:4]
  next_12_ckd_test <- ckd_shuffled[41:52]
  last_4_controls_test <- control_shuffled[5:8]
  
  random_order_total <- c(first_40_ckd, next_4_controls, next_12_ckd_test, last_4_controls_test)
  
  # Split the shuffled labels into training and testing sets:
  statustrain <- statuses[random_order_total[1:44]]
  statustest <- statuses[random_order_total[45:60]]
  
  # Apply this shuffling to the adjusted dataset:
  combineddataset_shuffled <- combineddataset[, random_order_total, drop = FALSE]
  Train <- combineddataset_shuffled[, 1:44]
  Test <- combineddataset_shuffled[, 45:60]
  
  # Update the columns for the Train and Test datasets:
  Train_cols <- colnames(Train)
  Test_cols <- colnames(Test)
  
  Train <- as.matrix(Train)
  Test <- as.matrix(Test)
  
  #----------------------------------------------------
  
  # Fetch the GEO dataset by accession number
  geo_data <- getGEO("GPL6480")
  
  # Access the expression data (data table)
  data_table <- geo_data@dataTable@table
  
  # Create a data frame from the platform data
  database <- data.frame(ID = geo_data@dataTable@table$ID, GENE = geo_data@dataTable@table$GENE, GENE_SYMBOL = geo_data@dataTable@table$GENE_SYMBOL)
  
  #----------------------------------------------------
  
  # Normalize Data
  Train_normalized <- normalize.quantiles(Train)
  Test_normalized <- normalize.quantiles(Test)
  
  # Set the row and column names for the normalized data
  rownames(Train_normalized) <- rownames(Train)
  colnames(Train_normalized) <- colnames(Train)
  
  rownames(Test_normalized) <- rownames(Test)
  colnames(Test_normalized) <- colnames(Test)
  
  # Log2 transform the normalized data
  Train_normalized_log2 <- log2(Train_normalized)
  Test_normalized_log2 <- log2(Test_normalized)
  
  # Add disease status as col names for the Train_normalized_log2 dataset
  colnames(Train_normalized_log2) <- paste(colnames(Train_normalized_log2), statustrain, sep = "")
  
  # Add disease status as col names for the Test_normalized_log2 dataset
  colnames(Test_normalized_log2) <- paste(colnames(Test_normalized_log2), statustest, sep = "")
  
  # num of columns 
  num_samples <- ncol(Train_normalized_log2)
  nSamplesTest <- ncol(Test_normalized_log2)
  
  # creating batches 
  batch_indicator_train <- factor(rep(c("Batch1", "Batch2"), each = num_samples / 2))
  batch_indicator_test <- factor(rep(c("Batch1Test", "Batch2test"), each = nSamplesTest / 2))
  
  set.seed(123)  # For reproducibility
  random_order <- sample(1:num_samples)
  randomTest <- sample(1:nSamplesTest)
  batch_indicator_train <- batch_indicator_train[random_order]
  batch_indicator_test <- batch_indicator_test[randomTest] 
  
  # randomizing the test labels with the dataset (future use)
  status_Testshuffled = statustest[randomTest]
  
  # randomizing the train labels with the dataset (future use)
  status_Trainshuffled = statustrain[random_order]
  
  # Perform ComBat on the train and test datasets separately
  corrected_train <- ComBat(dat = Train_normalized_log2[, random_order], batch = batch_indicator_train)
  corrected_test <- ComBat(dat = Test_normalized_log2[, randomTest], batch = batch_indicator_test)
  
  # Mapping the rownames of corrected_train to gene symbols
  corrected_train_gene_symbols <- corrected_train
  rownames(corrected_train_gene_symbols) <- database$GENE_SYMBOL[match(rownames(corrected_train), database$ID)]
  
  # Mapping the rownames of corrected_test to gene symbols
  corrected_test_gene_symbols <- corrected_test
  rownames(corrected_test_gene_symbols) <- database$GENE_SYMBOL[match(rownames(corrected_test), database$ID)]
  
  # Calculate variance for each gene
  gene_variance <- apply(corrected_train_gene_symbols, 1, var)
  
  # For each duplicated gene symbol, retain the one with the highest variance
  corrected_train_gene_symbols <- corrected_train_gene_symbols[order(-gene_variance),]
  corrected_train_gene_symbols <- corrected_train_gene_symbols[!duplicated(rownames(corrected_train_gene_symbols)),]
  
  # Repeat for test data
  gene_variance_test <- apply(corrected_test_gene_symbols, 1, var)
  corrected_test_gene_symbols <- corrected_test_gene_symbols[order(-gene_variance_test),]
  corrected_test_gene_symbols <- corrected_test_gene_symbols[!duplicated(rownames(corrected_test_gene_symbols)),]
  
  #----------------------------------------------------
  
  # Create a design matrix for the linear model
  design_matrix <- model.matrix(~ status_Trainshuffled)
  
  # Perform differential expression analysis using limma
  fit <- lmFit(corrected_train_gene_symbols, design_matrix)
  fit <- eBayes(fit)
  
  # Define the criteria for DEGs
  # Were looking for genes with FDR < 0.05 and absolute log2 fold change > 2
  
  p_value_filter <- fit$p.value[, "status_Trainshuffled"] < 0.01
  coefficient_filter <- abs(fit$coefficients[, "status_Trainshuffled"]) > 1
  
  de_genes_indices <- which(p_value_filter & coefficient_filter)
  
  # Extract the differentially expressed genes from the corrected_train dataset
  de_genes_expression <- corrected_train_gene_symbols[de_genes_indices, ]
  
  #cat("Number of DEGs:", length(de_genes_indices), "\n")
  
  #----------------------------------------------------
  
  DEG_coefficient <- 3
  
  # Create a data frame for the volcano plot
  volcano_data <- data.frame(
    Gene = rownames(de_genes_expression),
    Coefficient = fit$coefficients[rownames(de_genes_expression), "status_Trainshuffled"],
    PValue = fit$p.value[rownames(de_genes_expression), "status_Trainshuffled"]
  )
  
  # Define the color scheme
  volcano_data$Color <- ifelse(volcano_data$Coefficient > DEG_coefficient, "green",
                               ifelse(volcano_data$Coefficient < -DEG_coefficient, "red", "grey"))
  
  # Create the volcano plot
  volcano_plot <- ggplot(volcano_data, aes(x = Coefficient, y = -log10(PValue), color = Color)) +
    geom_point() +
    scale_color_identity() +
    labs(x = "Log2 Fold Change", y = "-log10(P-value)",
         title = "Volcano Plot of Differentially Expressed Genes") +
    theme_minimal()
  
  #Labeling Genes
  significant_genes <- subset(volcano_data, PValue < 0.01 & (Coefficient > DEG_coefficient | Coefficient < -DEG_coefficient))
  print(volcano_plot + geom_text(data = significant_genes, aes(label = Gene), vjust = 1, hjust = 1, size = 3))
  
  # Filter genes that are up-regulated or down-regulated
  regulated_genes <- de_genes_expression[abs(volcano_data$Coefficient) > DEG_coefficient, ]
  
  # Select the expression values of the regulated genes
  regulated_gene_expression <- regulated_genes[, 1:ncol(regulated_genes)]
  
  # Create a heatmap
  pheatmap(regulated_gene_expression,
           color = colorRampPalette(c("blue", "white", "red"))(100),
           clustering_distance_cols = "euclidean",
           clustering_distance_rows = "euclidean",
           main = "Heatmap of Up/Down-Regulated Genes",
           fontsize = 8)
  
  #----------------------------------------------------
  
  
  # Extracting DEGs based on coefficient value
  DEGs <- volcano_data$Gene[abs(volcano_data$Coefficient) > DEG_coefficient]
  
  # Transpose the datasets
  transposed_train_data <- t(corrected_train_gene_symbols)
  transposed_test_data <- t(corrected_test_gene_symbols)
  
  # Filter to DEGs
  train_subset <- transposed_train_data[, DEGs]
  test_subset <- transposed_test_data[, DEGs]
  
  # Check for any missing genes in the test dataset and fill them with median values from the training data
  missing_genes <- setdiff(DEGs, colnames(test_subset))
  for (gene in missing_genes) {
    test_subset[, gene] <- median(train_subset[, gene], na.rm = TRUE)
  }
  
  # Train the Random Forest model
  set.seed(123)  # For reproducibility
  rf_model <- randomForest(x = train_subset, y = as.factor(status_Trainshuffled), importance=TRUE)
  
  # Extracting important genes from the model
  important_genes <- importance(rf_model)
  
  # Rank genes by their importance
  important_genes_sorted <- important_genes[order(-important_genes[, "MeanDecreaseGini"]),]
  
  # Add the importance to the data frame for the current iteration
  col_name <- paste("Iteration", i)
  gene_importance_df[col_name] <- important_genes_sorted[match(gene_importance_df$Gene, rownames(important_genes_sorted)), "MeanDecreaseGini"]

  # Predict using the Random Forest model on the test data
  rf_predictions <- predict(rf_model, test_subset)
  
  # Calculate the accuracy for the current iteration
  accuracy <- sum(rf_predictions == as.factor(status_Testshuffled)) / length(status_Testshuffled)
  
  # Print the predictions vs. actual
  cat("Iteration", i, "\n")
  cat("Predictions: ", rf_predictions, "\n")
  cat("Actual: ", as.factor(status_Testshuffled), "\n\n")
  
  # Store the accuracy
  accuracy_list[i] <- accuracy
  
  }

# After the loop is done, We sort the genes by their accumulated importance
gene_importance_df$TotalImportance <- rowSums(gene_importance_df[, -1], na.rm=TRUE)
gene_importance_df$AverageImportance <- gene_importance_df$TotalImportance / num_iterations
sorted_genes_by_average_importance <- gene_importance_df[order(-gene_importance_df$AverageImportance),]

# Sort genes by their total importance
sorted_genes_by_importance <- gene_importance_df[order(-gene_importance_df$TotalImportance),]
rownames(sorted_genes_by_average_importance) <- sorted_genes_by_average_importance$Gene
sorted_genes_by_average_importance$Gene <- NULL
sorted_genes_by_average_importance$TotalImportance <- NULL

# Display top genes based on their importance (change n to the # of)
head(sorted_genes_by_average_importance, n=20)

# Average all the accuracy rates to find a combined accuracy
combined_accuracy <- mean(accuracy_list)
cat("Average Accuracy of Entire Model: ", combined_accuracy, "\n")

#----------------------------------------------------

# Extracting top 20 important genes for visualization
top_20_genes <- head(sorted_genes_by_average_importance, 20)
barplot(top_20_genes$AverageImportance, 
        las = 2, 
        names.arg = rownames(top_20_genes), 
        main = "Top 20 Important Genes", 
        col = "steelblue", 
        cex.names = 0.7, 
        horiz = TRUE)


