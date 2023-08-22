library(GEOquery)
library(preprocessCore)
library(sva)
library(randomForest) 
library(ggplot2)
library(limma)
library(AnnotationDbi)
library(pheatmap)
library(caret)
library(knitr)
library(glmnet)
library(pROC)

# Get DataSet from GEO Database
gse <- getGEO("GSE66494", GSEMatrix = TRUE)

# Get the Labels (possible future use)
gseData <- pData(gse[[1]])

# creating labels for the train dataset
statustrain <- c(
  rep(1, 47),
  rep(0, 5)
)

# creating labels for the test dataset
statustest <- c(
  rep(1, 5),
  rep(0, 3)
)

# Import Dataset and Save
Expr_GSM1623299 <- exprs(gse[[1]])

# Set column names for Train, Control, and Test datasets
sample_names <- colnames(Expr_GSM1623299)
Train_cols <- sample_names[2:53]
Test_cols <- sample_names[54:ncol(Expr_GSM1623299)]

# Splice Data
Train <- Expr_GSM1623299[, Train_cols, drop = FALSE]
Test <- Expr_GSM1623299[, Test_cols, drop = FALSE]

# Fetch the GEO dataset by accession number
geo_data <- getGEO("GPL6480")

# Access the expression data (data table)
data_table <- geo_data@dataTable@table

# Create a data frame from the platform data
database <- data.frame(ID = geo_data@dataTable@table$ID, GENE = geo_data@dataTable@table$GENE, GENE_SYMBOL = geo_data@dataTable@table$GENE_SYMBOL)

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

set.seed(123) # For reproducibility
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

View(corrected_train)
View(corrected_test)

# Create a design matrix for the linear model
design_matrix <- model.matrix(~ status_Trainshuffled)

# Perform differential expression analysis using limma
fit <- lmFit(corrected_train, design_matrix)
fit <- eBayes(fit)

# Define the criteria for DEGs
# We're looking for genes with FDR < 0.05 and absolute log2 fold change > 1
p_value_filter <- fit$p.value[, "status_Trainshuffled"] < 0.01
coefficient_filter <- abs(fit$coefficients[, "status_Trainshuffled"]) > 1

de_genes_indices <- which(p_value_filter & coefficient_filter)

# Extract the differentially expressed genes from the corrected_train dataset
de_genes_expression <- corrected_train[de_genes_indices, ]

cat("Number of DEGs:", length(de_genes_indices), "\n")

# Print the gene names and coefficients of DEGs
for (i in 1:nrow(de_genes_expression)) {
  gene_name <- rownames(de_genes_expression)[i]
  coefficient <- fit$coefficients[gene_name, "status_Trainshuffled"]
  cat("Gene:", gene_name, " | Coefficient:", coefficient, "\n")
}

# Create a data frame for the volcano plot
volcano_data <- data.frame(
  Gene = rownames(de_genes_expression),
  Coefficient = fit$coefficients[rownames(de_genes_expression), "status_Trainshuffled"],
  PValue = fit$p.value[rownames(de_genes_expression), "status_Trainshuffled"]
)

# Define the color scheme
volcano_data$Color <- ifelse(volcano_data$Coefficient > 4, "green",
                             ifelse(volcano_data$Coefficient < -4, "red", "grey"))

# Create the volcano plot with labeled genes
volcano_plot <- ggplot(volcano_data, aes(x = Coefficient, y = -log10(PValue), color = Color)) +
  geom_point() +
  geom_text(data = volcano_data[abs(volcano_data$Coefficient) > 4, ], 
            aes(label = Gene), hjust = 0.5, vjust = -1, size = 2.5, nudge_y = 0.2) +  # Adding gene labels
  scale_color_identity() +
  labs(x = "Log2 Fold Change", y = "-log10(P-value)",
       title = "Volcano Plot of Differentially Expressed Genes") +
  theme_minimal()

# Add labels to the green and red dots
volcano_plot_with_labels <- volcano_plot +
  geom_text(data = subset(volcano_data, Color %in% c("green", "red")),
            aes(label = Gene), hjust = 0, vjust = 0)

# Print the volcano plot with labeled genes
print(volcano_plot)

# Filter genes that are up-regulated or down-regulated
regulated_genes <- de_genes_expression[abs(volcano_data$Coefficient) > 4, ]

# Select the expression values of the regulated genes
regulated_gene_expression <- regulated_genes[, 1:ncol(regulated_genes)]

# Create a heatmap
pheatmap(regulated_gene_expression,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         clustering_distance_cols = "euclidean",
         clustering_distance_rows = "euclidean",
         main = "Heatmap of Up/Down-Regulated Genes",
         fontsize = 8)


# Sort DEGs by p-value and fold change
sorted_de_genes <- de_genes_expression[order(fit$p.value[de_genes_indices, "status_Trainshuffled"], decreasing = FALSE), ]

# Select the top 30 genes with the smallest p-values
top_de_genes <- sorted_de_genes[1:30, ]

# Print the number of selected top DEGs
cat("Number of Top DEGs:", nrow(top_de_genes), "\n")

# Print the gene symbols of the top DEGs
gene_symbols <- rownames(top_de_genes)
cat("Top Differentially Expressed Genes:\n")
cat(paste(gene_symbols, collapse = ", "), "\n")

# This column should contain gene symbols or IDs
gene_symbols <- rownames(de_genes_expression)

# If the gene symbols are the row names, you can directly use them
gene_names <- gene_symbols

# Feature Selection using Lasso Regression (glmnet package)
lasso_model <- cv.glmnet(as.matrix(t(corrected_train)), status_Trainshuffled, alpha = 1)
lasso_selected_genes <- coef(lasso_model, s = "lambda.min")

# Convert coefficient matrix to a vector of selected gene names
selected_gene_indices <- which(lasso_selected_genes != 0)
corrected_train_subset <- corrected_train[selected_gene_indices, ]

corrected_train_subset[is.na(corrected_train_subset)] <- 0

# Feature Selection using Recursive Feature Elimination (RFE) and Random Forest
# Note: This section requires the training dataset to be loaded as "corrected_train"
set.seed(123)  # For reproducibility
rfe_control <- rfeControl(functions = rfFuncs, method = "cv", number = 10)
num_features <- ncol(corrected_train_subset)
rfe_result <- rfe(t(corrected_train_subset), status_Trainshuffled, sizes = 1:num_features, rfeControl = rfe_control)

# Get the optimal features selected by RFE
selected_rfe_genes <- selected_rfe_genes[rfe_result$optVariables]

# Random Forest for Disease Risk Prediction Model
# Grid Search for the best mtry parameter
set.seed(123)  # For reproducibility
mtry_grid <- expand.grid(mtry = c(2, 5, 10, 20))
rf_model <- train(corrected_train[, selected_rfe_genes], status_Trainshuffled, method = "rf", trControl = trainControl(method = "cv", number = 10), tuneGrid = mtry_grid)

# Optimal mtry value
optimal_mtry <- rf_model$bestTune$mtry

# Building the Random Forest Model
set.seed(123)  # For reproducibility
final_rf_model <- randomForest(corrected_train[, selected_rfe_genes], status_Trainshuffled, mtry = optimal_mtry)

# Model Evaluation on Training Data
predicted_status <- predict(final_rf_model, corrected_train[, selected_rfe_genes])
confusion_matrix <- confusionMatrix(predicted_status, status_Trainshuffled)
accuracy <- confusion_matrix$overall["Accuracy"]
auc <- roc(status_Trainshuffled, final_rf_model$votes[, "1"])$auc

# Model Evaluation on Validation Data
# Assuming you have loaded and preprocessed validation datasets
predicted_status_val <- predict(final_rf_model, corrected_validation[, selected_rfe_genes])
confusion_matrix_val <- confusionMatrix(predicted_status_val, status_validation)
accuracy_val <- confusion_matrix_val$overall["Accuracy"]
auc_val <- roc(status_validation, final_rf_model$votes_val[, "1"])$auc

# Print evaluation results
cat("Training Data:\n")
cat("Accuracy:", accuracy, "\n")
cat("AUC:", auc, "\n\n")

cat("Validation Data:\n")
cat("Accuracy:", accuracy_val, "\n")
cat("AUC:", auc_val, "\n")
