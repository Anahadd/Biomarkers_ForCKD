options(scipen = 999)

library(GEOquery)
library(preprocessCore)
library(sva)
library(randomForest)  # Add the randomForest package
library(pheatmap)
library(ggplot2)
library(limma)
library(glmnet)
library(caret)
library(dplyr)
library(caret)

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

#----------------------

# Fetch the GEO dataset by accession number
geo_data <- getGEO("GPL6480")

# Access the expression data (data table)
data_table <- geo_data@dataTable@table

# Create a data frame from the platform data
database <- data.frame(ID = geo_data@dataTable@table$ID, GENE = geo_data@dataTable@table$GENE, GENE_SYMBOL = geo_data@dataTable@table$GENE_SYMBOL)

#------------------------

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

#---------------------

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

#--------------------

# Create a data frame for the volcano plot
volcano_data <- data.frame(
  Gene = rownames(de_genes_expression),
  Coefficient = fit$coefficients[rownames(de_genes_expression), "status_Trainshuffled"],
  PValue = fit$p.value[rownames(de_genes_expression), "status_Trainshuffled"]
)

# Define the color scheme
volcano_data$Color <- ifelse(volcano_data$Coefficient > 3.5, "green",
                             ifelse(volcano_data$Coefficient < -3.5, "red", "grey"))

# Create the volcano plot
volcano_plot <- ggplot(volcano_data, aes(x = Coefficient, y = -log10(PValue), color = Color)) +
  geom_point() +
  scale_color_identity() +
  labs(x = "Log2 Fold Change", y = "-log10(P-value)",
       title = "Volcano Plot of Differentially Expressed Genes") +
  theme_minimal()

# Add labels to the green and red dots
#volcano_plot_with_labels <- volcano_plot +
#geom_text(data = subset(volcano_data, Color %in% c("green", "red")),
#aes(label = Gene), hjust = 0, vjust = 0)

# Print the volcano plot
print(volcano_plot)

# Filter genes that are up-regulated or down-regulated
regulated_genes <- de_genes_expression[abs(volcano_data$Coefficient) > 3.5, ]

# Select the expression values of the regulated genes
regulated_gene_expression <- regulated_genes[, 1:ncol(regulated_genes)]

# Create a heatmap
pheatmap(regulated_gene_expression,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         clustering_distance_cols = "euclidean",
         clustering_distance_rows = "euclidean",
         main = "Heatmap of Up/Down-Regulated Genes",
         fontsize = 8)

#---------------------------------

# Extracting DEGs based on coefficient value
DEGs <- volcano_data$Gene[abs(volcano_data$Coefficient) > 3.5]

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

# Test the Random Forest model
rf_predictions <- predict(rf_model, test_subset)

# Print model predictions and actual status side by side for validation set
predictions_df <- data.frame(Predicted = rf_predictions, Actual = status_Testshuffled)
print(predictions_df)

# Calculate Accuracy
accuracy <- sum(rf_predictions == as.factor(status_Testshuffled)) / length(status_Testshuffled)
print(paste("Accuracy on test data:", round(accuracy * 100, 2), "%"))

#----------------------------------------------------

# If you want to see the importance of each gene in the model:
important_genes <- importance(rf_model)
important_genes_sorted <- important_genes[order(-important_genes[, "MeanDecreaseGini"]),]

# Extracting top 20 important genes for visualization
top_20_genes <- head(important_genes_sorted, 20)
barplot(top_20_genes[,"MeanDecreaseGini"], 
        las = 2, 
        names.arg = rownames(top_20_genes), 
        main = "Top 20 Important Genes", 
        col = "steelblue", 
        cex.names = 0.7, 
        horiz = TRUE)

#----------------------------------------------------

# Set up k-fold cross-validation
control <- trainControl(method="cv", 
                        number=10,  # 10-fold CV
                        search="grid",
                        classProbs = TRUE,
                        summaryFunction=twoClassSummary)

# Ensure response variable is a factor
status_Trainshuffled <- as.factor(status_Trainshuffled)
levels(status_Trainshuffled) <- c("Class0", "Class1")

# Train the model with cross-validation
set.seed(123)
model_cv <- train(train_subset, status_Trainshuffled, 
                  method="rf", 
                  trControl=control, 
                  metric="Accuracy",
                  importance=TRUE)

# Predict using the cross-validated model
rf_predictions_CV <- predict(model_cv, test_subset)

# Convert to factor and ensure consistent levels
rf_predictions_CV <- factor(rf_predictions_CV, levels=c("Class0", "Class1"))

# Convert 0 to Class0 and 1 to Class1 in status_Testshuffled
status_Testshuffled <- ifelse(status_Testshuffled == 0, "Class0", "Class1")

# Print model predictions and actual status side by side for validation set
predictions_df_CV <- data.frame(Predicted = rf_predictions_CV, Actual = status_Testshuffled)
print(predictions_df_CV)

# Now calculate the accuracy
accuracy_cv <- sum(rf_predictions_CV == status_Testshuffled) / length(status_Testshuffled)
print(paste("Accuracy on test data with CV:", round(accuracy_cv * 100, 2), "%"))

#----------------------------------------------------

# Extracting important genes for the cross-validated model
important_genes_cv <- model_cv$finalModel$importance
important_genes_sorted_cv <- important_genes_cv[order(-important_genes_cv[, "MeanDecreaseGini"]),]

# Extracting top 20 important genes for visualization
top_20_genes_cv <- head(important_genes_sorted_cv, 20)
barplot(top_20_genes_cv[,"MeanDecreaseGini"], 
        las = 2, 
        names.arg = rownames(top_20_genes_cv), 
        main = "Top 20 Important Genes (Cross-validated Model)", 
        col = "slateblue3", 
        cex.names = 0.7, 
        horiz = TRUE)

