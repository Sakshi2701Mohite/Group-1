# Reading the given tsv file and set the gene identifiers as row names
data <- read.delim("https://raw.githubusercontent.com/ghazkha/Assessment4/refs/heads/main/gene_expression.tsv", row.names = 1)

# Displaying table of values for the first six genes.
head(data)

# New column which is the mean of the other columns
data$mean <- rowMeans(data)

# Displaying a table of mean values for the first six genes.
head(data)

# Listing the 10 genes with the highest mean expression
top_10_genes <- head(data[order(-data$mean),],10)

# Output
top_10_genes

# Determining the number of genes with a mean <10
geneBelow10 <- sum(data$mean < 10)

# Output
geneBelow10

# Creating the mean expression column
data$mean_expression <- rowMeans(data[, -1], na.rm = TRUE)

# Creating a histogram plot of the mean values
hist(data$mean_expression, 
     main = "Histogram of Mean Gene Expression Values", 
     xlab = "Mean Expression", 
     col = "blue", 
     border = "black", 
     breaks = 20)

# Imported the csv file
growth_data <- read.csv("https://raw.githubusercontent.com/ghazkha/Assessment4/refs/heads/main/growth_data.csv")

# Checking the column names
colnames(growth_data)

# Checking the structure of the data
str(growth_data)

# Calculating Mean and SD for control site at the start
mean_control_start <- mean(growth_data$Circumf_2005_cm, na.rm = TRUE)
sd_control_start <- sd(growth_data$Circumf_2005_cm, na.rm = TRUE)

# Displaying the output
cat("Control Site - Start: Mean =", mean_control_start, "SD =", sd_control_start, "\n")

# Calculating Mean and SD for control site at the end
mean_control_end <- mean(growth_data$Circumf_2010_cm, na.rm = TRUE)
sd_control_end <- sd(growth_data$Circumf_2010_cm, na.rm = TRUE)

# Displaying the output
cat("Control Site - End: Mean =", mean_control_end, "SD =", sd_control_end, "\n")

# Mean and SD for treatment site at the start
mean_treatment_start <- mean(growth_data$Circumf_2015_cm, na.rm = TRUE)
sd_treatment_start <- sd(growth_data$Circumf_2015_cm, na.rm = TRUE)

# Print the results
cat("Treatment Site - Start: Mean =", mean_treatment_start, "SD =", sd_treatment_start, "\n")

mean_treatment_end <- mean(growth_data$Circumf_2020_cm, na.rm = TRUE)
sd_treatment_end <- sd(growth_data$Circumf_2020_cm, na.rm = TRUE)

# Print the results
cat("Treatment Site - End: Mean =", mean_treatment_end, "SD =", sd_treatment_end, "\n")
