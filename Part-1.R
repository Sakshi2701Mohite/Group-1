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
