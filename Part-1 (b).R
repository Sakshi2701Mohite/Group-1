# Imported the csv file
growth_data <- read.csv("https://raw.githubusercontent.com/ghazkha/Assessment4/refs/heads/main/growth_data.csv")

# Checking the column names
colnames(growth_data)

# Checking the structure of the data
str(growth_data)

# Installing the ggplot2 library
install.packages("ggplot2")

# Loading the ggplot2 library
library(ggplot2)

# Reshape the data for plotting of control site
circumference_data <- data.frame(
  site = rep(growth_data$Site, 2),
  period = rep(c("start", "end"), each = nrow(growth_data)),
  circumference = c(growth_data$Circumf_2005_cm, growth_data$Circumf_2010_cm)
)

# Assuming 'Site' is the correct column for site data
# Reshape the data for plotting
circumference_data <- data.frame(
  site = rep(growth_data$Site, 2),
  period = rep(c("start", "end"), each = nrow(growth_data)),
  circumference = c(growth_data$Circumf_2005_cm, growth_data$Circumf_2010_cm)
)

# Verify the reshaped data
head(circumference_data)

# Create the boxplot
ggplot(circumference_data, aes(x = interaction(site, period), y = circumference, fill = site)) +
  geom_boxplot() +
  labs(title = "Tree Circumference at Start and End of Study",
       x = "Site and Period", y = "Circumference (cm)") +
  theme_minimal()


# Reshape the data for plotting of treatment site
circumference_data <- data.frame(
  site = rep(growth_data$Site, 2),
  period = rep(c("start", "end"), each = nrow(growth_data)),
  circumference = c(growth_data$Circumf_2015_cm, growth_data$Circumf_2020_cm)
)

# Assuming 'Site' is the correct column for site data
# Reshape the data for plotting
circumference_data <- data.frame(
  site = rep(growth_data$Site, 2),
  period = rep(c("start", "end"), each = nrow(growth_data)),
  circumference = c(growth_data$Circumf_2015_cm, growth_data$Circumf_2020_cm)
)

# Verify the reshaped data
head(circumference_data)

# Create the boxplot
ggplot(circumference_data, aes(x = interaction(site, period), y = circumference, fill = site)) +
  geom_boxplot() +
  labs(title = "Tree Circumference at Start and End of Study",
       x = "Site and Period", y = "Circumference (cm)") +
  theme_minimal()

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

# Calculate the growth over the last 10 years
growth_control <- mean_control_end - mean_control_start
growth_treatment <- mean_treatment_end - mean_treatment_start

# Mean growth at each site
mean_growth_control <- mean(growth_control, na.rm = TRUE)
mean_growth_treatment <- mean(growth_treatment, na.rm = TRUE)

# Calculating the 10-year control site growth
growth_data$Control_site <- growth_data$Circumf_2015_cm - growth_data$Circumf_2005_cm

# Calculating the 10-year treatment site growth
growth_data$treatment_site <- growth_data$Circumf_2020_cm - growth_data$Circumf_2010_cm

# Print the results
growth_data$Control_site
growth_data$treatment_site

# Perform t-test
t_test_result <- t.test(growth_data$Control_site, growth_data$treatment_site)

# Print the results
print(t_test_result)
