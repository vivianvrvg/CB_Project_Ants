#########################################################################
## Title: ANT FORAGING BEHAVIOUR - Statistics                          ##
## Course: Project Computational Biology                               ##
## Authors: Camille Timmermans & Vivian Van Reybrouck Van Gelder       ##
#########################################################################

# -----------------------------
# 1. Load libraries
# -----------------------------

library(ggplot2)
library(dplyr)
library(car)       # For Levene's Test
library(ggsignif)  # For significance annotations
library(gridExtra) # For arranging multiple plots

# -----------------------------
# 2. Load and prepare data
# -----------------------------

# Load data from CSV files
with_pheromones <- read.csv("ant_with_pheromones.csv")
no_pheromones <- read.csv("ant_without_pheromones.csv")

# Assign group labels
with_pheromones$group <- "with_pheromones"
no_pheromones$group <- "no_pheromones"

# Combine the datasets into one data frame
df <- rbind(with_pheromones, no_pheromones)

# View the first few rows to verify
head(df)

# -----------------------------
# 3. Verify column names
# -----------------------------

# Check column names to ensure 'remaining_food' exists
print(colnames(df))
# Expected Output: "run", "remaining_food", "group"

# -----------------------------
# 4. Convert 'group' to factor
# -----------------------------

df$group <- factor(df$group, levels = c("with_pheromones", "no_pheromones"))

# -----------------------------
# 5. Summary statistics 
# -----------------------------

summary_stats <- df %>%
  group_by(group) %>%
  summarise(
    count = n(),
    mean_remaining_food = mean(remaining_food, na.rm = TRUE),
    sd_remaining_food = sd(remaining_food, na.rm = TRUE),
    median_remaining_food = median(remaining_food, na.rm = TRUE),
    IQR_remaining_food = IQR(remaining_food, na.rm = TRUE),
    se_remaining_food = sd_remaining_food / sqrt(count)
  )

print(summary_stats)

# -----------------------------
# 6. Visual inspection of raw data
# -----------------------------

# Histogram and Density Plot
ggplot(df, aes(x = remaining_food, fill = group)) +
  geom_histogram(alpha = 0.6, position = "identity", bins = 15) +
  geom_density(alpha = 0.2) +
  theme_minimal() +
  labs(title = "Distribution of Food Left by Group",
       x = "Food Left",
       y = "Frequency") +
  scale_fill_manual(values = c("with_pheromones" = "#1f77b4", "no_pheromones" = "#ff7f0e")) +
  theme(legend.title = element_blank())

# -----------------------------
# 7. Statistical testing
# -----------------------------

# Shapiro-Wilk Test for normality
shapiro_with <- shapiro.test(with_pheromones$remaining_food)
shapiro_without <- shapiro.test(no_pheromones$remaining_food)

print(shapiro_with)       # p = 0.007022 --> significant deviation from normality
print(shapiro_without)    # p = 0.5381 --> no significant deviation from normality

# Levene's Test for homogeneity of variances
levene_result <- leveneTest(remaining_food ~ group, data = df)
print(levene_result)     # p = 0.09946 --> equal variances

# Wilcoxon Rank-Sum Test (non-parametric)
wilcox_test_result <- wilcox.test(remaining_food ~ group, data = df)
print(wilcox_test_result)   # p-value < 2.2e-16 --> significant difference between groups

# -----------------------------
# 8. Boxplot 
# -----------------------------

plot_summary <- summary_stats
print(plot_summary)

boxplot_final <- ggplot(df, aes(x = group, y = remaining_food, fill = group)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +  # Boxplot without outliers
  geom_jitter(width = 0.2, alpha = 0.6, color = "black") +  # Jittered data points
  geom_errorbar(data = plot_summary, 
                aes(x = group, 
                    ymin = mean_remaining_food - se_remaining_food, 
                    ymax = mean_remaining_food + se_remaining_food), 
                inherit.aes = FALSE,  # Ensure it doesn't inherit aesthetics from the global ggplot call
                width = 0.2, color = "black") +  # Error bars
  geom_point(data = plot_summary, 
             aes(x = group, y = mean_remaining_food), 
             inherit.aes = FALSE,  # Prevent inherited aesthetics
             shape = 23, size = 3, fill = "white") +  # Mean points
  theme_minimal() +
  labs(title = "Comparison of Food Left by Ants With and Without Pheromones",
       x = "Group",
       y = "Food Left") +
  scale_fill_manual(values = c("with_pheromones" = "#1f77b4", "no_pheromones" = "#ff7f0e")) +  # Colors
  theme(
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 14)
  )


print(boxplot_final)

# -----------------------------
# 9. Line plot
# -----------------------------

# Create plot for "with_pheromones"
plot_with <- ggplot(subset(df, group == "with_pheromones"), aes(x = run, y = remaining_food)) +
  geom_line(color = "#1f77b4", size = 1) +  # Blue line
  theme_minimal() +
  labs(title = "With Pheromones",
       x = "Run Number",
       y = "Food Left") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 12)
  )

# Create plot for "no_pheromones"
plot_without <- ggplot(subset(df, group == "no_pheromones"), aes(x = run, y = remaining_food)) +
  geom_line(color = "#ff7f0e", size = 1) +  # Orange line
  theme_minimal() +
  labs(title = "Without Pheromones",
       x = "Run Number",
       y = "Food Left") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 12)
  )

# Arrange the two line plots side by side
grid.arrange(plot_with, plot_without, ncol = 2)
