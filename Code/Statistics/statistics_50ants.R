#########################################################################
## Title: ANT FORAGING BEHAVIOUR - Statistics baseline                 ##
## Course: Project Computational Biology                               ##
#########################################################################

# -----------------------------
# 1. Load libraries
# -----------------------------

library(ggplot2)
library(dplyr)
library(car)       # For Levene's Test
library(ggsignif)  # For significance annotations
library(gridExtra) # For arranging multiple plots
library(grid)

# -----------------------------
# 2. Load and prepare data
# -----------------------------

# Load data from CSV files
ant50_with_pheromones <- read.csv("50ants_with_pheromones.csv")
ant50_no_pheromones <- read.csv("50ants_no_pheromones.csv")

# Assign group labels
ant50_with_pheromones$group <- "with_pheromones"
ant50_no_pheromones$group <- "no_pheromones"

# Combine the datasets into one data frame
ant50 <- rbind(ant50_with_pheromones, ant50_no_pheromones)

# View the first few rows to verify
head(ant50)
str(ant50)

# -----------------------------
# 3. Verify column names
# -----------------------------

print(colnames(ant50))
# Expected Output: "run", "total_steps", "group"

# -----------------------------
# 4. Convert 'group' to factor
# -----------------------------

ant50$group <- factor(ant50$group, levels = c("with_pheromones", "no_pheromones"))
head(ant50)
str(ant50)

# -----------------------------
# 5. Summary statistics 
# -----------------------------

summary_stats_ant50 <- ant50 %>%
  group_by(group) %>%
  summarise(
    count = n(),
    mean_total_steps = mean(total_steps, na.rm = TRUE),
    sd_total_steps = sd(total_steps, na.rm = TRUE),
    median_total_steps = median(total_steps, na.rm = TRUE),
    IQR_total_steps = IQR(total_steps, na.rm = TRUE),
    se_total_steps = sd(total_steps, na.rm = TRUE) / sqrt(count)  
  )

print(summary_stats_ant50)

boxplot(ant50$total_steps ~ ant50$group, outline = TRUE)

# -----------------------------
# 6. Visual inspection of raw data
# -----------------------------

# Histogram and Density Plot
ggplot(ant50, aes(x = total_steps, fill = group)) +
  geom_histogram(alpha = 0.6, position = "identity", bins = 15) +
  geom_density(alpha = 0.2) +
  theme_minimal() +
  labs(title = "Total steps until food depletion",
       x = "Total steps",
       y = "Frequency") +
  scale_fill_manual(values = c("with_pheromones" = "#1f77b4", "no_pheromones" = "#ff7f0e")) +
  theme(legend.title = element_blank())

# -----------------------------
# 7. Statistical testing
# -----------------------------

# Shapiro-Wilk Test for normality
shapiro_with_ant50 <- shapiro.test(ant50_with_pheromones$total_steps)
shapiro_without_ant50 <- shapiro.test(ant50_no_pheromones$total_steps)

print(shapiro_with_ant50)       # p = 1.797e-11 --> significant deviation from normality
print(shapiro_without_ant50)    # p < 2.2e-16 --> significant deviation from normality

# Levene's Test for homogeneity of variances
levene_result_ant50 <- leveneTest(total_steps ~ group, data = ant50)
print(levene_result_ant50)     # p < 2.2e-16 --> unequal variances

# Wilcoxon Rank-Sum Test (non-parametric)
wilcox_test_result_ant50 <- wilcox.test(total_steps ~ group, data = ant50, exact = FALSE)
print(wilcox_test_result_ant50)  # p-value < 2.2e-16 --> significant difference between groups

# Group sizes
n1 <- sum(ant50$group == "with_pheromones")
n2 <- sum(ant50$group == "no_pheromones")
N <- n1 + n2

# Compute mean and standard deviation of W
mu_W <- (n1 * n2) / 2
sigma_W <- sqrt((n1 * n2 * (n1 + n2 + 1)) / 12)

# Compute Z
W <- wilcox_test_result_ant50$statistic
Z <- (W - mu_W) / sigma_W
Z # -37.94997

# Compute r
r <- Z / sqrt(N)
print(r) # -0.8485872 --> strong effect

# -----------------------------
# 8. Boxplot 
# -----------------------------

plot_summary_ant50 <- summary_stats_ant50
print(plot_summary_ant50)

boxplot_final_ant50 <- ggplot(ant50, aes(x = group, y = total_steps, fill = group)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +  # Boxplot without outliers
  geom_jitter(width = 0.2, alpha = 0.6, color = "black") +  # Jittered data points
  geom_errorbar(
    data = plot_summary_ant50, 
    aes(
      x = group, 
      ymin = mean_total_steps - se_total_steps, 
      ymax = mean_total_steps + se_total_steps
    ), 
    inherit.aes = FALSE,  # Ensure it doesn't inherit aesthetics from the global ggplot call
    width = 0.2, color = "black"
  ) +  # Error bars
  geom_point(
    data = plot_summary_ant50, 
    aes(x = group, y = mean_total_steps), 
    inherit.aes = FALSE,  # Prevent inherited aesthetics
    shape = 23, size = 3, fill = "white"
  ) +  # Mean points
  scale_x_discrete(labels = c(
    "with_pheromones" = "With pheromone trails",
    "no_pheromones" = "Without pheromone trails"
  )) +  # Rename x-axis categories
  scale_fill_manual(values = c(
    "with_pheromones" = "#1f77b4", 
    "no_pheromones" = "#ff7f0e"
  )) +  # Set custom colors
  theme_bw() +
  labs(y = "Total steps until food depletion") +
  theme(
    legend.position = "none",                      # Remove legend
    plot.title = element_blank(),                  # Remove title
    axis.title.x = element_blank(),                # Remove x-axis label
    axis.title.y = element_text(size = 14),        # Keep and enlarge y-axis label
    axis.text = element_text(size = 12),           # Enlarge axis tick values
    strip.text = element_text(size = 14, face = "bold")  # Enlarge facet condition labels
  )

print(boxplot_final_ant50)

# -----------------------------
# 9. Line plot
# -----------------------------

# Calculate common y-axis range
y_range <- range(ant50$total_steps, na.rm = TRUE)

# Create plot for "with_pheromones"
plot_with <- ggplot(subset(ant50, group == "with_pheromones"), aes(x = run, y = total_steps)) +
  geom_line(color = "#1f77b4", size = 1) +  # Blue line
  theme_minimal() +
  labs(title = "With Pheromones",
       x = "Run Number",
       y = "Total Steps") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 12)
  ) +
  ylim(y_range)

# Create plot for "no_pheromones"
plot_without <- ggplot(subset(ant50, group == "no_pheromones"), aes(x = run, y = total_steps)) +
  geom_line(color = "#ff7f0e", size = 1) +  # Orange line
  theme_minimal() +
  labs(title = "Without Pheromones",
       x = "Run Number",
       y = "Total Steps") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 12)
  ) +
  ylim(y_range)

# Arrange the two line plots side by side
grid.arrange(plot_with, plot_without, ncol = 2)

# ---------------------------------
# 9. Line plot for every 100 steps
# ---------------------------------

# Load the datasets
ant50_food_track_with <- read.csv('food_tracking_pheromones_50ants.csv')
ant50_food_track_no <- read.csv('food_tracking_no_pheromones_50ants.csv')

# Add a column to identify the condition (with or without pheromones)
ant50_food_track_with$condition <- "With pheromone trails"
ant50_food_track_no$condition <- "Without pheromone trails"

# Combine both datasets into one for faceting
ant50_food_track_all <- rbind(ant50_food_track_with, ant50_food_track_no)

# Determine the x-axis limits (based on the maximum number of steps in either dataset)
max_steps <- max(max(ant50_food_track_with$after_nr_steps), max(ant50_food_track_no$after_nr_steps))

# Create the plot
plot <- ggplot(ant50_food_track_all, aes(x = after_nr_steps, y = remaining_food, group = run)) +
  geom_line(aes(color = condition)) +  
  labs(
    x = "Number of steps taken by 50 ants",
    y = "Remaining food in arena (20x20)"
  ) +
  scale_x_continuous(limits = c(0, max_steps)) +  
  scale_color_manual(values = c(
    "With pheromone trails" = "#1f77b4",
    "Without pheromone trails" = "#ff7f0e"
  )) +  # Set the custom colors
  facet_wrap(~condition, scales = "free_y", ncol = 2, labeller = label_value) +
  theme_bw() +  # Use a basic theme
  theme(
    axis.title = element_text(size = 14),      # Larger, bold axis labels
    axis.text = element_text(size = 12),                     # Larger axis values
    strip.text = element_text(size = 14),      # Larger facet condition labels
    legend.position = "none"                                 # Remove legend
  )

print(plot)


