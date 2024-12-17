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
with_pheromones <- read.csv("ant_with_pheromones.csv")
no_pheromones <- read.csv("ant_no_pheromones.csv")

# Assign group labels
with_pheromones$group <- "with_pheromones"
no_pheromones$group <- "no_pheromones"

# Combine the datasets into one data frame
df <- rbind(with_pheromones, no_pheromones)

# View the first few rows to verify
head(df)
str(df)

# -----------------------------
# 3. Verify column names
# -----------------------------

print(colnames(df))
# Expected Output: "run", "total_steps", "group"

# -----------------------------
# 4. Convert 'group' to factor
# -----------------------------

df$group <- factor(df$group, levels = c("with_pheromones", "no_pheromones"))
head(df)
str(df)

# -----------------------------
# 5. Summary statistics 
# -----------------------------

summary_stats <- df %>%
  group_by(group) %>%
  summarise(
    count = n(),
    mean_total_steps = mean(total_steps, na.rm = TRUE),
    sd_total_steps = sd(total_steps, na.rm = TRUE),
    median_total_steps = median(total_steps, na.rm = TRUE),
    IQR_total_steps = IQR(total_steps, na.rm = TRUE),
    se_total_steps = sd(total_steps, na.rm = TRUE) / sqrt(count)  
  )

print(summary_stats)

boxplot(df$total_steps ~ df$group, outline = TRUE)

# -----------------------------
# 6. Visual inspection of raw data
# -----------------------------

# Histogram and Density Plot
ggplot(df, aes(x = total_steps, fill = group)) +
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
shapiro_with <- shapiro.test(with_pheromones$total_steps)
shapiro_without <- shapiro.test(no_pheromones$total_steps)

print(shapiro_with)       # p < 2.2e-16 --> significant deviation from normality
print(shapiro_without)    # p < 2.2e-16 --> significant deviation from normality

# Levene's Test for homogeneity of variances
levene_result <- leveneTest(total_steps ~ group, data = df)
print(levene_result)     # p < 2.2e-16 --> unequal variances

# Wilcoxon Rank-Sum Test (non-parametric)
wilcox_test_result <- wilcox.test(total_steps ~ group, data = df, exact = FALSE)
print(wilcox_test_result)  # p-value < 2.2e-16 --> significant difference between groups

# Extract Z statistic (Z is derived from the p-value and W statistic)
Z <- qnorm(wilcox_test_result$p.value / 2, lower.tail = FALSE)  # Two-tailed test
Z 

# Total sample size
n1 <- sum(df$group == "with_pheromones")
n2 <- sum(df$group == "no_pheromones")
N <- n1 + n2

# Compute r
r <- Z / sqrt(N)
print(r)  # Standardized effect size: 0.7612419


# -----------------------------
# 8. Boxplot 
# -----------------------------

plot_summary <- summary_stats
print(plot_summary)

boxplot_final <- ggplot(df, aes(x = group, y = total_steps, fill = group)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +  # Boxplot without outliers
  geom_jitter(width = 0.2, alpha = 0.6, color = "black") +  # Jittered data points
  geom_errorbar(
    data = plot_summary, 
    aes(
      x = group, 
      ymin = mean_total_steps - se_total_steps, 
      ymax = mean_total_steps + se_total_steps
    ), 
    inherit.aes = FALSE,  # Ensure it doesn't inherit aesthetics from the global ggplot call
    width = 0.2, color = "black"
  ) +  # Error bars
  geom_point(
    data = plot_summary, 
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

print(boxplot_final)

# -----------------------------
# 9. Line plot
# -----------------------------

# Calculate common y-axis range
y_range <- range(df$total_steps, na.rm = TRUE)

# Create plot for "with_pheromones"
plot_with <- ggplot(subset(df, group == "with_pheromones"), aes(x = run, y = total_steps)) +
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
plot_without <- ggplot(subset(df, group == "no_pheromones"), aes(x = run, y = total_steps)) +
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
food_track_with <- read.csv('food_tracking_pheromones.csv')
food_track_no <- read.csv('food_tracking_no_pheromones.csv')

# Add a column to identify the condition (with or without pheromones)
food_track_with$condition <- "With pheromone trails"
food_track_no$condition <- "Without pheromone trails"

# Combine both datasets into one for faceting
food_track_all <- rbind(food_track_with, food_track_no)

# Determine the x-axis limits (based on the maximum number of steps in either dataset)
max_steps <- max(max(food_track_with$after_nr_steps), max(food_track_no$after_nr_steps))

# Create the plot
plot <- ggplot(food_track_all, aes(x = after_nr_steps, y = remaining_food, group = run)) +
  geom_line(aes(color = condition)) +  
  labs(
    x = "Number of steps taken by 10 ants",
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


