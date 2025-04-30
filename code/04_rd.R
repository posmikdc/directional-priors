# ==============================================================================
# Title: Manifold RD Estimator Simulations
# Original Author: Daniel Posmik
# E-Mail: daniel_posmik@brown.edu
# ==============================================================================

# Setup ========================================================================
library(dplyr)
library(ggplot2)
library(circular)
library(latex2exp)

# Set working directory
setwd("~/Documents/brown/classes/spring25/php2530-bayes/directional-priors")

# Set Seed
set.seed(1)

# Simulation Parameters ========================================================
n_sims <- 500       # Number of simulation runs
n_pre <- 500        # Number of pre-treatment observations
n_post <- 500       # Number of post-treatment observations
t_star <- n_pre     # Cutoff point (time of treatment)

# Pre-treatment parameters (unit circle)
mean_pre <- pi/2
kappa_pre <- 5      # Concentration parameter (inverse of variance)
radius_pre <- 1     # Radius of pre-treatment manifold (unit circle)

# Post-treatment parameters (larger circle)
mean_post <- pi
kappa_post <- 5     # Same concentration for fair comparison
radius_post <- 2    # Radius of post-treatment manifold

# True treatment effect
true_effect <- mean_post - mean_pre

# Derived parameters
curvature_pre <- 1/radius_pre^2   # Gaussian curvature for pre-treatment sphere
curvature_post <- 1/radius_post^2 # Gaussian curvature for post-treatment sphere

# Function to generate data from von Mises distribution and add manifold properties
generate_data <- function(n, mean_angle, kappa, radius, time, is_post_treatment) {
  # Generate angles from von Mises distribution
  angles <- as.numeric(rvonmises(n, mu = circular(mean_angle), kappa = kappa))
  
  # Calculate Gaussian curvature
  curvature <- 1/radius^2
  
  # Create data frame
  data <- data.frame(
    time = time,
    angle = angles,
    radius = radius,
    curvature = curvature,
    post_treatment = as.integer(is_post_treatment),
    x = radius * cos(angles),  # Cartesian coordinates for visualization
    y = radius * sin(angles)
  )
  
  return(data)
}

# Estimator Function ==========================================================
# Standard RD estimator
rd_standard <- function(data, bandwidth = NULL) {
  if (!is.null(bandwidth)) {
    # Use bandwidth for local linear regression
    pre_data <- data %>% 
      filter(post_treatment == 0) %>%
      filter(time >= t_star - bandwidth)
    
    post_data <- data %>% 
      filter(post_treatment == 1) %>%
      filter(time <= t_star + bandwidth)
  } else {
    # Use all data
    pre_data <- data %>% filter(post_treatment == 0)
    post_data <- data %>% filter(post_treatment == 1)
  }
  
  # Calculate means
  mean_pre <- mean(pre_data$angle)
  mean_post <- mean(post_data$angle)
  
  # Treatment effect is difference in means
  effect <- mean_post - mean_pre
  
  return(effect)
}

# Run Simulation ===============================================================
run_simulation <- function() {
  results <- data.frame(
    sim = 1:n_sims,
    rd_estimate = numeric(n_sims)
  )
  
  for (i in 1:n_sims) {
    # Generate data
    pre_treatment <- generate_data(
      n = n_pre, 
      mean_angle = mean_pre, 
      kappa = kappa_pre, 
      radius = radius_pre, 
      time = 1:n_pre, 
      is_post_treatment = FALSE
    )
    
    post_treatment <- generate_data(
      n = n_post, 
      mean_angle = mean_post, 
      kappa = kappa_post, 
      radius = radius_post, 
      time = (t_star + 1):(t_star + n_post), 
      is_post_treatment = TRUE
    )
    
    # Combine datasets
    full_data <- rbind(pre_treatment, post_treatment)
    
    # Calculate treatment effect
    results$rd_estimate[i] <- rd_standard(full_data)
  }
  
  return(results)
}

# Additional simulation with varying kappa values
run_simulation_varying_kappa <- function() {
  kappa_values <- c(0.5, 1, 2, 5, 10, 20)
  results <- data.frame(
    kappa = rep(kappa_values, each = n_sims),
    sim = rep(1:n_sims, times = length(kappa_values)),
    rd_estimate = numeric(n_sims * length(kappa_values))
  )
  
  row_idx <- 1
  for (k in kappa_values) {
    for (i in 1:n_sims) {
      # Generate data with current kappa value
      pre_treatment <- generate_data(
        n = n_pre, 
        mean_angle = mean_pre, 
        kappa = k,  # Varying kappa
        radius = radius_pre, 
        time = 1:n_pre, 
        is_post_treatment = FALSE
      )
      
      post_treatment <- generate_data(
        n = n_post, 
        mean_angle = mean_post, 
        kappa = k,  # Varying kappa
        radius = radius_post, 
        time = (t_star + 1):(t_star + n_post), 
        is_post_treatment = TRUE
      )
      
      # Combine datasets
      full_data <- rbind(pre_treatment, post_treatment)
      
      # Calculate treatment effect
      results$rd_estimate[row_idx] <- rd_standard(full_data)
      
      row_idx <- row_idx + 1
    }
  }
  
  return(results)
}

# Execute simulations
results <- run_simulation()
results_varying_kappa <- run_simulation_varying_kappa()

# Visualizations ===============================================================

# 1. Visualize a single data sample
visualize_sample_data <- function() {
  # Generate one sample
  pre_treatment <- generate_data(
    n = n_pre, 
    mean_angle = mean_pre, 
    kappa = kappa_pre, 
    radius = radius_pre, 
    time = 1:n_pre, 
    is_post_treatment = FALSE
  )
  
  post_treatment <- generate_data(
    n = n_post, 
    mean_angle = mean_post, 
    kappa = kappa_post, 
    radius = radius_post, 
    time = (t_star + 1):(t_star + n_post), 
    is_post_treatment = TRUE
  )
  
  # Calculate mean angles for pre and post treatment
  mean_pre_x <- radius_pre * cos(mean_pre)
  mean_pre_y <- radius_pre * sin(mean_pre)
  mean_post_x <- radius_post * cos(mean_post)
  mean_post_y <- radius_post * sin(mean_post)
  
  # Calculate angle difference (treatment effect) in radians
  effect_radians <- mean_post - mean_pre
  
  # Plot
  ggplot() +
    # Pre-treatment data
    geom_point(data = pre_treatment, 
               aes(x = x, y = y), 
               color = "black", alpha = 0.4, size = 1) +
    # Post-treatment data
    geom_point(data = post_treatment, 
               aes(x = x, y = y), 
               color = "darkorange", alpha = 0.4, size = 1) +
    # Pre-treatment circle
    geom_path(data = data.frame(
      angle = seq(0, 2*pi, length.out = 100),
      x = radius_pre * cos(seq(0, 2*pi, length.out = 100)),
      y = radius_pre * sin(seq(0, 2*pi, length.out = 100))
    ), aes(x = x, y = y), color = "black", linetype = "dashed") +
    # Post-treatment circle
    geom_path(data = data.frame(
      angle = seq(0, 2*pi, length.out = 100),
      x = radius_post * cos(seq(0, 2*pi, length.out = 100)),
      y = radius_post * sin(seq(0, 2*pi, length.out = 100))
    ), aes(x = x, y = y), color = "red", linetype = "dashed") +
    # Mean points
    geom_point(aes(x = mean_pre_x, y = mean_pre_y), color = "darkgrey", size = 3) +
    geom_point(aes(x = mean_post_x, y = mean_post_y), color = "red", size = 3) +
    # Connect mean points to origin
    geom_segment(aes(x = 0, y = 0, xend = mean_pre_x, yend = mean_pre_y), 
                 arrow = arrow(length = unit(0.2, "cm")), color = "black") +
    geom_segment(aes(x = 0, y = 0, xend = mean_post_x, yend = mean_post_y), 
                 arrow = arrow(length = unit(0.2, "cm")), color = "red") +
    # Annotations
    annotate("text", x = mean_pre_x * 1.1, y = mean_pre_y * 1.2, 
             label = TeX("$\\bar{\\theta}^{~t<t^*}$"), color = "black", size = 4) +
    annotate("text", x = mean_post_x * 1.1, y = mean_post_y * 1.1, 
             label = TeX("$\\bar{\\theta}^{~t \\geq t^*}$"), color = "red", size = 4) +
    # Add annotation for effect size (π/2)
    annotate("text", x = 0, y = -0.4, 
             label = TeX(paste0("$\\tau = $", round(effect_radians, 2), " (deg: ", 
                                round(effect_radians * 180/pi, 2), "$^{\\circ}$)")), 
             size = 4) +
    # Theme and labels
    theme_minimal() +
    labs(
      title = "Canonical RD: Pre- and Post-Treatment Data",
      subtitle = "Black: Pre-treatment (r=1), Orange: Post-treatment (r=2)",
      x = "x",
      y = "y"
    ) +
    coord_fixed(ratio = 1) +  # Equal aspect ratio
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 12),
      axis.title = element_text(size = 12)
    )
}

# 2. Distribution of treatment effect estimates
plot_effect_distribution <- function(results) {
  # Calculate mean and standard deviation
  mean_estimate <- mean(results$rd_estimate)
  sd_estimate <- sd(results$rd_estimate)
  
  # Create density plot
  ggplot(results, aes(x = rd_estimate)) +
    geom_density(fill = "darkblue", alpha = 0.5) +
    geom_vline(xintercept = true_effect, linetype = "dashed", color = "black") +
    # Add mean line
    geom_vline(xintercept = mean_estimate, color = "blue", linetype = "dotted") +
    # Annotations
    annotate("text", x = true_effect, y = 5, 
             label = paste("True Effect =", round(true_effect, 2)), 
             angle = 90, vjust = -0.5) +
    # Formatting
    labs(
      title = "Distribution of RD Treatment Effect Estimates",
      subtitle = paste0("True effect = ", round(true_effect, 2), 
                        ", Mean = ", round(mean_estimate, 2), 
                        " (SD = ", round(sd_estimate, 2), ")"),
      x = "Estimated Treatment Effect",
      y = "Density"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 10)
    )
}

# 3. Error metrics
plot_error_metrics <- function(results) {
  # Calculate bias, MSE, MAE
  results <- results %>%
    mutate(
      bias = rd_estimate - true_effect,
      sq_error = (rd_estimate - true_effect)^2,
      abs_error = abs(rd_estimate - true_effect)
    )
  
  # Calculate summary statistics
  bias <- mean(results$bias)
  mse <- mean(results$sq_error)
  mae <- mean(results$abs_error)
  sd_val <- sd(results$rd_estimate)
  
  # Create data frame for plotting
  error_summary <- data.frame(
    Metric = factor(c("Bias", "MSE", "MAE", "SD"), 
                    levels = c("Bias", "MSE", "MAE", "SD")),
    Value = c(bias, mse, mae, sd_val)
  )
  
  # Create bar plot
  ggplot(error_summary, aes(x = Metric, y = Value, fill = Metric)) +
    geom_col() +
    geom_text(aes(label = round(Value, 3)), vjust = -0.5, size = 4) +
    scale_fill_manual(values = rep("darkblue", 4)) +
    labs(
      title = "RD Estimator Performance Metrics",
      x = "Metric",
      y = "Value"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      legend.position = "none"
    )
}

# 4. Performance across varying concentration parameters (kappa)
plot_kappa_performance <- function(results_varying_kappa) {
  # Calculate bias and MSE for each kappa
  summary_by_kappa <- results_varying_kappa %>%
    group_by(kappa) %>%
    summarize(
      bias = mean(rd_estimate - true_effect),
      abs_bias = mean(abs(rd_estimate - true_effect)),
      mse = mean((rd_estimate - true_effect)^2),
      sd = sd(rd_estimate)
    )
  
  # Plot bias by kappa
  p1 <- ggplot(summary_by_kappa, aes(x = kappa, y = abs_bias)) +
    geom_line(size = 1, color = "darkblue") +
    geom_point(size = 3, color = "darkblue") +
    scale_x_log10() +  # Log scale for kappa
    labs(
      title = "Absolute Bias by Concentration Parameter",
      x = "Concentration Parameter (κ, log scale)",
      y = "Absolute Bias"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold")
    )
  
  # Plot MSE by kappa
  p2 <- ggplot(summary_by_kappa, aes(x = kappa, y = mse)) +
    geom_line(size = 1, color = "darkblue") +
    geom_point(size = 3, color = "darkblue") +
    scale_x_log10() +  # Log scale for kappa
    labs(
      title = "Mean Squared Error by Concentration Parameter",
      x = "Concentration Parameter (κ, log scale)",
      y = "MSE"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold")
    )
  
  # Plot SD by kappa
  p3 <- ggplot(summary_by_kappa, aes(x = kappa, y = sd)) +
    geom_line(size = 1, color = "darkblue") +
    geom_point(size = 3, color = "darkblue") +
    scale_x_log10() +  # Log scale for kappa
    labs(
      title = "Standard Deviation by Concentration Parameter",
      x = "Concentration Parameter (κ, log scale)",
      y = "Standard Deviation"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold")
    )
  
  return(list(p1 = p1, p2 = p2, p3 = p3))
}

# Generate all plots
p_sample <- visualize_sample_data()
p_dist <- plot_effect_distribution(results)
p_metrics <- plot_error_metrics(results)
p_kappa <- plot_kappa_performance(results_varying_kappa)

# Print summary statistics
cat("\n=== Canonical RD Simulation Results ===\n\n")
cat(sprintf("True treatment effect: %.4f\n", true_effect))
cat(sprintf("Number of simulations: %d\n", n_sims))
cat(sprintf("Mean estimate: %.4f\n", mean(results$rd_estimate)))
cat(sprintf("Standard deviation: %.4f\n", sd(results$rd_estimate)))
cat(sprintf("Bias: %.4f\n", mean(results$rd_estimate) - true_effect))
cat(sprintf("MSE: %.4f\n", mean((results$rd_estimate - true_effect)^2)))
cat(sprintf("MAE: %.4f\n", mean(abs(results$rd_estimate - true_effect))))

# Display plots
p1 <- gridExtra::grid.arrange(p_sample, p_dist, ncol = 2)

p2 <- gridExtra::grid.arrange(p_metrics,
                              p_kappa$p1,
                              p_kappa$p2,
                              p_kappa$p3,
                              nrow = 2, 
                              ncol = 2)

ggsave("fig/rd_sim.png", p1, width = 12, height = 6, dpi = 300)
ggsave("fig/rd_diagnostics.png", p2, width = 12, height = 6, dpi = 300)
