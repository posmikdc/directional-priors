# ==============================================================================
# Title: Bayesian Inference with Spherical Normal Distribution
# Original Author: Daniel Posmik
# E-Mail: daniel_posmik@brown.edu
# ==============================================================================

# Set-Up =======================================================================
library(dplyr)
library(tidyverse)
library(circular)
library(plotly)
library(mvtnorm)

# Set working directory
setwd("~/Documents/brown/classes/spring25/php2530-bayes/directional-priors")

# Import data
df_proj <- read.csv("data/02_proj-results.csv")

# Posterior inference ==========================================================
# Function to calculate posterior quantities
calculate_posterior <- function(prior_theta, prior_kappa = 1) {
  prior_var <- 1/(prior_kappa)^2
  data_kappa <- mean(df_proj$kappa)
  data_var <- 1/(data_kappa)^2
  data_theta <- mean(df_proj$theta_global) + pi
  
  # Posterior estimates
  post_theta <- (prior_var)/(prior_var + data_var) * prior_theta + 
    (data_var)/(prior_var + data_var) * data_theta + pi
  post_var <- 1/(1/prior_var + 1/data_var)
  
  # Ensure post_theta is within [0, 2π]
  post_theta <- post_theta %% (2*pi)
  
  return(list(
    post_theta = post_theta,
    post_var = post_var
  ))
}

# Prior means to visualize
prior_means <- c(pi/2, pi, 3*pi/2)
prior_results <- lapply(prior_means, calculate_posterior)

# Create Rose Plot Style Visualization =========================================
create_rose_plot_style <- function(prior_theta, posterior_result, data_angles, show_legend = FALSE) {
  # Number of bins for the rose plot
  n_bins <- 24
  bin_width <- 2*pi/n_bins
  
  # Create bins for data
  breaks <- seq(0, 2*pi, length.out = n_bins + 1)
  bin_mids <- breaks[-1] - bin_width/2
  
  # Count data in each bin
  data_bins <- table(cut(data_angles %% (2*pi), breaks = breaks, include.lowest = TRUE))
  data_bins <- as.numeric(data_bins)
  
  # Generate wrapped normal samples for posterior
  n_samples <- 10000
  wrapped_samples <- rnorm(n_samples, 
                           mean = posterior_result$post_theta,
                           sd = sqrt(posterior_result$post_var)) %% (2*pi)
  
  # Count posterior samples in each bin
  post_bins <- table(cut(wrapped_samples, breaks = breaks, include.lowest = TRUE))
  post_bins <- as.numeric(post_bins)
  
  # Normalize bin counts
  data_bins_norm <- data_bins / max(data_bins)
  post_bins_norm <- post_bins / max(post_bins)
  
  # Create data for the rose plot
  plot_data <- data.frame(
    bin_mid = bin_mids,
    bin_start = breaks[-length(breaks)],
    bin_end = breaks[-1],
    data_count = data_bins_norm,
    post_count = post_bins_norm
  )
  
  # Create a function to generate rose plot wedge paths
  create_wedge_path <- function(start_angle, end_angle, radius) {
    theta_seq <- seq(start_angle, end_angle, length.out = 10)
    
    # Create the path from center to outer edge to center
    path_x <- c(0, radius * cos(theta_seq), 0)
    path_y <- c(0, radius * sin(theta_seq), 0)
    
    return(list(x = path_x, y = path_y))
  }
  
  # Calculate data mean angle
  data_mean <- circular::mean.circular(circular(data_angles, units="radians"))
  
  # Initialize plotly
  p <- plot_ly()
  
  # Add reference circles
  for (r in seq(0.2, 1, by = 0.2)) {
    circle_x <- r * cos(seq(0, 2*pi, length.out = 100))
    circle_y <- r * sin(seq(0, 2*pi, length.out = 100))
    
    p <- p %>% add_trace(
      x = circle_x, y = circle_y,
      type = 'scatter', mode = 'lines',
      line = list(color = 'gray', width = 0.5),
      showlegend = FALSE
    )
  }
  
  # Add data wedges (now in light red)
  for (i in 1:nrow(plot_data)) {
    if (plot_data$data_count[i] > 0) {
      wedge <- create_wedge_path(
        plot_data$bin_start[i], 
        plot_data$bin_end[i],
        plot_data$data_count[i]
      )
      
      p <- p %>% add_trace(
        x = wedge$x, y = wedge$y,
        fill = 'toself',
        fillcolor = 'rgba(255, 200, 200, 0.7)',  # Light red
        line = list(color = 'darkred', width = 0.5),
        showlegend = show_legend && i == 1,
        name = 'Data Distribution'
      )
    }
  }
  
  # Add posterior wedges (now in light blue)
  for (i in 1:nrow(plot_data)) {
    if (plot_data$post_count[i] > 0) {
      wedge <- create_wedge_path(
        plot_data$bin_start[i], 
        plot_data$bin_end[i],
        plot_data$post_count[i]
      )
      
      p <- p %>% add_trace(
        x = wedge$x, y = wedge$y,
        fill = 'toself',
        fillcolor = 'rgba(200, 200, 255, 0.7)',  # Light blue
        line = list(color = 'darkblue', width = 0.5),
        showlegend = show_legend && i == 1,
        name = 'Posterior Distribution'
      )
    }
  }
  
  # Add angle labels around the circle
  angle_labels <- seq(0, 330, by = 30)
  label_radius <- 1.2
  annotations <- list()
  
  for (angle_deg in angle_labels) {
    angle_rad <- angle_deg * pi / 180
    annotations <- c(annotations, list(
      list(
        x = label_radius * cos(angle_rad),
        y = label_radius * sin(angle_rad),
        text = as.character(angle_deg),
        showarrow = FALSE,
        font = list(size = 10)
      )
    ))
  }
  
  # Add prior mean
  p <- p %>% add_trace(
    x = c(0, 1.1 * cos(prior_theta)),
    y = c(0, 1.1 * sin(prior_theta)),
    type = 'scatter', mode = 'lines',
    line = list(color = 'green', width = 2, dash = 'dash'),
    showlegend = show_legend,
    name = 'Prior Mean'
  )
  
  # Add data mean
  p <- p %>% add_trace(
    x = c(0, 1.1 * cos(data_mean)),
    y = c(0, 1.1 * sin(data_mean)),
    type = 'scatter', mode = 'lines',
    line = list(color = 'darkred', width = 2, dash = 'dash'),
    showlegend = show_legend,
    name = 'Data Mean'
  )
  
  # Add posterior mean
  p <- p %>% add_trace(
    x = c(0, 1.1 * cos(posterior_result$post_theta)),
    y = c(0, 1.1 * sin(posterior_result$post_theta)),
    type = 'scatter', mode = 'lines',
    line = list(color = 'darkblue', width = 2, dash = 'dash'),
    showlegend = show_legend,
    name = 'Posterior Mean'
  )
  
  # Layout adjustments
  p <- p %>% layout(
    title = paste("Prior Mean =", round(prior_theta * 180/pi, 1), "°"),
    xaxis = list(
      title = "",
      zeroline = FALSE,
      showgrid = FALSE,
      showticklabels = FALSE,
      range = c(-1.3, 1.3)
    ),
    yaxis = list(
      title = "",
      zeroline = FALSE,
      showgrid = FALSE,
      showticklabels = FALSE,
      range = c(-1.3, 1.3),
      scaleanchor = "x"
    ),
    annotations = annotations,
    showlegend = FALSE,
    margin = list(l = 40, r = 40, b = 40, t = 60)
  )
  
  return(p)
}

# Create visualizations ========================================================
# Prepare the data
data_angles <- df_proj$theta_global

# Generate the plots - only show legend for plot1
plot1 <- create_rose_plot_style(prior_means[1], prior_results[[1]], data_angles, show_legend = TRUE)
plot2 <- create_rose_plot_style(prior_means[2], prior_results[[2]], data_angles, show_legend = FALSE)
#plot3 <- create_rose_plot_style(prior_means[3], prior_results[[3]], data_angles, show_legend = FALSE)

# Create the combined plot with shared legend
combined_plot <- plotly::subplot(plot1, plot2, nrows = 1, shareY = TRUE) %>%
  layout(
    showlegend = TRUE,
    legend = list(
      orientation = "h",      # Horizontal legend
      xanchor = "center",     # Center the legend horizontally
      x = 0.5,                # Position at the center of the plot
      y = 0                   # Position at the bottom
    ),
    margin = list(b = 100)    # Add bottom margin for the legend
  )

# Display the combined plot
combined_plot

ggsave("fig/posterior-roseplot.png", combined_plot, width = 12, height = 6, dpi = 300)
