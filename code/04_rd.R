# ==============================================================================
# Title: Manifold RD Simulation
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
# Simulation ===================================================================

# Parameters
n <- 500
mean_post <- pi
mean_pre <- pi/2
var_post <- var_pre <- 1
kappa_post <- 0.25
kappa_pre <- 1  

# Generate samples
set.seed(2025)
# Pre-treatment: mean at pi/2, variance 1
pre_treatment <- rvonmises(n, mu = circular(mean_pre), kappa = kappa_pre)
# Post-treatment: mean at pi, variance 1
post_treatment <- rvonmises(n, mu = circular(mean_post), kappa = kappa_post)

# Calculate treatment effect distribution parameters
mean_effect <- mean_post - mean_pre
# Variance formula as specified
variance_effect <- sqrt((kappa_pre * cos(mean_pre) + kappa_post * cos(mean_post))^2 + 
                          (kappa_pre * sin(mean_pre) + kappa_post * sin(mean_post))^2)

# Calculate treatment effect samples
# For circular data, we compute the difference between post and pre angles
# We need to wrap the angles to ensure they stay in [0, 2π]
treatment_effect <- (as.numeric(post_treatment) - as.numeric(pre_treatment)) %% (2*pi)

# Function to create rose plot with plotly
create_rose_plot_style <- function(angles, mean_angle, title, kappa_or_var = 1, is_variance = FALSE, 
                                   plot_color = 'rgba(200, 200, 255, 0.7)', line_color = 'darkblue') {
  # Number of bins for the rose plot
  n_bins <- 24
  bin_width <- 2*pi/n_bins
  
  # Create bins for data
  breaks <- seq(0, 2*pi, length.out = n_bins + 1)
  bin_mids <- breaks[-1] - bin_width/2
  
  # Count data in each bin
  data_bins <- table(cut(angles %% (2*pi), breaks = breaks, include.lowest = TRUE))
  data_bins <- as.numeric(data_bins)
  
  # Normalize bin counts
  data_bins_norm <- data_bins / max(data_bins)
  
  # Create data for the rose plot
  plot_data <- data.frame(
    bin_mid = bin_mids,
    bin_start = breaks[-length(breaks)],
    bin_end = breaks[-1],
    data_count = data_bins_norm
  )
  
  # Create a function to generate rose plot wedge paths
  create_wedge_path <- function(start_angle, end_angle, radius) {
    theta_seq <- seq(start_angle, end_angle, length.out = 10)
    
    # Create the path from center to outer edge to center
    path_x <- c(0, radius * cos(theta_seq), 0)
    path_y <- c(0, radius * sin(theta_seq), 0)
    
    return(list(x = path_x, y = path_y))
  }
  
  # Calculate data mean angle if not provided
  if (is.null(mean_angle)) {
    mean_angle <- circular::mean.circular(circular(angles, units="radians"))
  }
  
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
  
  # Add data wedges
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
        fillcolor = plot_color,
        line = list(color = line_color, width = 0.5),
        showlegend = i == 1,  # Only show legend for the first wedge
        name = title
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
  
  # Add mean direction arrow
  p <- p %>% add_trace(
    x = c(0, 1.1 * cos(mean_angle)),
    y = c(0, 1.1 * sin(mean_angle)),
    type = 'scatter', mode = 'lines',
    line = list(color = line_color, width = 2, dash = 'dash'),
    showlegend = FALSE,
    name = 'Mean'
  )
  
  # Add parameter annotations at the top
  mean_text <- paste("μ =", round(mean_angle * 180/pi, 1), "°")
  param_text <- ifelse(is_variance, 
                       paste("σ =", round(kappa_or_var, 2)),
                       paste("κ =", round(kappa_or_var, 2)))
  
  # Add parameter annotations at the top
  annotations <- c(annotations, list(
    list(
      x = 0, y = -1.65,
      text = paste(mean_text, ", ", param_text),
      showarrow = FALSE,
      xanchor = "center",
      font = list(size = 11)
    )
  ))
  
  # Layout adjustments
  p <- p %>% layout(
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
    showlegend = TRUE,
    margin = list(l = 40, r = 40, b = 40, t = 40) # Added top margin for parameter text
  )
  
  return(p)
}

# Create rose plots with different colors
plot1 <- create_rose_plot_style(
  as.numeric(pre_treatment), 
  mean_pre, 
  "Pre-treatment Distribution", 
  kappa_pre,
  plot_color = 'rgba(173, 216, 230, 0.7)',  # Light blue
  line_color = 'navy'
)

plot2 <- create_rose_plot_style(
  as.numeric(post_treatment), 
  mean_post, 
  "Post-treatment Distribution", 
  kappa_post,
  plot_color = 'rgba(144, 238, 144, 0.7)',  # Light green
  line_color = 'darkgreen'
)

plot3 <- create_rose_plot_style(
  treatment_effect, 
  mean_effect, 
  "Treatment Effect Distribution", 
  variance_effect, 
  is_variance = TRUE,
  plot_color = 'rgba(255, 222, 173, 0.7)',  # Light orange
  line_color = 'darkorange'
)

# Create the combined plot with shared legend at the bottom
combined_plot <- plotly::subplot(plot1, plot2, plot3, nrows = 1, shareY = TRUE) %>%
  layout(
    showlegend = TRUE,
    legend = list(
      orientation = "h",      # Horizontal legend
      xanchor = "center",     # Center the legend horizontally
      x = 0.5,                # Position at the center of the plot
      y = 0                   # Position at the bottom
    ),
    margin = list(b = 80)     # Add bottom margin for the legend
  )

# Display the combined plot
combined_plot

# Create rose plots with different colors
plot1 <- create_rose_plot_style(
  as.numeric(pre_treatment), 
  mean_pre, 
  "Pre-treatment Distribution", 
  kappa_pre,
  plot_color = 'rgba(173, 216, 230, 0.7)',  # Light blue
  line_color = 'navy'
)

plot2 <- create_rose_plot_style(
  as.numeric(post_treatment), 
  mean_post, 
  "Post-treatment Distribution", 
  kappa_post,
  plot_color = 'rgba(144, 238, 144, 0.7)',  # Light green
  line_color = 'darkgreen'
)

plot3 <- create_rose_plot_style(
  treatment_effect, 
  mean_effect, 
  "Treatment Effect Distribution", 
  variance_effect, 
  is_variance = TRUE,
  plot_color = 'rgba(255, 222, 173, 0.7)',  # Light orange
  line_color = 'darkorange'
)

# Create the combined plot with shared legend at the bottom
combined_plot <- plotly::subplot(plot1, plot2, plot3, nrows = 1, shareY = TRUE) %>%
  layout(
    showlegend = TRUE,
    legend = list(
      orientation = "h",      # Horizontal legend
      xanchor = "center",     # Center the legend horizontally
      x = 0.5,                # Position at the center of the plot
      y = 0.2                   # Position at the bottom
    ),
    margin = list(b = 80)     # Add bottom margin for the legend
  )

# Display the combined plot
combined_plot