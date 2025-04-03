# ==============================================================================
# Title: Projection of Points
# Original Author: Daniel Posmik
# E-Mail: daniel_posmik@brown.edu
# ==============================================================================

# Set-Up =======================================================================
library(dplyr)
library(tidyverse)
library(circular)
library(plotly)

# Set working directory
setwd("~/Documents/brown/classes/spring25/php2530-bayes/directional-priors")

# Import data
df <- read.csv("data/01_sim-results.csv")

# Projecting data on unit sphere ===============================================
# Function to compute distance to unit sphere
sphere_dist <- function(p, noisy_point) {
  # p is normalized during optimization
  p_norm <- p / sqrt(sum(p^2))
  # Compute Euclidean distance
  return(sqrt(sum((noisy_point - p_norm)^2)))
}

# Project points to unit sphere
project_to_sphere <- function(df) {
  results <- data.frame(
    X_proj = numeric(nrow(df)),
    Y_proj = numeric(nrow(df)),
    Z_proj = numeric(nrow(df)),
    distance = numeric(nrow(df)),
    theta = numeric(nrow(df)),
    phi = numeric(nrow(df))
  )
  
  for (i in 1:nrow(df)) {
    noisy_point <- c(df$X_noisy[i], df$Y_noisy[i], df$Z_noisy[i])
    
    # Use noisy point as initial guess
    opt <- optim(
      par = noisy_point,
      fn = sphere_dist,
      noisy_point = noisy_point,
      method = "BFGS"
    )
    
    # Normalize to ensure exactly on sphere
    proj_point <- opt$par / sqrt(sum(opt$par^2))
    results$X_proj[i] <- proj_point[1]
    results$Y_proj[i] <- proj_point[2]
    results$Z_proj[i] <- proj_point[3]
    
    # Calculate distance
    results$distance[i] <- sqrt(sum((noisy_point - proj_point)^2))
    
    # Calculate projection vector (from noisy to projected)
    proj_vector <- proj_point - noisy_point
    
    # Convert to spherical coordinates (in local frame centered at noisy point)
    # r = distance (already calculated)
    # theta = azimuthal angle in xy-plane from x-axis (0 to 2π)
    # phi = polar angle from z-axis (0 to π)
    
    # Calculate theta (azimuthal angle in xy-plane)
    results$theta[i] <- atan2(proj_vector[2], proj_vector[1])
    # Ensure theta is in [0, 2π)
    if (results$theta[i] < 0) {
      results$theta[i] <- results$theta[i] + 2*pi
    }
    
    # Calculate phi (polar angle from z-axis)
    results$phi[i] <- acos(proj_vector[3] / results$distance[i])
  }
  
  return(results)
}

# Apply projection function to data
projections <- project_to_sphere(df)
projections$kappa <- 1

# Combine with original data
df_proj <- bind_cols(df, projections)

# Add global spherical coordinates
df_proj <- df_proj %>%
  mutate(
    # Calculate theta_global (azimuthal angle in xy-plane)
    theta_global = atan2(Y_proj, X_proj),
    # Ensure theta_global is in [0, 2π)
    theta_global = ifelse(theta_global < 0, theta_global + 2*pi, theta_global),
    # Calculate phi_global (polar angle from z-axis)
    phi_global = acos(Z_proj)
  )

# Save data 
write.csv(df_proj, file = "data/02_proj-results.csv")

# Visualization ================================================================
# Define inputs
theta <- seq(0, 2*pi, length.out = 30)
phi <- seq(0, pi, length.out = 30)
theta_mat <- outer(theta, rep(1, length(phi)))
phi_mat <- outer(rep(1, length(theta)), phi)
x_sphere <- cos(theta_mat) * sin(phi_mat)
y_sphere <- sin(theta_mat) * sin(phi_mat)
z_sphere <- cos(phi_mat)

# Show projections
plot_ly() %>%
  # Add unit sphere (surface)
  add_surface(
    x = x_sphere,
    y = y_sphere,
    z = z_sphere,
    opacity = 0.3,
    colorscale = list(c(0, 1), c("lightblue", "lightblue")),
    showscale = FALSE,
    name = "Unit Sphere"
  ) %>%
  # Add noisy points
  add_markers(
    x = df_proj$X_noisy, 
    y = df_proj$Y_noisy, 
    z = df_proj$Z_noisy,
    color = I("red"),
    size = 2,
    name = "Noisy Points"
  ) %>%
  # Add projected points
  add_markers(
    x = df_proj$X_proj, 
    y = df_proj$Y_proj, 
    z = df_proj$Z_proj,
    color = I("blue"),
    size = 2,
    name = "Projected Points"
  ) %>%
  # Add lines connecting original to projected
  add_trace(
    type = "scatter3d",
    mode = "lines",
    x = unlist(lapply(1:nrow(df_proj), function(i) c(df_proj$X_noisy[i], df_proj$X_proj[i], NA))),
    y = unlist(lapply(1:nrow(df_proj), function(i) c(df_proj$Y_noisy[i], df_proj$Y_proj[i], NA))),
    z = unlist(lapply(1:nrow(df_proj), function(i) c(df_proj$Z_noisy[i], df_proj$Z_proj[i], NA))),
    line = list(color = "gray", width = 1),
    showlegend = FALSE
  ) %>%
  layout(
    title = "Projection of Noisy Points onto Unit Sphere",
    scene = list(
      xaxis = list(title = "X"),
      yaxis = list(title = "Y"),
      zaxis = list(title = "Z"),
      aspectmode = "cube"
    )
  )

# Show origin vectors
plot_ly() %>%
  # Add unit sphere (surface)
  add_surface(
    x = x_sphere,
    y = y_sphere,
    z = z_sphere,
    opacity = 0.2,
    colorscale = list(c(0, 1), c("lightblue", "lightblue")),
    showscale = FALSE,
    name = "Unit Sphere"
  ) %>%
  # Add origin point
  add_markers(
    x = 0, 
    y = 0, 
    z = 0,
    color = I("black"),
    size = 4,
    name = "Origin"
  ) %>%
  # Add projected points
  add_markers(
    x = df_proj$X_proj, 
    y = df_proj$Y_proj, 
    z = df_proj$Z_proj,
    color = I("blue"),
    size = 3,
    name = "Projected Points"
  ) %>%
  # Add vectors from origin to projected points
  add_trace(
    type = "scatter3d",
    mode = "lines",
    x = unlist(lapply(1:nrow(df_proj), function(i) c(0, df_proj$X_proj[i], NA))),
    y = unlist(lapply(1:nrow(df_proj), function(i) c(0, df_proj$Y_proj[i], NA))),
    z = unlist(lapply(1:nrow(df_proj), function(i) c(0, df_proj$Z_proj[i], NA))),
    line = list(color = "green", width = 2),
    name = "Vectors from Origin"
  ) %>%
  # Add axes
  add_trace(
    type = "scatter3d",
    mode = "lines",
    x = c(0, 1.5, NA, 0, NA, 0),
    y = c(0, 0, NA, 0, 1.5, NA),
    z = c(0, 0, NA, 0, 0, 1.5),
    line = list(color = "gray", width = 1, dash = "dash"),
    name = "Coordinate Axes"
  ) %>%
  # Add labels for axes
  add_trace(
    type = "scatter3d",
    mode = "text",
    x = c(1.6, 0, 0),
    y = c(0, 1.6, 0),
    z = c(0, 0, 1.6),
    text = c("X", "Y", "Z"),
    textfont = list(size = 12),
    showlegend = FALSE
  ) %>%
  layout(
    title = "Vectors from Origin to Projected Points on Unit Sphere",
    scene = list(
      xaxis = list(title = "", range = c(-1.5, 1.5)),
      yaxis = list(title = "", range = c(-1.5, 1.5)),
      zaxis = list(title = "", range = c(-1.5, 1.5)),
      aspectmode = "cube"
    )
  )
