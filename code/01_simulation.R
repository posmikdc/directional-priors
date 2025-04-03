# ==============================================================================
# Title: Simulation of Data
# Author: Daniel Posmik
# E-Mail: daniel_posmik@brown.edu
# ==============================================================================

# Set-Up =======================================================================
library(dplyr)
library(tidyverse)
library(circular)
library(plotly)

# Set working directory
setwd("~/Documents/brown/classes/spring25/php2530-bayes/directional-priors")

# Set random seed for reproducibility
set.seed(1)

# Generate True Data ===========================================================
# Sample size
n <- 500

# Sample from circular normal distribution (wrapped normal)
theta <- rwrappednormal(n, mu = circular(0), rho = 0.8)
phi <- rwrappednormal(n, mu = circular(pi/2), rho = 0.8)

# Convert spherical coordinates to Cartesian
X <- sin(phi) * cos(theta)
Y <- sin(phi) * sin(theta)
Z <- cos(phi)

# Create data frame
df <- data.frame(X, Y, Z)

# Simulate Data ================================================================
# Noise parameters 
E <- 0.25

# Add noisy versions of each point
# X_noisy, Y_noisy, Z_noisy are sampled from N(original point, E)
df$X_noisy <- df$X + rnorm(n, mean = 0, sd = E)
df$Y_noisy <- df$Y + rnorm(n, mean = 0, sd = E)
df$Z_noisy <- df$Z + rnorm(n, mean = 0, sd = E)

# Save data
write.csv(df, file = "data/01_sim-results.csv")

# Visualization on unit sphere
sphere_points <- 20
theta_sphere <- seq(0, 2*pi, length.out = sphere_points)
phi_sphere <- seq(0, pi, length.out = sphere_points)
theta_matrix <- matrix(rep(theta_sphere, sphere_points), sphere_points, sphere_points)
phi_matrix <- matrix(rep(phi_sphere, each = sphere_points), sphere_points, sphere_points)
x_sphere <- sin(phi_matrix) * cos(theta_matrix)
y_sphere <- sin(phi_matrix) * sin(theta_matrix)
z_sphere <- cos(phi_matrix)

# 3D plotly visualization
fig <- plot_ly() %>%
  # Add unit sphere
  add_surface(x = x_sphere, y = y_sphere, z = z_sphere, 
              opacity = 0.2, colorscale = list(c(0, 1), c("lightblue", "lightblue")),
              showscale = FALSE) %>%
  # Add original data points
  add_markers(data = df, x = ~X, y = ~Y, z = ~Z, 
              color = I("red"), size = 6, 
              name = "Original Points") %>%
  # Add noisy data points
  add_markers(data = df, x = ~X_noisy, y = ~Y_noisy, z = ~Z_noisy, 
              color = I("blue"), size = 6, opacity = 0.7,
              name = "Noisy Points") %>%
  # Layout
  layout(title = paste("3D Circular Normal with Noise (E =", E, ")"),
         scene = list(
           aspectmode = "cube",
           xaxis = list(title = "X"),
           yaxis = list(title = "Y"),
           zaxis = list(title = "Z")
         ))

# Display the plot
fig







