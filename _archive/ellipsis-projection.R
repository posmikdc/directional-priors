#===============================================================================
# Title: Prior Information in Manifold Projections
# Author: Daniel Posmik
# E-Mail: daniel_posmik@brown.edu
#===============================================================================

# Set-Up =======================================================================
library(dplyr)
library(plotly)

# Simulation ===================================================================
# Set random seed for reproducibility
set.seed(123)

# Parameters for the elliptic curve
a <- 3  # semi-major axis
b <- 1  # semi-minor axis

# Generate 5 points on an elliptic curve
n <- 5
theta <- seq(0, 2*pi, length.out = n+1)[1:n]  # equally spaced angles

# Generate points on the perfect ellipse
x_ellipse <- a * cos(theta)
y_ellipse <- b * sin(theta)

# Add Gaussian noise N(0,1) to create observed points
noise_sd <- 0.2
x_observed <- x_ellipse + rnorm(n, 0, noise_sd)
y_observed <- y_ellipse + rnorm(n, 0, noise_sd)

# Generate uniform(0,1) values for the z coordinate
z_values <- runif(n)

# Create the data frame
data <- data.frame(
  x = x_observed,
  y = y_observed,
  z = z_values
)

# Projection ===================================================================
# Project a point onto the ellipse and return coordinates and angles
project_to_ellipse <- function(x, y, z, a, b) {
  # Starting guess for theta (x-y plane angle)
  theta_guess <- atan2(y/b, x/a)
  
  # Function to minimize: squared distance from (x,y,z) to a point on the ellipse
  distance_function <- function(theta) {
    ellipse_x <- a * cos(theta)
    ellipse_y <- b * sin(theta)
    ellipse_z <- 0 # Since our ellipse only lives in x,y space
    return(sqrt((x - ellipse_x)^2 + (y - ellipse_y)^2 + (z - ellipse_z)^2))
  }
  
  # Optimize to find the best theta
  result <- optimize(distance_function, c(0, 2*pi))
  theta_optimal <- result$minimum
  
  # Calculate projected coordinates
  x_proj <- a * cos(theta_optimal)
  y_proj <- b * sin(theta_optimal)
  z_proj <- 0 # The ellipse lives exclusively in the x-y plane
  
  # Calculate distance between original and projected point
  dist <- sqrt((x - x_proj)^2 + (y - y_proj)^2 + (z - z_proj)^2)
  
  # Calculate the second angle (phi) - elevation from x-y plane
  # This is the angle between the x-y plane and the line to the original point
  phi <- atan2(z, sqrt(x^2 + y^2))
  
  # Calculate theta_centered - the angle from point's projection to ellipse point
  # This is the angle between (x,y,0) to (x_proj,y_proj,0) centered at (x,y,0)
  # Using atan2 to get the standard unit circle angle orientation
  theta_centered <- atan2(y_proj - y, x_proj - x)
  
  # Convert to standardized angle where 0° is along positive x-axis and angles increase counterclockwise
  # Ensure the angle is in [0, 360) range
  theta_centered_deg <- (theta_centered * 180/pi) %% 360
  
  # Calculate the derivative (tangent vector) at the projection point
  # For an ellipse parametrized as (a*cos(t), b*sin(t)), the derivative is (-a*sin(t), b*cos(t))
  dx_dt <- -a * sin(theta_optimal)
  dy_dt <- b * cos(theta_optimal)
  
  # Calculate the magnitude of the derivative vector
  derivative_magnitude <- sqrt(dx_dt^2 + dy_dt^2)
  
  # Normalize the derivative to get a unit tangent vector
  dx_dt_normalized <- dx_dt / derivative_magnitude
  dy_dt_normalized <- dy_dt / derivative_magnitude
  
  # Also calculate the normal vector (perpendicular to the tangent)
  # For a 2D curve, the normal is simply the tangent rotated by 90 degrees
  nx <- -dy_dt_normalized
  ny <- dx_dt_normalized
  
  # Return the projected point coordinates, angles, distances, and derivatives
  return(list(
    x_proj = x_proj,
    y_proj = y_proj,
    z_proj = z_proj,
    theta = (theta_optimal * 180/pi) %% 360,  # Azimuthal angle in x-y plane (in degrees, standardized to [0,360))
    theta_centered = theta_centered_deg, # Centered angle (in degrees, standardized to [0,360))
    phi = phi * 180/pi,              # Elevation angle from x-y plane (in degrees)
    distance = dist,
    theta_param = theta_optimal,         # The parameter value (in radians)
    dx_dt = dx_dt,                       # x-component of derivative
    dy_dt = dy_dt,                       # y-component of derivative
    derivative_magnitude = derivative_magnitude, # Scalar magnitude of the derivative
    dx_dt_norm = dx_dt_normalized,       # x-component of normalized derivative (unit tangent)
    dy_dt_norm = dy_dt_normalized,       # y-component of normalized derivative (unit tangent)
    nx = nx,                             # x-component of normal vector
    ny = ny                              # y-component of normal vector
  ))
}

# Apply the projection to all points
results <- list()
for (i in 1:nrow(data)) {
  results[[i]] <- project_to_ellipse(data$x[i], data$y[i], data$z[i], a, b)
}

# Create a data frame with all results
projections <- data.frame(
  x_orig = data$x,
  y_orig = data$y,
  z_orig = data$z,
  x_proj = sapply(results, function(r) r$x_proj),
  y_proj = sapply(results, function(r) r$y_proj),
  z_proj = sapply(results, function(r) r$z_proj),
  theta = sapply(results, function(r) r$theta),
  theta_centered = sapply(results, function(r) r$theta_centered),
  phi = sapply(results, function(r) r$phi),
  distance = sapply(results, function(r) r$distance),
  theta_param = sapply(results, function(r) r$theta_param),
  dx_dt = sapply(results, function(r) r$dx_dt),
  dy_dt = sapply(results, function(r) r$dy_dt),
  derivative_magnitude = sapply(results, function(r) r$derivative_magnitude),
  dx_dt_norm = sapply(results, function(r) r$dx_dt_norm),
  dy_dt_norm = sapply(results, function(r) r$dy_dt_norm),
  nx = sapply(results, function(r) r$nx),
  ny = sapply(results, function(r) r$ny)
)

# Print the results
head(projections)

# Print specifically the derivative magnitudes
cat("\nDerivative Magnitudes at Projection Points:\n")
for (i in 1:nrow(projections)) {
  cat(sprintf("Point %d: (%.2f, %.2f) - Derivative Magnitude: %.4f\n", 
              i, 
              projections$x_proj[i], 
              projections$y_proj[i], 
              projections$derivative_magnitude[i]))
}

# Visualization ================================================================
# Create a dense set of points for the perfect ellipse for visualization
ellipse_theta <- seq(0, 2*pi, length.out = 100)
ellipse_x <- a * cos(ellipse_theta)
ellipse_y <- b * sin(ellipse_theta)
ellipse_z <- rep(0, 100)  # The ellipse is in the x-y plane

# Create the 3D visualization
fig <- plot_ly()

# Add the perfect ellipse
fig <- fig %>% add_trace(
  x = ellipse_x, y = ellipse_y, z = ellipse_z,
  type = 'scatter3d', mode = 'lines',
  line = list(color = 'black', width = 3),
  name = 'Perfect Ellipse'
)

# Add the original points
fig <- fig %>% add_trace(
  x = projections$x_orig, y = projections$y_orig, z = projections$z_orig,
  type = 'scatter3d', mode = 'markers',
  marker = list(color = 'red', size = 8),
  name = 'Original Points'
)

# Add the projected points on the ellipse
fig <- fig %>% add_trace(
  x = projections$x_proj, y = projections$y_proj, z = projections$z_proj,
  type = 'scatter3d', mode = 'markers',
  marker = list(color = 'blue', size = 8),
  name = 'Projected Points on Ellipse'
)

# Add the projection of the original points onto the x-y plane
fig <- fig %>% add_trace(
  x = projections$x_orig, y = projections$y_orig, z = rep(0, nrow(projections)),
  type = 'scatter3d', mode = 'markers',
  marker = list(color = 'green', size = 8),
  name = 'Original Points Projected to x-y Plane'
)

# Add the triangles that show the projection and theta_centered angle
for (i in 1:nrow(projections)) {
  # Create the triangle: original point, projected point on ellipse, and projection of original to x-y plane
  x_vals <- c(projections$x_orig[i], projections$x_proj[i], projections$x_orig[i], projections$x_orig[i])
  y_vals <- c(projections$y_orig[i], projections$y_proj[i], projections$y_orig[i], projections$y_orig[i])
  z_vals <- c(projections$z_orig[i], projections$z_proj[i], 0, projections$z_orig[i])
  
  # Add the triangle as a 3D line
  fig <- fig %>% add_trace(
    x = x_vals, y = y_vals, z = z_vals,
    type = 'scatter3d', mode = 'lines',
    line = list(color = 'purple', width = 2),
    name = paste('Triangle', i),
    showlegend = (i == 1)  # Only show one legend entry for all triangles
  )
  
  # Add a line showing the theta_centered angle in the x-y plane
  # We'll create an arc to represent the angle
  
  # Calculate the radius for the arc (1/4 of the distance to the projected point for visibility)
  r <- sqrt((projections$x_proj[i] - projections$x_orig[i])^2 + 
              (projections$y_proj[i] - projections$y_orig[i])^2) * 0.25
  
  # Create points for an arc representing the theta_centered angle
  arc_points <- 20  # Number of points to use for the arc
  
  # Start angle is 0 (along x-axis from the origin point)
  start_angle <- 0
  # End angle is the theta_centered value (in radians)
  end_angle <- projections$theta_centered[i] * pi/180
  
  # Always draw arcs counterclockwise (standard mathematical convention)
  arc_angles <- seq(start_angle, end_angle, length.out = arc_points)
  
  # Calculate the arc coordinates
  arc_x <- projections$x_orig[i] + r * cos(arc_angles)
  arc_y <- projections$y_orig[i] + r * sin(arc_angles)
  arc_z <- rep(0, arc_points)  # The arc is in the x-y plane
  
  # Add text annotation for the theta_centered angle value
  mid_idx <- floor(arc_points/2)
  fig <- fig %>% add_trace(
    x = c(arc_x[mid_idx]), y = c(arc_y[mid_idx]), z = c(arc_z[mid_idx]),
    type = 'scatter3d', mode = 'text',
    text = paste(round(projections$theta_centered[i], 1), "°"),
    textfont = list(color = 'black', size = 12),
    showlegend = FALSE
  )
}

# Add vertical lines connecting original points to their projections on the x-y plane
for (i in 1:nrow(projections)) {
  fig <- fig %>% add_trace(
    x = c(projections$x_orig[i], projections$x_orig[i]),
    y = c(projections$y_orig[i], projections$y_orig[i]),
    z = c(projections$z_orig[i], 0),
    type = 'scatter3d', mode = 'lines',
    line = list(color = 'gray', width = 1, dash = 'dash'),
    showlegend = (i == 1),  # Only show one legend entry
    name = 'Vertical Projection'
  )
}

# Add lines from projections on the x-y plane to the projected points on the ellipse
for (i in 1:nrow(projections)) {
  fig <- fig %>% add_trace(
    x = c(projections$x_orig[i], projections$x_proj[i]),
    y = c(projections$y_orig[i], projections$y_proj[i]),
    z = c(0, 0),
    type = 'scatter3d', mode = 'lines',
    line = list(color = 'darkgreen', width = 2),
    showlegend = (i == 1),  # Only show one legend entry
    name = 'Projection to Ellipse'
  )
}

# Add tangent vectors at the projected points
# Scale factor for better visualization
scale_factor <- 0.5

for (i in 1:nrow(projections)) {
  # Tangent vector (derivative)
  fig <- fig %>% add_trace(
    x = c(projections$x_proj[i], projections$x_proj[i] + scale_factor * projections$dx_dt_norm[i]),
    y = c(projections$y_proj[i], projections$y_proj[i] + scale_factor * projections$dy_dt_norm[i]),
    z = c(0, 0),
    type = 'scatter3d', mode = 'lines',
    line = list(color = 'red', width = 4),
    showlegend = (i == 1),  # Only show one legend entry
    name = 'Tangent Vector (Derivative)'
  )
  
  # Normal vector
  fig <- fig %>% add_trace(
    x = c(projections$x_proj[i], projections$x_proj[i] + scale_factor * projections$nx[i]),
    y = c(projections$y_proj[i], projections$y_proj[i] + scale_factor * projections$ny[i]),
    z = c(0, 0),
    type = 'scatter3d', mode = 'lines',
    line = list(color = 'blue', width = 4),
    showlegend = (i == 1),  # Only show one legend entry
    name = 'Normal Vector'
  )
}

# Set up the layout
fig <- fig %>% layout(
  scene = list(
    xaxis = list(title = "X"),
    yaxis = list(title = "Y"),
    zaxis = list(title = "Z"),
    aspectmode = 'data'
  ),
  title = "3D Visualization of Points Projected onto an Ellipse"
)

# Print the visualization
fig

# Projection Update with Priors ================================================
# Define the prior on theta 
theta_prior <- 90

n <- nrow(projections)

calc_posterior_theta <- function(theta_centered, derivative_magnitude, theta_prior){
  for (i in 1:n){
    likelihood <- Each point (x,y,z) parameterized as (r, theta, phi)
    prior <- rnorm(1, mean = theta_centered, sd = sqrt(1/derivative_magnitude))
    posterior <- 
  }
}
