# ==============================================================================
# Title: Circular Package
# Author: Daniel Posmik
# E-Mail: daniel_posmik@brown.edu
# ==============================================================================

# Set-Up =======================================================================
library(dplyr)
library(tidyverse)
library(circular)

# Simulation ===================================================================
# Set random seed for reproducibility
set.seed(1)

# Parameters
n <- 500
mu <- circular(pi, units = "radians")
rho_values <- c(0.3, 0.6, 0.9)
colors <- c("blue", "red", "green")

# Create an empty data frame to store results
results <- data.frame(
  angle = numeric(0),
  rho = numeric(0)
)

# Generate random wrapped normal data with different concentrations
for (i in 1:length(rho_values)) {
  samples <- circular::rwrappednormal(n, mu = mu, rho = rho_values[i])
  results <- rbind(results, data.frame(
    angle = as.numeric(samples),
    rho = rep(as.character(rho_values[i]), n)
  ))
}

# Circular histogram
par(mfrow = c(2, 2))

# 1. Rose diagram for all concentrations
plot(circular(results$angle),
  stack = TRUE, bins = 36,
  main = "Rose Diagram of Wrapped Normal Distributions",
  col = colors[as.factor(results$rho)]
)
legend("topright",
  legend = paste("rho =", rho_values),
  fill = colors, cex = 0.8, bg = "white"
)

# 2. Individual rose diagrams for each concentration
for (i in 1:length(rho_values)) {
  plot(circular(results$angle[results$rho == as.character(rho_values[i])]),
    stack = TRUE, bins = 36,
    main = paste("Rose Diagram (rho =", rho_values[i], ")"),
    col = colors[i]
  )
}

# Linear plot with density curves
par(mfrow = c(1, 1))

# Convert to data frame for ggplot
df <- data.frame(
  value = results$angle,
  concentration = factor(results$rho, levels = as.character(rho_values))
)

# Create density plot
ggplot(df, aes(x = value, fill = concentration)) +
  geom_density(alpha = 0.5) +
  scale_x_continuous(
    limits = c(0, 2 * pi),
    breaks = c(0, pi / 2, pi, 3 * pi / 2, 2 * pi),
    labels = c("0", "π/2", "π", "3π/2", "2π")
  ) +
  labs(
    title = "Density of Wrapped Normal Distributions",
    x = "Angle (radians)",
    y = "Density",
    fill = "Concentration (ρ)"
  ) +
  theme_minimal() +
  scale_fill_manual(values = colors)

# 3D Visualization on the unit circle
library(plotly)

# Create points on the unit circle
theta <- seq(0, 2 * pi, length.out = 100)
circle_x <- cos(theta)
circle_y <- sin(theta)

# Create points for the distributions
dist_points <- list()
for (i in 1:length(rho_values)) {
  angles <- results$angle[results$rho == as.character(rho_values[i])]
  dist_points[[i]] <- data.frame(
    x = cos(angles),
    y = sin(angles),
    rho = rep(rho_values[i], length(angles))
  )
}

# Combine all points
all_points <- do.call(rbind, dist_points)

# Create the plot
p <- plot_ly() %>%
  # Add unit circle
  add_trace(
    x = circle_x, y = circle_y,
    type = "scatter", mode = "lines",
    line = list(color = "black", width = 1),
    name = "Unit Circle"
  )

# Add points for each concentration
for (i in 1:length(rho_values)) {
  p <- p %>% add_trace(
    x = dist_points[[i]]$x, y = dist_points[[i]]$y,
    type = "scatter", mode = "markers",
    marker = list(color = colors[i], size = 3),
    name = paste("ρ =", rho_values[i])
  )
}

# Add mean direction
p <- p %>% add_trace(
  x = c(0, cos(as.numeric(mu))), y = c(0, sin(as.numeric(mu))),
  type = "scatter", mode = "lines",
  line = list(color = "black", width = 2, dash = "dash"),
  name = "Mean Direction"
)

# Configure layout
p <- p %>% layout(
  title = "Wrapped Normal Distributions on Unit Circle",
  xaxis = list(title = "", zeroline = TRUE, scaleanchor = "y"),
  yaxis = list(title = "", zeroline = TRUE),
  showlegend = TRUE
)

print(p)
