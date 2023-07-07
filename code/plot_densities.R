plot_densities <- function(mu0, sigma0, mu1, sigma1, x_low, x_high) {
  # Create a sequence of x values
  x <- seq(from = x_low, to = x_high, by = 0.01)

  # Compute density values for the two normal distributions
  density0 <- dnorm(x, mean = mu0, sd = sigma0)
  density1 <- dnorm(x, mean = mu1, sd = sigma1)

  # Plot the densities
  plot(x, density0, type = "l", col = "blue",
       xlab = "x", ylab = "Density",
       main = "Normal Densities", ylim = c(0, max(density0, density1)))
  lines(x, density1, type = "l", col = "red")

  # Add a legend
  legend("topright", legend = c("Density 0", "Density 1"),
         lty = c(1, 1), col = c("blue", "red"))
}

# Test the function
plot_densities(mu0 = 0, sigma0 = 1, mu1 = 2, sigma1 = 1.5, x_low = -5, x_high = 5)
