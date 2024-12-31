rm=(list=ls())
# Parameters
Lx <- 1000     # Length of the domain (m)
Ly <- 1000     # Width of the domain (m)
T <- 100       # Transmissivity (m²/day)
W <- 0.05      # Pumping rate (m³/day per unit area)
xwell <- 500   # Well location (m)
ywell <- 500   # Well location (m)
h_fixed <- 50  # Fixed head value at the boundary (m)

# Grid size
Nx <- 50       # Number of grid points in x direction
Ny <- 50       # Number of grid points in y direction
dx <- Lx / Nx  # Grid spacing in x direction
dy <- Ly / Ny  # Grid spacing in y direction

# Initialize the hydraulic head grid
h <- matrix(NA, nrow = Ny, ncol = Nx)

# Apply boundary conditions
h[1, ] <- h_fixed         # Top boundary (No flow)
h[Ny, ] <- h_fixed       # Bottom boundary (No flow)
h[, 1] <- h_fixed         # Left boundary (Fixed head)
h[, Nx] <- h_fixed        # Right boundary (Fixed head)

# Initial guess for the hydraulic head inside the domain
h[2:(Ny-1), 2:(Nx-1)] <- h_fixed

# Gauss-Seidel Method to solve for hydraulic head
max_iter <- 1000          # Maximum number of iterations
tolerance <- 1e-4         # Convergence tolerance

for (iter in 1:max_iter) {
  h_old <- h  # Save the old hydraulic head
  
  # Update the hydraulic head for each grid point using the Gauss-Seidel method
  for (i in 2:(Ny-1)) {
    for (j in 2:(Nx-1)) {
      # Applying the 2D groundwater flow equation (finite difference)
      # The well term W is applied to the grid cell at the well location
      if (i == ywell/dy && j == xwell/dx) {
        h[i, j] <- (h[i+1, j] + h[i-1, j] + h[i, j+1] + h[i, j-1] - W * dx * dy / T) / 4
      } else {
        h[i, j] <- (h[i+1, j] + h[i-1, j] + h[i, j+1] + h[i, j-1]) / 4
      }
    }
  }
  
  # Check for convergence (absolute difference between old and new head)
  if (max(abs(h - h_old)) < tolerance) {
    cat("Converged in", iter, "iterations\n")
    break
  }
}



# Contour plot

contour(x = seq(0, Lx, length.out = Nx), y = seq(0, Ly, length.out = Ny), 
        z = h, xlab = "X (m)", ylab = "Y (m)", main = "Hydraulic Head Distribution",
        col = terrain.colors(100), lwd = 2)


# Heatmap plot
image(x = seq(0, Lx, length.out = Nx), y = seq(0, Ly, length.out = Ny), 
      z = h, xlab = "X (m)", ylab = "Y (m)", main = "Hydraulic Head Heatmap",
      col = terrain.colors(100))

# Analyzing the effect of pumping (Increase W to 0.1 m³/day/unit area)

W_new <- 0.1  # New pumping rate

# Reinitialize and solve again for the new pumping rate
h_new <- matrix(NA, nrow = Ny, ncol = Nx)
h_new[1, ] <- h_fixed
h_new[Ny, ] <- h_fixed
h_new[, 1] <- h_fixed
h_new[, Nx] <- h_fixed
h_new[2:(Ny-1), 2:(Nx-1)] <- h_fixed

for (iter in 1:max_iter) {
  h_old <- h_new  # Save the old hydraulic head
  
  # Update the hydraulic head for each grid point
  for (i in 2:(Ny-1)) {
    for (j in 2:(Nx-1)) {
      # Apply new pumping rate term
      if (i == ywell/dy && j == xwell/dx) {
        h_new[i, j] <- (h_new[i+1, j] + h_new[i-1, j] + h_new[i, j+1] + h_new[i, j-1] - W_new * dx * dy / T) / 4
      } else {
        h_new[i, j] <- (h_new[i+1, j] + h_new[i-1, j] + h_new[i, j+1] + h_new[i, j-1]) / 4
      }
    }
  }
  
  # Check for convergence
  if (max(abs(h_new - h_old)) < tolerance) {
    cat("Converged in", iter, "iterations for W =", W_new, "\n")
    break
  }
}

# Plot the new hydraulic head distribution
contour(x = seq(0, Lx, length.out = Nx), y = seq(0, Ly, length.out = Ny), 
        z = h_new, xlab = "X (m)", ylab = "Y (m)", main = paste("Hydraulic Head Distribution (W =", W_new, ")"),
        col = terrain.colors(100), lwd = 2)

image(x = seq(0, Lx, length.out = Nx), y = seq(0, Ly, length.out = Ny), 
      z = h_new, xlab = "X (m)", ylab = "Y (m)", main = paste("Hydraulic Head Heatmap (W =", W_new, ")"),
      col = terrain.colors(100))

