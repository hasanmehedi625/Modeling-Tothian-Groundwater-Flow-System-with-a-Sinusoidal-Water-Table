# Code setup --------------------------------------------------------------
rm(list = ls())

library('fields')
library('viridis')


# Declaration of Variables ------------------------------------------------

w = 40 # pumping rate (discharge per unit area) of the well.
t = 100 # Ttransmissivity of the aquifer (m²/day).
h = 50

Lx = 1000 # Length of the domain in the x direction (m)
Lz = 1000 # Depth of the domain in the z direction (m)
dx = 10 # Grid spacing in the x direction
dz = 10 # Grid spacing in the z direction

nx = (Lx/dx) # Number of grid points in the x direction
nz = (Lz/dz) # Number of grid points in the z direction (depth)
iw = 50
jw = 50
# Create the x and z grids
x = seq(0, Lx, length.out = nx)  #length.out = nx
z = seq(0, Lz, length.out = nz) 



# Calculation -------------------------------------------------------------

# Initialize the hydraulic head matrix (H)
H = matrix(0, nrow = nz + 2, ncol = nx + 2)

# set the well condiotion
# H[(iw + 1), (jw + 1)] = 80

# Set boundary conditions: Sinusoidal head at the surface (water table), no-flow at the sides and bottom
H[2:(nz+1),1] = 50
H[2:(nz+1),(nx + 2)] = 50

# Define convergence tolerance and maximum number of iterations
tolerance = 1e-3  # Convergence criterion
max_iter = 50000  # Maximum number of iterations


# Gauss-Seidel Iteration to solve Laplace's equation

for(iter in 1: max_iter){
  H[1, ] = H[2, ] # No-flow boundary at the bottom
  H[(nz+2), ] = H[(nz+1), ] # No-flow boundary at the bottom
  H_old = H
  
  # Update the head values using Gauss-Seidel iteration
  for(i in 1:nz+1) {
    for(j in 2:(nx + 1)){
      if (i == (iw + 1) && j == (jw + 1)){
        # set the well condiotion
        H[i, j] = (((H[i-1, j] + H[i+1, j]) * dz^2 + (H[i, j-1] + H[i, j+1]) * dx^2) / (2 * (dx^2 + dz^2))) - (w/t)
      }
      else {
        H[i, j] = ((H[i-1, j] + H[i+1, j]) * dz^2 + (H[i, j-1] + H[i, j+1]) * dx^2) / (2 * (dx^2 + dz^2))
      }
      
    }
  }
  # Check for convergence
  if (max(abs(H - H_old)) < tolerance) {
    cat("Converged in", iter, "iterations.\n")
    break
  }
}  



Hsim = H[2:(nz + 1),2:(nx+1)]



# Plotting ----------------------------------------------------------------
# Visualize the hydraulic head distribution

Hplot =t(Hsim)
Hplot = Hplot[,ncol(Hplot):1]

my_colors <- viridis::viridis(100, alpha = 0.2)  # Set alpha for transparency

image.plot(x, z, Hplot,
           # ylim = c(0, 10500),
           col = my_colors,
           xlab = "Distance (m)", 
           ylab = "Depth (m)", 
           main = "Hydraulic Head Distribution (Gauss-Seidel)")

# Add contour lines for hydraulic head with defined intervals
contour_levels <- seq(min(Hplot, na.rm = TRUE), max(Hplot, na.rm = TRUE), 
                      by = .5)

contour(x, z, Hplot, add = TRUE, col = "black", levels = contour_levels)

