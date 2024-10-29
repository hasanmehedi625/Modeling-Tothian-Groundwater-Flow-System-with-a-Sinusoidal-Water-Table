rm(list=ls())

library("fields")

# Grid direction ------------------------------------------------
Lx = 20000 # length of the problem domain along x axis
Lz = 10000 # height of the problem domain along z axis

#nx=80 # number of grids in x direction
#nz = 40 # number of grids in z direction

dx = 250 #Lx/nx #length of grids in x axis
dz = 250 #Lz/nz # ?grid spacing in z direction

x=seq(0, Lx, dx)

z=seq(0, Lz, dz)

nx =length(x)
nz =length(z)


# Parameters --------------------------------------------------------------

a = 50 # amplitude of the wave along meter
l = 5000 # It is (lambda) the period of the wave
s = 0.02 # It is slope of the water table, it is actually tan delta
b = 2*pi/l #frequency


# Calculation of the water table ------------------------------------------

ht = Lz + x*s + a*{sin(b*x/cos(s))}/cos(s)

plot(x, ht, type = "l", col = "red",
     xlab = "Distance", ylab = "head over time",
     main = "Head over Distance")



# Setup of calculation grid and boundary condition ------------------------

H = matrix(0, nrow = nz+1, ncol = nx+2)


H[1, 2:(nx+1)] = ht


tolerance = 0.001 #Level of differences between two cells
max_iteration = 50000 #Maximum number of iterations


# Gauss-Seidel Iteration to solve Laplace's equation ----------------------

for (iter in 1:max_iteration){
  H[(nz+1), ] = H[nz,] #no-flow boundary at bottom
  H[, 1] = H[,2] #no-flow boundary at the left side
  H[, (nx+2)] = H[, (nx+1)] # no-flow right side
  H_old = H
  
  for (i in 2: nz){
    for(j in 2: (nx+1)){
      H[i,j]= (((H[i-1,j]+ H[i+1,j])*dz^2)+ ((H[i, j-1]+H[i, j+1])*dx^2))/(2*(dx^2+dz^2))
    }
  }
  
  
  if (max(abs(H-H_old))<tolerance){
    cat("Converged in", iter, "iterations.\n")
    break
  }
}

Hsim = H[1:nz, 2:(nx+1)]


Hplot=t(Hsim)
Hplot = Hplot[,ncol(Hplot):1]

my_colors <- viridis::viridis(100, alpha = 0.2)  # Set alpha for transparency

image.plot(x, z, Hplot, col = my_colors,xlab = "Distance (m)", 
           ylab = "Depth (m)", 
           main = "Hydraulic Head Distribution (Gauss-Seidel)")
# Add contour lines for hydraulic head with defined intervals
contour_levels <- seq(min(Hplot, na.rm = TRUE), max(Hplot, na.rm = TRUE), 
                      by = 10)
contour(x, z, Hplot, add = TRUE, col = "black", levels = contour_levels)
lines(x,ht,col="red")
