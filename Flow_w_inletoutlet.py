import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
import os

# Physical parameters
rho = 1.225      # kg/m^3, air density
nu = 1.5e-5      # m^2/s, kinematic viscosity of air

# Domain parameters
Lx = 1.0         # Length of domain in x-direction
Ly = 1.0         # Length of domain in y-direction
dx = 0.02
dy = 0.02
nx = int(Lx / dx) + 1
ny = int(Ly / dy) + 1

x = np.linspace(0, Lx, nx)
y = np.linspace(0, Ly, ny)
X, Y = np.meshgrid(x, y)

# Time parameters
dt = 0.001       # Time step
nt = 500         # Number of time steps

# Initialize variables
u = np.zeros((ny, nx))
v = np.zeros((ny, nx))
p = np.zeros((ny, nx))

# Obstacle parameters
obstacle_width = 0.2   # Width of the obstacle
obstacle_height = 0.2  # Height of the obstacle
obstacle_x = Lx / 2 - obstacle_width / 2  # Centered in x
obstacle_y = 0

# Create masks for obstacle, inlet, and outlet
obstacle = np.zeros_like(X, dtype=bool)
obstacle_inlet = np.zeros_like(X, dtype=bool)
obstacle_outlet = np.zeros_like(X, dtype=bool)

for i in range(ny):
    for j in range(nx):
        # Entire obstacle
        if (obstacle_x <= X[i, j] <= obstacle_x + obstacle_width and
            obstacle_y <= Y[i, j] <= obstacle_y + obstacle_height):
            obstacle[i, j] = True

        # Inlet (left edge of the obstacle)
        if (obstacle_x <= X[i, j] <= obstacle_x + 0.1 * obstacle_width and
            obstacle_y <= Y[i, j] <= obstacle_y + obstacle_height):
            obstacle_inlet[i, j] = True

        # Outlet (right edge of the obstacle)
        if (obstacle_x + 0.9 * obstacle_width <= X[i, j] <= obstacle_x + obstacle_width and
            obstacle_y <= Y[i, j] <= obstacle_y + obstacle_height):
            obstacle_outlet[i, j] = True


# Boundary conditions function
def apply_boundary_conditions(u, v, p):
    u[:, 0] = 0.1         # Inflow at left boundary
    u[:, -1] = u[:, -2]   # Outflow at right boundary
    u[0, :] = -u[1, :]    # No-slip at bottom wall
    u[-1, :] = u[-2, :]   # Free shear at top boundary

    v[:, 0] = 0           # No vertical velocity at inflow
    v[:, -1] = v[:, -2]   # Outflow at right boundary
    v[0, :] = 0           # No-slip at bottom wall
    v[-1, :] = v[-2, :]   # Free shear at top boundary

    # Obstacle boundary conditions
    u[obstacle] = 0
    v[obstacle] = 0

    # Inlet: Prescribe velocity
    u[obstacle_inlet] = 0.2  # Specify inlet velocity
    v[obstacle_inlet] = 0

    # Outlet: Zero-gradient conditions
    u[obstacle_outlet] = u[obstacle_outlet]  # Maintain outlet velocity
    v[obstacle_outlet] = v[obstacle_outlet]  # Maintain outlet velocity

    # Pressure boundary conditions
    p[:, 0] = p[:, 1]     # Zero pressure gradient at inflow
    p[:, -1] = p[:, -2]   # Zero pressure gradient at outflow
    p[0, :] = p[1, :]     # Zero pressure gradient at bottom wall
    p[-1, :] = p[-2, :]   # Zero pressure gradient at top boundary

    # Enforce pressure at obstacle
    p[obstacle] = p[obstacle]

    return u, v, p


# Create a directory to save plots
save_dir = 'velocity_field_plots'
if not os.path.exists(save_dir):
    os.makedirs(save_dir)

# Time-stepping loop
for n in range(nt):
    un = u.copy()
    vn = v.copy()
    pn = p.copy()

    # Apply boundary conditions
    u, v, p = apply_boundary_conditions(u, v, p)

    # Calculate the RHS of the pressure Poisson equation
    b = np.zeros_like(p)
    b[1:-1,1:-1] = (rho * ((1/dt) *
                   ((un[1:-1,2:] - un[1:-1,0:-2])/(2*dx) +
                    (vn[2:,1:-1] - vn[0:-2,1:-1])/(2*dy)) -
                   ((un[1:-1,2:] - un[1:-1,0:-2])/(2*dx))**2 -
                   2 * ((un[2:,1:-1] - un[0:-2,1:-1])/(2*dy) *
                        (vn[1:-1,2:] - vn[1:-1,0:-2])/(2*dx)) -
                   ((vn[2:,1:-1] - vn[0:-2,1:-1])/(2*dy))**2))

    # Pressure Poisson equation solver
    for it in range(50):  # Number of Poisson iterations per time step
        pn = p.copy()
        p[1:-1,1:-1] = (((pn[1:-1,2:] + pn[1:-1,0:-2]) * dy**2 +
                         (pn[2:,1:-1] + pn[0:-2,1:-1]) * dx**2) /
                        (2 * (dx**2 + dy**2)) -
                        dx**2 * dy**2 / (2 * (dx**2 + dy**2)) * b[1:-1,1:-1])

        # Apply pressure boundary conditions
        p = apply_boundary_conditions(u, v, p)[2]

    # Velocity field updates
    u[1:-1,1:-1] = (un[1:-1,1:-1] -
                    dt / dx * un[1:-1,1:-1] * (un[1:-1,1:-1] - un[1:-1,0:-2]) -
                    dt / dy * vn[1:-1,1:-1] * (un[1:-1,1:-1] - un[0:-2,1:-1]) -
                    dt / (2 * rho * dx) * (p[1:-1,2:] - p[1:-1,0:-2]) +
                    nu * (dt / dx**2 * (un[1:-1,2:] - 2 * un[1:-1,1:-1] + un[1:-1,0:-2]) +
                          dt / dy**2 * (un[2:,1:-1] - 2 * un[1:-1,1:-1] + un[0:-2,1:-1])))

    v[1:-1,1:-1] = (vn[1:-1,1:-1] -
                    dt / dx * un[1:-1,1:-1] * (vn[1:-1,1:-1] - vn[1:-1,0:-2]) -
                    dt / dy * vn[1:-1,1:-1] * (vn[1:-1,1:-1] - vn[0:-2,1:-1]) -
                    dt / (2 * rho * dy) * (p[2:,1:-1] - p[0:-2,1:-1]) +
                    nu * (dt / dx**2 * (vn[1:-1,2:] - 2 * vn[1:-1,1:-1] + vn[1:-1,0:-2]) +
                          dt / dy**2 * (vn[2:,1:-1] - 2 * vn[1:-1,1:-1] + vn[0:-2,1:-1])))

    # Apply boundary conditions
    u, v, p = apply_boundary_conditions(u, v, p)

    # Enforce obstacle conditions (no-slip)
    u[obstacle] = 0
    v[obstacle] = 0

# Save plots at every 50 time steps
#took this out of the loop for faster testing 
if n % 50 == 0:
    plt.figure(figsize=(11, 7))
    plt.contourf(X, Y, np.sqrt(u**2 + v**2), alpha=0.5, cmap="viridis")  # Velocity magnitude
    plt.colorbar(label='Velocity magnitude')
    plt.quiver(X[::3, ::3], Y[::3, ::3], u[::3, ::3], v[::3, ::3], scale=10)  # Velocity vectors
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title('Velocity field at time step {}'.format(n))
    
    # Plot full obstacle
    plt.plot([obstacle_x, obstacle_x + obstacle_width, obstacle_x + obstacle_width, obstacle_x, obstacle_x],
             [obstacle_y, obstacle_y, obstacle_y + obstacle_height, obstacle_y + obstacle_height, obstacle_y],
             'k-', linewidth=2, label="Obstacle")  # Black outline for obstacle

    # Highlight inlet region 
    inlet_start_y = obstacle_y + 0.45 * obstacle_height  # Start point (middle of the boundary)
    inlet_end_y = obstacle_y + 0.55 * obstacle_height    # End point (10% span)
    plt.plot([obstacle_x, obstacle_x],  # x-coordinates
             [inlet_start_y, inlet_end_y],  # y-coordinates
             'b-', linewidth=2, label="Inlet")  # Inlet as blue line

    # Highlight outlet region 
    outlet_start_y = obstacle_y + 0.45 * obstacle_height  # Start point (middle of the boundary)
    outlet_end_y = obstacle_y + 0.55 * obstacle_height    # End point (10% span)
    plt.plot([obstacle_x + obstacle_width, obstacle_x + obstacle_width],  # x-coordinates
             [outlet_start_y, outlet_end_y],  # y-coordinates
             'r-', linewidth=2, label="Outlet")  # Outlet as red line

    plt.legend()
    filename = os.path.join(save_dir, 'velocity_field_step_{}.png'.format(n))  # Save plot
    plt.savefig(filename)
    plt.close()  # Close the figure to free up memory

# Final visualization
plt.figure(figsize=(11, 7))
plt.contourf(X, Y, np.sqrt(u**2 + v**2), alpha=0.5, cmap="viridis")  # Velocity magnitude
plt.colorbar(label='Velocity magnitude')
plt.quiver(X[::3, ::3], Y[::3, ::3], u[::3, ::3], v[::3, ::3], scale=10)  # Velocity vectors
plt.xlabel('X')
plt.ylabel('Y')
plt.title('Velocity field after {} time steps'.format(nt))

# Plot full obstacle as black outline
plt.plot([obstacle_x, obstacle_x + obstacle_width, obstacle_x + obstacle_width, obstacle_x, obstacle_x],
         [obstacle_y, obstacle_y, obstacle_y + obstacle_height, obstacle_y + obstacle_height, obstacle_y],
         'k-', linewidth=2, label="Obstacle")  # Black outline for obstacle

# Highlight inlet region 
inlet_start_y = obstacle_y + 0.45 * obstacle_height  # Start point 
inlet_end_y = obstacle_y + 0.55 * obstacle_height    # End point 
plt.plot([obstacle_x, obstacle_x],  # x-coordinates
         [inlet_start_y, inlet_end_y],  # y-coordinates
         'b-', linewidth=2, label="Inlet")  

# Highlight outlet region 
outlet_start_y = obstacle_y + 0.45 * obstacle_height  # Start point 
outlet_end_y = obstacle_y + 0.55 * obstacle_height    # End point 
plt.plot([obstacle_x + obstacle_width, obstacle_x + obstacle_width],  # x-coordinates
         [outlet_start_y, outlet_end_y],  # y-coordinates
         'r-', linewidth=2, label="Outlet") 

plt.legend()
plt.savefig(os.path.join(save_dir, 'velocity_field_final.png'))  # Save final plot
plt.show()
