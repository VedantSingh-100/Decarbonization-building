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
obstacle_y = Ly / 2 - obstacle_height / 2  # Centered in y

# Create obstacle mask
obstacle = np.zeros_like(X, dtype=bool)
for i in range(ny):
    for j in range(nx):
        if (obstacle_x <= X[i, j] <= obstacle_x + obstacle_width and
            obstacle_y <= Y[i, j] <= obstacle_y + obstacle_height):
            obstacle[i, j] = True

# Boundary conditions function
def apply_boundary_conditions(u, v, p):
    # Apply boundary conditions for u
    u[:, 0] = 0.1         # Inflow at left boundary
    u[:, -1] = u[:, -2]   # Outflow at right boundary
    u[0, :] = -u[1, :]    # No-slip at bottom wall
    u[-1, :] = u[-2, :]   # Free shear at top boundary

    # Apply boundary conditions for v
    v[:, 0] = 0           # No vertical velocity at inflow
    v[:, -1] = v[:, -2]   # Outflow at right boundary
    v[0, :] = 0           # No-slip at bottom wall
    v[-1, :] = v[-2, :]   # Free shear at top boundary

    # Apply boundary conditions at the obstacle (no-slip walls)
    u[obstacle] = 0
    v[obstacle] = 0

    # Pressure boundary conditions
    p[:, 0] = p[:, 1]     # Zero pressure gradient at inflow
    p[:, -1] = p[:, -2]   # Zero pressure gradient at outflow
    p[0, :] = p[1, :]     # Zero pressure gradient at bottom wall
    p[-1, :] = p[-2, :]   # Zero pressure gradient at top boundary

    # Enforce pressure at obstacle (could be Neumann or Dirichlet)
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
    if n % 50 == 0:
        plt.figure(figsize=(11,7))
        plt.contourf(X, Y, np.sqrt(u**2 + v**2), alpha=0.5)
        plt.colorbar(label='Velocity magnitude')
        plt.quiver(X[::3, ::3], Y[::3, ::3], u[::3, ::3], v[::3, ::3])
        plt.xlabel('X')
        plt.ylabel('Y')
        plt.title('Velocity field at time step {}'.format(n))
        # Plot obstacle
        plt.plot([obstacle_x, obstacle_x + obstacle_width, obstacle_x + obstacle_width, obstacle_x, obstacle_x],
                 [obstacle_y, obstacle_y, obstacle_y + obstacle_height, obstacle_y + obstacle_height, obstacle_y],
                 'k-', linewidth=2)
        # Save the figure
        filename = os.path.join(save_dir, 'velocity_field_step_{}.png'.format(n))
        plt.savefig(filename)
        plt.close()  # Close the figure to free up memory

# Final visualization
plt.figure(figsize=(11,7))
plt.contourf(X, Y, np.sqrt(u**2 + v**2), alpha=0.5)
plt.colorbar(label='Velocity magnitude')
plt.quiver(X[::3, ::3], Y[::3, ::3], u[::3, ::3], v[::3, ::3])
plt.xlabel('X')
plt.ylabel('Y')
plt.title('Velocity field after {} time steps'.format(nt))
# Plot obstacle
plt.plot([obstacle_x, obstacle_x + obstacle_width, obstacle_x + obstacle_width, obstacle_x, obstacle_x],
         [obstacle_y, obstacle_y, obstacle_y + obstacle_height, obstacle_y + obstacle_height, obstacle_y],
         'k-', linewidth=2)
plt.savefig(os.path.join(save_dir, 'velocity_field_final.png'))
plt.show()
