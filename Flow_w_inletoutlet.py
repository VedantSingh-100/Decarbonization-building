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

def compute_pressure_gradients(p, dx, dy):
    """Compute the pressure gradients."""
    dpdx = np.gradient(p, axis=1) / dx
    dpdy = np.gradient(p, axis=0) / dy
    return dpdx, dpdy


def plot_velocity_pressure_vectors(X, Y, u, v, p, dx, dy, obstacle_coords, n=None, save_path=None):
    """Plot velocity and pressure vectors."""
    plt.figure(figsize=(11, 7))
    velocity_mag = np.sqrt(u**2 + v**2)
    dpdx, dpdy = compute_pressure_gradients(p, dx, dy)
    pressure_mag = np.sqrt(dpdx**2 + dpdy**2)

    # Contour plot for velocity magnitude
    plt.contourf(X, Y, velocity_mag, levels=50, alpha=0.8, cmap="viridis")
    plt.colorbar(label='Velocity Magnitude')
    
    # Velocity vectors
    skip = 5
    plt.quiver(X[::skip, ::skip], Y[::skip, ::skip], u[::skip, ::skip], v[::skip, ::skip],
               color='white', scale=50, pivot='middle', alpha=0.8, label="Velocity Vectors")
    
    # Pressure vectors
    plt.quiver(X[::skip, ::skip], Y[::skip, ::skip], dpdx[::skip, ::skip], dpdy[::skip, ::skip],
               color='red', scale=100, pivot='middle', alpha=0.8, label="Pressure Vectors")

    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title(f'Velocity and Pressure Vectors {"at time step {}".format(n) if n else ""}')
    
    # Draw obstacle
    plt.plot(*obstacle_coords, 'k-', linewidth=2, label="Obstacle")
    plt.legend(loc='upper right')
    
    if save_path:
        plt.savefig(save_path)
        plt.close()
    else:
        plt.show()


def update_particle_positions(particle_positions, u, v, dt):
    """Update particle positions based on velocity field."""
    x, y = particle_positions[:, 0], particle_positions[:, 1]
    x_idx = np.clip((x / dx).astype(int), 0, u.shape[1] - 1)
    y_idx = np.clip((y / dy).astype(int), 0, u.shape[0] - 1)
    
    # Get velocity at particle positions
    u_particles = u[y_idx, x_idx]
    v_particles = v[y_idx, x_idx]
    
    # Update positions
    x_new = x + u_particles * dt
    y_new = y + v_particles * dt
    
    # Reflect particles that hit the boundaries
    x_new = np.clip(x_new, 0, Lx)
    y_new = np.clip(y_new, 0, Ly)
    return np.vstack([x_new, y_new]).T


def plot_particles(X, Y, u, v, particle_positions, obstacle_coords, n=None, save_path=None):
    """Plot particles over the velocity field."""
    plt.figure(figsize=(11, 7))
    velocity_mag = np.sqrt(u**2 + v**2)
    
    # Contour plot for velocity magnitude
    plt.contourf(X, Y, velocity_mag, levels=50, alpha=0.5, cmap="viridis")
    plt.colorbar(label='Velocity Magnitude')
    
    # Plot particles
    plt.scatter(particle_positions[:, 0], particle_positions[:, 1], c='red', s=10, label='Particles')

    # Velocity vectors for reference
    skip = 5
    plt.quiver(X[::skip, ::skip], Y[::skip, ::skip], u[::skip, ::skip], v[::skip, ::skip],
               color='white', scale=50, pivot='middle', alpha=0.6, label="Velocity Vectors")
    
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title(f'Particles in Flow Field {"at time step {}".format(n) if n else ""}')
    
    # Draw obstacle
    plt.plot(*obstacle_coords, 'k-', linewidth=2, label="Obstacle")
    plt.legend(loc='upper right')
    
    if save_path:
        plt.savefig(save_path)
        plt.close()
    else:
        plt.show()


# Define obstacle coordinates for plotting
obstacle_coords = ([obstacle_x, obstacle_x + obstacle_width, obstacle_x + obstacle_width, obstacle_x, obstacle_x],
                   [obstacle_y, obstacle_y, obstacle_y + obstacle_height, obstacle_y + obstacle_height, obstacle_y])

# Generate initial particle positions (e.g., grid of particles)
particle_positions = np.random.rand(100, 2) * [Lx, Ly]  # 100 particles randomly placed

# Example usage
plot_velocity_pressure_vectors(X, Y, u, v, p, dx, dy, obstacle_coords, n=nt,
                                save_path=os.path.join(save_dir, 'velocity_pressure_vectors.png'))

# Update and plot particles
for t in range(10):  # Example: update and plot for 10 steps
    particle_positions = update_particle_positions(particle_positions, u, v, dt)
    plot_particles(X, Y, u, v, particle_positions, obstacle_coords, n=t,
                   save_path=os.path.join(save_dir, f'particles_step_{t}.png'))
