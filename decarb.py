import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.patches import Rectangle


""" vvv Number of iterations and wind profile vvv  and other parameters"""

nt = 10000 #timesteps
apply_wind_profile = True #apply power-law profile (true) or uniform profile)
apply_obstacle_mask = True  # Set to False to disable the obstacle mask
initialize_flow = True

u_r = 8.49      # m/s (reference wind speed)
z_r = 10        # feet (reference height of measurement)

u_vel = u_r    # m/s (uniform wind profile speed)
v_vel = 0     # m/s (uniform wind profile speed)

Lx = 500   # Length of domain in x-direction (VARY THIS IF FLOW TOUCHES BOUNDARY)
Ly = 130     # Length of domain in y-direction (INCREASE THIS IF FLOW TOUCHES BOUNDARY)

dx = 1      # TRY TO KEEP CONSTANT at 1, only change if divergence
dy = 1      # TRY TO KEEP CONSTANT at 1, only change if divergence

dt = 0

diffusion_coefficient = 2*10**-5

n_switch = 500

time = 0
""" ^^^ Number of iterations and wind profile^^^ """

""" Obstacle definition """
### Define obstacle region in the middle of the domain ###

obstacle_x_min = Lx/4
obstacle_x_max = obstacle_x_min+10
obstacle_y_min = 0
obstacle_y_max = 20     # KEEP TOTAL AREA = 200 m2

mult = 2/dx # multiplier to subtract the appropriate amount of cells from inlet wall face to define actual inlet

vol_flow = 180  # m3/s  ###KEEP CONSTANT AT 180 FOR EXPERIMENTS
u_inlet = -vol_flow/(obstacle_y_max - mult*dy)
u_outlet = u_inlet

""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
def add_obstacle_outline(ax, obstacle_x_min, obstacle_x_max, obstacle_y_min, obstacle_y_max):
    width = obstacle_x_max - obstacle_x_min
    height = obstacle_y_max - obstacle_y_min
    obstacle = Rectangle((obstacle_x_min, obstacle_y_min),  width, height, linewidth=1, edgecolor='black', facecolor='black', linestyle='-')
    ax.add_patch(obstacle)

def plot_particles(particle_x_positions,particle_y_positions):
    plt.figure(figsize=(8, 6))
    plt.scatter(particle_x_positions, particle_y_positions, color='red', label='Particles')  # Adjust y
    plt.xlim(0, Lx)
    plt.ylim(0, Ly)
    plt.xlabel('x [m]')
    plt.ylabel('y [m]')
    plt.title('Particle Trajectories')
    plt.gca().set_aspect('equal', adjustable='box')
    plt.legend()
    for xc in x:
        plt.axvline(x=xc, color='black', linestyle='--', linewidth=0.5)
    for yc in y:
        plt.axhline(y=yc, color='black', linestyle='--', linewidth=0.5)
    plt.show()
    
def calculate_dt(u,v):
    u_max = np.max(np.abs(u))
    v_max = np.max(np.abs(v))
    vel_dlength_max = max(u_max/dx,v_max/dy)
    courrant_num_adjuster = 0.3
    dt = 1/vel_dlength_max
    
    # Viscous/diffusive dt
    dt_diff = dx**2/(mu/rho)
        
    dt_diff_sc = dx**2 / (4.0 * diffusion_coefficient)
    
    dt_min = min(dt, dt_diff, dt_diff_sc)

    # Determine which dt is chosen and print an explanation
    if dt_min == dt:
        dt_type = "convective dt"
    elif dt_min == dt_diff:
        dt_type = "viscous/diffusive dt"
    elif dt_min == dt_diff_sc:
        dt_type = "scalar diffusion dt"
    else:
        dt_type = "unknown dt type"
    print(f'n: {n}, dt: {dt_min} ({dt_type})')
    return courrant_num_adjuster*min(dt,dt_diff,dt_diff_sc)

def calculate_courrant(u,v,dt):
    u_max = np.max(np.abs(u))
    v_max = np.max(np.abs(v))
    nu_u = u_max*dt/dx
    nu_v = v_max*dt/dy
    nu = np.max([nu_u,nu_v])
    return nu

def apply_boundary_conditions(u, v, p, c):
    # Left boundary, inflow BC
    if apply_wind_profile:
        y_padded = np.linspace(y[0] - dy / 2, y[-1] + dy / 2, len(u[:, 0]))
        u[:, 0] = u_r * (y_padded / z_r) ** alpha
    else:
        u[:,0] = u_vel      # velocity at left-side boundary
        v[1,:] = v_vel      # velocity at left-side boundary
    
    # Left boundary, inflow BC (continued)
    p[:,0] = p[:,1]     # dP/dx|x=0 = 0
    v[:,0] = 0         # v|x=0 = 0
    
    # Right boundary, outflow BC
    u[:,-1] = u[:,-2]   # u(M+1,1) = u(M,1)
    p[:,-1] = 0         # P|x=Lx = 0 (0 gauge outlet)
    # p[:,-1] = p[:,-2]  # P(M+1,1) = -P(M,1)
    v[:,-1] = v[:,-2]   # dv/dx = 0
    
    # Top boundary, pressure outlet
    u[-1, :] = u[-2, :]  # Zero gradient for horizontal velocity
    v[-1, :] = v[-2, :]  # Zero gradient for vertical velocity
    p[-1, :] = 0         # Zero gauge pressure
     
    
    # Bottom boundary, wall
    u[0,:] = -u[1,:]    # no slip at bottom wall
    p[0,:] = p[1,:]     # bottom wall, horizontal wall dP/dy = 0 
    v[0,:] = 0          # bottom wall (ground), v = 0

    return u, v, p


def add_quiver(ax, X, Y, u, v, scale=1000, skip=10, color='grey', width=0.001):
    skip_slices = (slice(None, None, skip), slice(None, None, skip))
    ax.quiver(
        X[skip_slices], Y[skip_slices], u[skip_slices], v[skip_slices],
        color=color, scale=scale, width=width)

    
def plot_velocity(u, v):
    u_centered = 0.5 * (u[:-1, :-1] + u[1:, :-1])   # Average to cell centers in x-direction
    v_centered = 0.5 * (v[:, :-1] + v[:, 1:])       # Average to cell centers in y-direction
    velocity_magnitude = np.sqrt(u_centered**2 + v_centered**2)

    x_plot = x[:u_centered.shape[1]]
    y_plot = y[:u_centered.shape[0]]
    X_plot, Y_plot = np.meshgrid(x_plot, y_plot)

    fig, ax = plt.subplots(figsize=(12, 4)) 
    contour = ax.contourf(X_plot, Y_plot, velocity_magnitude, levels=np.linspace(0, np.max(velocity_magnitude), 21), cmap='viridis')
    fig.colorbar(contour, ax=ax, label='Velocity Magnitude [m/s]')

    add_quiver(ax, X_plot, Y_plot, u_centered, v_centered, width=0.002)

    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_title(f'Velocity magnitude at t = {time:0.02f}[s]')


    ax.set_xlim([0, Lx])
    ax.set_ylim([0, Ly])
    ax.set_aspect('equal') 
    add_obstacle_outline(ax, obstacle_x_min, obstacle_x_max, obstacle_y_min, obstacle_y_max)

    plt.subplots_adjust(left=0.05, right=0.95, top=0.9, bottom=0.1)  
    plt.show()

def plot_u_velocity(u, v):
    u_centered = 0.5 * (u[:-1, :-1] + u[1:, :-1])  
    x_plot = x[:u_centered.shape[1]]
    y_plot = y[:u_centered.shape[0]]
    X_plot, Y_plot = np.meshgrid(x_plot, y_plot)

    fig, ax = plt.subplots(figsize=(12, 4))  
    contour = ax.contourf(X_plot, Y_plot, u_centered, levels=np.linspace(np.min(u_centered), np.max(u_centered), 21), cmap='coolwarm')
    fig.colorbar(contour, ax=ax, label='u-Velocity Magnitude [m/s]')

    # add_quiver(ax, X_plot, Y_plot, u_centered, v_centered, width=0.002)

    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_title(f'u-Velocity at t = {time:0.02f}[s]')


    ax.set_xlim([0, Lx])
    ax.set_ylim([0, Ly])
    ax.set_aspect('equal')  
    add_obstacle_outline(ax, obstacle_x_min, obstacle_x_max, obstacle_y_min, obstacle_y_max)
    
    plt.subplots_adjust(left=0.05, right=0.95, top=0.9, bottom=0.1)  
    plt.show()


def plot_v_velocity(u, v):
    v_centered = 0.5 * (v[:, :-1] + v[:, 1:])    
    x_plot = x[:v_centered.shape[1]]
    y_plot = y[:v_centered.shape[0]]
    X_plot, Y_plot = np.meshgrid(x_plot, y_plot)

    fig, ax = plt.subplots(figsize=(12, 4)) 
    contour = ax.contourf(X_plot, Y_plot, v_centered, levels=np.linspace(np.min(v_centered), np.max(v_centered), 21), cmap='viridis')
    fig.colorbar(contour, ax=ax, label='v-Velocity Magnitude [m/s]')

    # add_quiver(ax, X_plot, Y_plot, u_centered, v_centered, width=0.002)

    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_title(f'v-Velocity at t = {time:0.02f}[s]')

    ax.set_xlim([0, Lx])
    ax.set_ylim([0, Ly])
    ax.set_aspect('equal') 
    add_obstacle_outline(ax, obstacle_x_min, obstacle_x_max, obstacle_y_min, obstacle_y_max)

    plt.subplots_adjust(left=0.05, right=0.95, top=0.9, bottom=0.1) 
    plt.show()

def plot_pressure(p,):
    # x_plot = x[:u.shape[1]]
    # y_plot = y[:u.shape[0]]
    # X_plot, Y_plot = np.meshgrid(x_plot, y_plot)
    u_centered = 0.5 * (u[:-1, :-1] + u[1:, :-1])   # Average to cell centers in x-direction
    v_centered = 0.5 * (v[:, :-1] + v[:, 1:])       # Average to cell centers in y-direction
    # velocity_magnitude = np.sqrt(u_centered**2 + v_centered**2)

    x_plot = x[:u_centered.shape[1]]
    y_plot = y[:u_centered.shape[0]]
    X_plot, Y_plot = np.meshgrid(x_plot, y_plot)
    fig, ax = plt.subplots(figsize=(12, 4)) 
    contour = ax.contourf(XP, YP, p, levels=np.linspace(np.min(p), np.max(p), 21), cmap='coolwarm')
    fig.colorbar(contour, ax=ax, label='Pressure [Pa]')

    add_quiver(ax, X_plot, Y_plot, u_centered, v_centered, width=0.002)

    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_title(f'Pressure at t = {time:0.02f}[s]')

    ax.set_xlim([0, Lx])
    ax.set_ylim([0, Ly])
    ax.set_aspect('equal')  
    add_obstacle_outline(ax, obstacle_x_min, obstacle_x_max, obstacle_y_min, obstacle_y_max)

    plt.subplots_adjust(left=0.05, right=0.95, top=0.9, bottom=0.1) 
    plt.show()

def plot_scalar(c, mass_fraction):
    c_clip = np.clip(c, 0, 1)
    # x_plot = x[:u.shape[1]]
    # y_plot = y[:u.shape[0]]
    # X_plot, Y_plot = np.meshgrid(x_plot, y_plot)
    u_centered = 0.5 * (u[:-1, :-1] + u[1:, :-1])   # Average to cell centers in x-direction
    v_centered = 0.5 * (v[:, :-1] + v[:, 1:])       # Average to cell centers in y-direction

    x_plot = x[:u_centered.shape[1]]
    y_plot = y[:u_centered.shape[0]]
    X_plot, Y_plot = np.meshgrid(x_plot, y_plot)

    fig, ax = plt.subplots(figsize=(12, 4)) 
    contour = ax.contourf(XP, YP, c_clip, levels=np.linspace(0, 1, 21), cmap='turbo')
    fig.colorbar(contour, ax=ax, label='Mass Fraction')

    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_title(f'Mass Fraction at t = {time:0.02f}[s], Re-entrainment = {mass_fraction * 100:.2f}%')

    ax.set_xlim([0, Lx])
    ax.set_ylim([0, Ly])
    ax.set_aspect('equal')  
    add_obstacle_outline(ax, obstacle_x_min, obstacle_x_max, obstacle_y_min, obstacle_y_max)
    add_quiver(ax, X_plot, Y_plot, u_centered, v_centered, width=0.002)

    plt.subplots_adjust(left=0.05, right=0.95, top=0.9, bottom=0.1) 
    plt.show()
    
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

"""Grid and field initialization"""""
print(f"dx:{dx} m, dy:{dy} m")

x = np.arange(0, Lx+dx, dx) # x-grid
y = np.arange(0, Ly+dy, dy) # y-grid

print(x)
n = 0

x_u = np.arange(0, (Lx+dx), dx)          # u-velocity in x-nodes
y_u = np.arange(0, (Ly+dy)+dy, dy)       # u-velocity in y-nodes
XU,YU = np.meshgrid(x_u,y_u)
u = np.zeros([len(y_u),len(x_u)])*2

# v-vel
x_v = np.arange(0, (Lx+dx)+dx, dx)       # v-velocity, x-nodes
y_v = np.arange(0,(Ly+dy), dy)           # v-velocity, y-nodes
XV,YV = np.meshgrid(x_v,y_v)
v = np.zeros([len(y_v),len(x_v)])


# pressure
x_p = np.arange(0,(Lx+dx)+dx, dx)        # Pressure x-nodes
y_p = np.arange(0,(Ly+dy)+dy, dy)        # Pressure y-nodes
XP, YP = np.meshgrid(x_p, y_p)
p = np.zeros([len(y_p),len(x_p)])

c = np.zeros_like(p) #scalar transport
midpoint = len(c) // 2

""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
"""wind"""
alpha = 1/7
mu = 1.81*10.**-5.  #Pa-s
rho = 1.184         #kg/ms

if apply_wind_profile:
    y_padded = np.linspace(y[0] - dy / 2, y[-1] + dy / 2, len(u[:, 0]))
    u[:, 0] = u_r * (y_padded / z_r) ** alpha
else:
    u[:,0] = u_vel      # velocity at left-side boundary


""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
"""BCs"""

# Outflow BC at right boundary
# Adding ghost nodes in x_u at right boundary
x_u = np.arange(0, (Lx+dx) + dx, dx)  # Adding extra column of ghost nodes at right side boundary
XU,YU = np.meshgrid(x_u,y_u)        # Recreating meshgrid for plotting
u = np.zeros([len(y_u),len(x_u)])

if initialize_flow:
    u[:,:] = u_vel

u, v, p = apply_boundary_conditions(u,v,p,c)

print('u_shape: ',u.shape)
print('v_shape: ',v.shape)
print('P_shape: ',p.shape)


""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
"""Internal Nodes"""

"""Fractional step method, version 1"""

"""Step 1, ignore P in the momentum eqn."""
u_n = u*0
v_n = v*0
p_n = p*0
c_n = c.copy()
c_old = c.copy()

u_frac = u*2
v_frac = v*0

interior_x_u = range(1, len(x_u) - 1)
interior_y_u = range(1, len(y_u) - 1)

interior_x_v = range(1, len(x_v) - 1)
interior_y_v = range(1, len(y_v) - 1)
b = np.zeros_like(p)


num_particles = 100
particles = np.zeros((num_particles, 2))  # Each row is [x, y]

# particles x-coords
particles[:, 0] = 0  # x-coordinate is fixed

# particle y-coords
particles[:, 1] = np.linspace(2, Ly-1, num_particles)


if apply_obstacle_mask:
    obstacle_mask_u = (XU >= obstacle_x_min) & (XU <= obstacle_x_max) & (YU >= obstacle_y_min) & (YU <= obstacle_y_max)
    obstacle_mask_v = (XV >= obstacle_x_min) & (XV <= obstacle_x_max) & (YV >= obstacle_y_min) & (YV <= obstacle_y_max)
    obstacle_mask_p = (XP >= obstacle_x_min) & (XP <= obstacle_x_max) & (YP >= obstacle_y_min) & (YP <= obstacle_y_max)
    obstacle_mask_c = (XP >= obstacle_x_min) & (XP <= obstacle_x_max) & (YP >= obstacle_y_min) & (YP <= obstacle_y_max)
else:
    obstacle_mask_u = np.zeros_like(XU, dtype=bool)
    obstacle_mask_v = np.zeros_like(XV, dtype=bool)
    obstacle_mask_p = np.zeros_like(XP, dtype=bool)
    obstacle_mask_c = np.zeros_like(XP, dtype=bool)


# Ensure initial obstacle velocities are zero
u[obstacle_mask_u] = 0
v[obstacle_mask_v] = 0

# Compute maximum velocity magnitude (u or v) dynamically
u_max = np.max(np.abs(u))  # Maximum horizontal velocity
v_max = np.max(np.abs(v))  # Maximum vertical velocity
max_velocity = max(u_max, v_max)

dt = min(dx, dy) / max_velocity * 0.9  # CFL condition with safety factor
# dt = min(dtx,dty)
startdt = dt
time = 0

mass_fraction_last = 0
mass_fraction_change = np.inf

pressure_iterations = []


for n in range(nt):
    
    c_old = c.copy()

    nu = calculate_courrant(u, v, dt)
    dt = calculate_dt(u, v)
    
    p_old = p_n.copy()
    
    
    for j in interior_y_u:
        for i in interior_x_u:
            u_j_i = u[j, i]
            u_j_i_plus_1 = u[j,i+1]
            u_j_i_minus_1 = u[j,i-1]
            v_j_i = v[j, i]
            v_j_i_plus_1 = v[j,i+1]
            v_j_minus_1_i = v[j-1,i]
            v_j_minus_1_i_plus_1 = v[j-1,i+1]
            u_j_plus_1_i = u[j+1,i]
            u_j_minus_1_i = u[j-1,i]

            adv_duudx = nu/dx * (np.abs(u_j_i + u_j_i_plus_1) / 2 * (u_j_i - u_j_i_plus_1) / 2
                                   - np.abs(u_j_i_minus_1 + u_j_i) / 2 * (u_j_i_minus_1 - u_j_i) / 2)

            duudx = (1 / dx * (((u_j_i + u_j_i_plus_1) / 2) ** 2 - ((u_j_i + u_j_i_minus_1) / 2) ** 2) + adv_duudx)

            adv_dvudy = nu / dy * (np.abs(v_j_i + v_j_i_plus_1) / 2 * (u_j_i - u_j_plus_1_i) / 2
                                   - np.abs(v_j_minus_1_i + v_j_minus_1_i_plus_1) / 2 * (u_j_minus_1_i - u_j_i) / 2)

            dvudy = (1 / dy * ((v_j_i + v_j_i_plus_1) / 2 * (u_j_i + u_j_plus_1_i) / 2
                - (v_j_minus_1_i + v_j_minus_1_i_plus_1) / 2 * (u_j_minus_1_i + u_j_i) / 2)
                + adv_dvudy)

            d2udx2 = (u_j_i_plus_1 - 2 * u_j_i + u_j_i_minus_1) / dx**2
            d2udy2 = (u_j_plus_1_i - 2 * u_j_i + u_j_minus_1_i) / dy**2

            u_frac[j, i] = dt * (mu/rho * (d2udx2 + d2udy2) - (duudx + dvudy)) + u_j_i

    for j in interior_y_v:
        for i in interior_x_v:
            v_j_i = v[j, i]
            v_j_i_plus_1 = v[j,i+1]
            v_j_i_minus_1 = v[j,i-1]
            u_j_i = u[j,i]
            u_j_plus_1_i = u[j+1,i]
            u_j_i_minus_1 = u[j,i-1]
            u_j_plus_1_i_minus_1 = u[j+1, i-1]
            v_j_plus_1_i = v[j+1, i]
            v_j_minus_1_i = v[j-1, i]

            adv_duvdx = nu / dx * (np.abs(u_j_i + u_j_plus_1_i) / 2 * (v_j_i - v_j_i_plus_1) / 2
                                   - np.abs(u_j_i_minus_1 + u_j_plus_1_i_minus_1) / 2 * (v_j_i_minus_1 - v_j_i) / 2)

            duvdx = (1 / dx * ((u_j_i + u_j_plus_1_i) / 2 * (v_j_i_plus_1 + v_j_i) / 2 
                               - (u_j_i_minus_1 + u_j_plus_1_i_minus_1) / 2 * (v_j_i + v_j_i_minus_1) / 2)
                                 + adv_duvdx)

            adv_dvudy = nu / dy * (np.abs(v_j_i + v_j_plus_1_i) / 2 * (v_j_i - v_j_plus_1_i) / 2
                - np.abs(v_j_minus_1_i + v_j_i) / 2 * (v_j_minus_1_i - v_j_i) / 2)

            dvvdy = (1 / dy
                * (((v_j_i + v_j_plus_1_i) / 2) ** 2 - ((v_j_i + v_j_minus_1_i) / 2) ** 2)
                + adv_dvudy)

            d2vdx2 = (v_j_i_plus_1 - 2 * v_j_i + v_j_i_minus_1) / dx**2
            d2vdy2 = (v_j_plus_1_i - 2 * v_j_i + v_j_minus_1_i) / dy**2

            v_frac[j, i] = dt * (mu/rho * (d2vdx2 + d2vdy2) - (duvdx + dvvdy)) + v_j_i

    # b calculation (pressure Poisson)
    for j in range(1, len(y_p) - 1):
        for i in range(1, len(x_p) - 1):
            b[j, i] = (1 / dt
                       * ((u_frac[j,i] - u_frac[j,i-1]) / dx + (v_frac[j,i] - v_frac[j-1,i]) / dy))

        
    tolerance = 1e-3  # Convergence tolerance
    max_iterations = 100000  # Maximum number of iterations to prevent infinite looping
    omega = 0.7  # relaxation factor
    pressure_convergence = []
    for it in range(max_iterations):  # Maximum number of Poisson iterations per time step
        pn = p.copy()
        p[1:-1,1:-1] = (1-omega)*pn[1:-1,1:-1] + omega * (
        ((pn[1:-1,2:] + pn[1:-1,:-2])*dy**2 + 
         (pn[2:,1:-1] + pn[:-2,1:-1])*dx**2) /
        (2*(dx**2+dy**2)) -
        dx**2*dy**2/(2*(dx**2+dy**2))*b[1:-1,1:-1]
    )

        # pressure boundary conditions
        p = apply_boundary_conditions(u_frac, v_frac, p, c)[2]
        max_diff = np.max(np.abs(p - pn))  # Maximum difference at any node
        pressure_convergence.append(max_diff)

        if max_diff < tolerance:
            print(f"n: {n}, Pressure convergence, iter: {it} iterations")
            pressure_iterations.append(it)
            break
    else:
        pressure_iterations.append(max_iterations)
        print(f"n: {n}, pressure, did not converge")
    p_n = p
    
    # Update u-velocity
    for j in interior_y_u:  # Loop over interior y indices for u
        for i in interior_x_u:  # Loop over interior x indices for u
            dpdx = (p[j, i+1] - p[j, i]) / dx  # Pressure gradient in x
            u_n[j, i] = u_frac[j, i] - dt * dpdx  # Update u
    
    # Update v-velocity
    for j in interior_y_v:  # Loop over interior y indices for v
        for i in interior_x_v:  # Loop over interior x indices for v
            dpdy = (p[j+1, i] - p[j, i]) / dy  # Pressure gradient in y
            v_n[j, i] = v_frac[j, i] - dt * dpdy  # Update v
    
    u_n, v_n, p_n = apply_boundary_conditions(u_n, v_n, p_n, c_n)
    
    
    # Apply obstacle condition
    # c_n[obstacle_mask_c] = 0    
    # After the loop, copy the new values:
    # c = c_n.copy()

    
    """"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    vel_tolerance = 1e-2
    p_tolerance = 1e-2
    
    # For u and p (both 16x42)
    mask_u_p = np.ones_like(u, dtype=bool)
    
    # For v (15x42)
    mask_v = np.ones_like(v, dtype=bool)
    
    ignore_offset = 5

    # For u and p:
    mask_u_p[:ignore_offset, :] = False
    mask_u_p[-ignore_offset:, :] = False
    mask_u_p[:, :ignore_offset] = False
    mask_u_p[:, -ignore_offset:] = False
    
    mask_u_p[obstacle_mask_p] = False  # For pressure
    mask_u_p[obstacle_mask_u] = False  # For u
    
    # For v (note v has one less row in the vertical direction)
    mask_v[:ignore_offset, :] = False
    mask_v[-ignore_offset:, :] = False
    mask_v[:, :ignore_offset] = False
    mask_v[:, -ignore_offset:] = False
    
    # If you have a corresponding obstacle mask for v:
    mask_v[obstacle_mask_v] = False

    u_diff = np.abs(u_n - u)
    v_diff = np.abs(v_n - v)
    p_diff = np.abs(p_n - p_old)
    
    # Apply masks
    u_diff[~mask_u_p] = 0
    p_diff[~mask_u_p] = 0
    v_diff[~mask_v] = 0
    
    u_max_diff = np.max(u_diff)
    v_max_diff = np.max(v_diff)
    p_max_diff = np.max(p_diff)

    
    max_vel_diff = max(u_max_diff, v_max_diff)  # combined velocity criterion
    
    u = u_n.copy()
    v = v_n.copy()
    p = p_n.copy()
    c = c_n.copy()
    
    """obstacle"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    if apply_obstacle_mask:
        # inlet on right side (downwind)
        inlet_mask_u = ((XU >= obstacle_x_max - dx) & (XU <= obstacle_x_max)
                        & (YU >= obstacle_y_min+mult/2*dy) & (YU <= obstacle_y_max-mult/2*dy))
        # outlet on left side (upwind)
        outlet_mask_u = ((XU >= obstacle_x_min) & (XU <= obstacle_x_min + dx) 
                         & (YU >= obstacle_y_min+mult/2*dy) & (YU <= obstacle_y_max-mult/2*dy))
        
        inlet_mask_v = ((XV >= obstacle_x_max - dx) & (XV <= obstacle_x_max) 
                        & (YV >= obstacle_y_min+mult/2*dy) & (YV <= obstacle_y_max-mult/2*dy))
        
        outlet_mask_v = ((XV >= obstacle_x_min) & (XV <= obstacle_x_min + dx)
                         & (YV >= obstacle_y_min+mult/2*dy) & (YV <= obstacle_y_max-mult/2*dy))
        
        outlet_mask_p = ((XP >= obstacle_x_min) & (XP <= obstacle_x_min + dx) 
                         & (YP >= obstacle_y_min+mult/2*dy) & (YP <= obstacle_y_max-mult/2*dy))
        
        inlet_mask_p = ((XP >= obstacle_x_max+dx) & (XP < obstacle_x_max+2*dx)
                        & (YP >= obstacle_y_min+mult/2*dy) & (YP < obstacle_y_max-mult/2*dy))

    
        # Zero out velocities in the obstacle region
        u_n[obstacle_mask_u] = 0
        v_n[obstacle_mask_v] = 0
    
        # Set inlet velocities
        u[inlet_mask_u] = u_inlet
        v[inlet_mask_v] = 0  # Assume no vertical velocity at the inlet
    
        # Set outlet velocities
        u[outlet_mask_u] = u_outlet
        v[outlet_mask_v] = 0  # Assume no vertical velocity at the outlet
        
        c_n[obstacle_mask_c] = 0
        
        if n > n_switch:
            c_n[outlet_mask_p] = 1.0
        
        c = c_n.copy()
    
        inlet_cells = np.where(inlet_mask_p)  # This gives arrays of j, i indices where mask is True
        int_tracer = 0.0
        int_flow = 0.0
        for j, i in zip(inlet_cells[0], inlet_cells[1]):
            cell_c = c[j, i]  # assuming the domain is to the left of the inlet face
            cell_u = u[j, i+1]
        
            flux_tracer = rho * cell_c * cell_u * dy
            flux_air = rho * cell_u * dy
        
            int_tracer += flux_tracer
            int_flow += flux_air
        
        if int_flow != 0:
            mass_fraction = int_tracer / int_flow
            mass_fraction_change = abs(mass_fraction - mass_fraction_last)
            # else:
            #     mass_fraction_change = None
        else:
            mass_fraction = 0.0
        c_n[inlet_mask_p] = 0.0
        c = c_n.copy()

        print(f"Mass fraction at building inlet: {mass_fraction}")
        
    for j in range(1, len(y_p)-2):
        for i in range(1, len(x_p)-1):
            
            # left face (between i and i-1):
            if u[j,i] > 0:
                c_left = c_old[j,i-1]
            else:
                c_left = c_old[j,i]
            flux_left = u[j,i] * c_left
            
            # right face (between i and i+1):
            if u[j,i+1] > 0:
                c_right = c_old[j,i]
            else:
                c_right = c_old[j,i+1]
            flux_right = u[j,i+1] * c_right
    
            ducdx = (flux_right - flux_left) / dx
    
            # bottom face flux (between j and j+1):
            if v[j,i] > 0:
                c_bottom = c_old[j-1,i]  
            else:
                c_bottom = c_old[j,i]    
            flux_bottom = v[j,i] * c_bottom
    
            # Top face flux (between j and j+1):
            if v[j+1,i] > 0:
                c_top = c_old[j,i]       
            else:
                c_top = c_old[j+1,i]
            flux_top = v[j+1,i] * c_top
    
            dvcdy = (flux_top - flux_bottom) / dy
    
            d2cdx2 = (c_old[j,i+1] - 2*c_old[j,i] + c_old[j,i-1]) / dx**2
            d2cdy2 = (c_old[j+1,i] - 2*c_old[j,i] + c_old[j-1,i]) / dy**2
    
            c_n[j,i] = (c_old[j,i] - dt*(ducdx + dvcdy) + dt*diffusion_coefficient*(d2cdx2 + d2cdy2))
    
    vol_flux = (int_flow+int_tracer)/rho
    print(f'vol_flux from outlet = {vol_flux}')


    if n % 5 == 0:
        plot_velocity(u_n,v_n)
        plot_pressure(p)
        plot_scalar(c,mass_fraction)
            
    time += dt

    selected_indices = np.where(inlet_mask_p)
    num_cells_selected = len(selected_indices[0])    
    
    print(f"n: {n}, p_max_diff: {p_max_diff:.5f}, max_vel_diff: {max_vel_diff:.5f}, mass_frac_diff: {mass_fraction_change}")

    mass_frac_tolerance = 1e-3
    if max_vel_diff < vel_tolerance and p_max_diff < p_tolerance and mass_fraction_change < mass_frac_tolerance and n > n_switch:
        # n_switch = n+1
        n_converge = n
        break
    
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
"""Plotting"""

# wind profile
y_padded = np.linspace(y[0] - dy / 2, y[-1] + dy / 2, len(u[:, 0]))
# u_z = u_r * (y_padded/z_r)**alpha

y_padded_corrected = np.maximum(y_padded, 1e-6)  # small positive floor
u[:, 0] = u_r * (y_padded_corrected / z_r)**alpha
height_bldg = 20 # meters

u[-1,:] = u[-2,:] # just for plotting
y_fine = np.linspace(y[0], y[-1], 1000)  
u_z_fine = -u_r * (y_fine / z_r) ** alpha  

height_bldg = 20  # meters

plt.figure()
plt.plot(u_z_fine, y_fine, label='Analytical Wind Profile', color='blue', linewidth=2)  # Smooth analytical
plt.plot(-u[:, 0], y_padded, label='Discretized Wind Profile', color='orange', linestyle='--', linewidth=2)  # Discretized
plt.axhline(y=height_bldg, color='gray', linestyle='--', linewidth=1, label=f'Bldg. Height ({height_bldg} m)')

plt.title("Wind Profile")
plt.xlabel("Velocity [m/s]")
plt.ylabel("Height [m]")
plt.gca().yaxis.set_label_position("right")
plt.gca().yaxis.tick_right()
plt.gca().set_aspect(0.25, adjustable='box')  
plt.xlim(right=0)
plt.legend()
plt.grid(True)
plt.show()



plt.figure(figsize=(12, 8))
# skip = 1

XP = XP-dx/2
YP = YP-dy/2
plt.scatter(XP, YP, color='green', marker='+', s=45, label='p Grid points')
# plt.scatter(XP[::skip], YP[::skip], color='green', marker='+', s=45, label='p Grid points')


YU = YU-dy/2
plt.scatter(XU, YU, color='blue', marker='>', s=45,label='u - Grid points')
# plt.scatter(XU[::skip], YU[::skip], color='blue', marker='>', s=45,label='u - Grid points')

XV = XV-dx/2
plt.scatter(XV, YV, color='red', marker='^', s=45,label='v - Grid points')
# plt.scatter(XV[::skip], YV[::skip], color='red', marker='^', s=45,label='v - Grid points')


plt.plot([0, Lx, Lx, 0, 0], [0, 0, Ly, Ly, 0], color='black', linestyle='-', linewidth=2, label="Domain Boundary")

plt.xlabel('x (m)')
plt.ylabel('y (m)')
plt.title('Grid')
plt.axis('equal')  # Ensure equal aspect ratio
plt.grid(False)  # Turn off default matplotlib grid
plt.xlim([0-1, Lx+1])
plt.ylim([0-1-2*dy,Ly+1+2*dy])
plt.axis('equal')  # Ensure equal aspect ratio
plt.legend()
plt.show()

enddt = dt

cmax = np.max(c_n)
cmax_loc = np.argmax(c_n)
rowmax, colmax = np.unravel_index(cmax_loc,c.shape)

plt.figure(figsize=(8, 6))
plt.plot(range(len(pressure_convergence)), pressure_convergence)
plt.yscale('log')
plt.xlabel("Iteration")
plt.ylabel("Maximum Pressure Difference")
plt.title("Pressure Convergence")
plt.grid(True)
plt.legend()
plt.savefig("Pressure Convergence Plot")
plt.show()

plt.figure(figsize=(8, 6))
plt.plot(range(len(pressure_iterations)), pressure_iterations)
plt.xlabel("Timestep (n)")
plt.ylabel("Pressure Solver Iterations")
plt.title("Pressure Solver Convergence per Timestep")
plt.grid(True)
plt.legend()
plt.savefig("Pressure Poisson Solver Convergence")
plt.show()

print(pressure_iterations)
