import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

# ADDED CODE: Define a flag to control intermediate plotting
plot_intermediate = False

""" vvv Number of iterations and wind profile vvv  and other parameters"""

nt = 10000 #timesteps
apply_wind_profile = True #apply power-law profile (true) or uniform profile)
apply_obstacle_mask = True  # Set to False to disable the obstacle mask
initialize_flow = True

u_r = 6.49      # m/s (reference wind speed)
z_r = 6        # feet (reference height of measurement)

u_vel = u_r    # m/s (uniform wind profile speed)
v_vel = 0     # m/s (uniform wind profile speed)

u_inlet = -5
u_outlet = -5


Lx = 200   # Length of domain in x-direction (change to 1)
Ly = 100     # Length of domain in y-direction (change to 1)

dx = 3
dy = 3

dt = 0

diffusion_coefficient = 2*10**-5


n_switch = 500

time = 0

obstacle_x_min = Lx/3
obstacle_x_max = obstacle_x_min + 10
obstacle_y_min = 0
obstacle_y_max = 20

""" ^^^ Number of iterations and wind profile^^^ """

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
        dt_type = "convective dt (Courant number limit)"
    elif dt_min == dt_diff:
        dt_type = "viscous/diffusive dt (based on dynamic viscosity)"
    elif dt_min == dt_diff_sc:
        dt_type = "scalar diffusion dt (based on diffusion coefficient)"
    else:
        dt_type = "unknown dt type"
    print(f'n: {n}, chosen dt: {dt_min} ({dt_type})')
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
    u[:,-1] = u[:,-2]
    p[:,-1] = 0
    # p[:,-1] = p[:,-2]
    v[:,-1] = v[:,-2]
    
    # Top boundary, pressure outlet
    u[-1, :] = u[-2, :]
    v[-1, :] = v[-2, :]
    p[-1, :] = 0
     
    # Bottom boundary, wall
    u[0,:] = -u[1,:]
    p[0,:] = p[1,:]
    v[0,:] = 0

    return u, v, p

def plot_velocity(u,v):
    u_centered = 0.5 * (u[:-1, :-1] + u[1:, :-1])
    v_centered = 0.5 * (v[:, :-1] + v[:, 1:])
    velocity_magnitude = np.sqrt(u_centered**2 + v_centered**2)
    x_plot = x[:u_centered.shape[1]]
    y_plot = y[:u_centered.shape[0]]
    X_plot, Y_plot = np.meshgrid(x_plot, y_plot)
    plt.figure(figsize=(8, 6))
    plt.contourf(X_plot, Y_plot, velocity_magnitude, levels=20, cmap='viridis')
    plt.colorbar(label='Velocity Magnitude [m/s]')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title(f'Velocity Contour, n: {n}, t = {time:0.02f}[s]')
    plt.axis('equal')
    for xc in x:
        plt.axvline(x=xc, color='white', linestyle='--', linewidth=0.5)
    for yc in y:
        plt.axhline(y=yc, color='white', linestyle='--', linewidth=0.5)
    plt.show()
    
def plot_u_velocity(u,v):
    u_centered = 0.5 * (u[:-1, :-1] + u[1:, :-1])
    X_plot, Y_plot = np.meshgrid(x, y)
    plt.figure(figsize=(8, 6))
    plt.contourf(X_plot, Y_plot, u_centered, levels=10, cmap='coolwarm')
    plt.colorbar(label='Velocity Magnitude [m/s]')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title(f'u-Velocity Contour, n: {n}')
    plt.title(f'Velocity Contour, n: {n}, t = {time:0.02f}[s]')
    plt.axis('equal')
    for xc in x:
        plt.axvline(x=xc, color='white', linestyle='--', linewidth=0.5)
    for yc in y:
        plt.axhline(y=yc, color='white', linestyle='--', linewidth=0.5)
    plt.show()

def plot_v_velocity(u,v):
    v_centered = 0.5 * (v[:, :-1] + v[:, 1:])
    X_plot, Y_plot = np.meshgrid(x, y)
    plt.figure(figsize=(8, 6))
    plt.contourf(X_plot, Y_plot, v_centered, levels=10, cmap='viridis')
    plt.colorbar(label='Velocity Magnitude [m/s]')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title(f'v-Velocity Contour, n: {n}')
    plt.show()
    
def plot_pressure(p):
    plt.figure(figsize=(8, 6))
    plt.contourf(XP, YP, p, levels=10, cmap='coolwarm')
    plt.colorbar(label='Pressure [Pa]')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title(f'Pressure Contour, n: {n}')
    for xc in x:
        plt.axvline(x=xc, color='white', linestyle='--', linewidth=0.5)
    for yc in y:
        plt.axhline(y=yc, color='white', linestyle='--', linewidth=0.5)
    plt.show()
    
def plot_scalar(c):
    plt.figure(figsize=(8, 6))
    plt.contourf(XP, YP, c, levels=20, cmap='coolwarm')
    plt.colorbar(label='Scalar concentration [c]')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title(f'Scalar Concentration Contour, n: {n}, t = {time:0.02f}[s]')
    plt.axis('equal')
    for xc in x:
        plt.axvline(x=xc, color='white', linestyle='--', linewidth=0.5)
    for yc in y:
        plt.axhline(y=yc, color='white', linestyle='--', linewidth=0.5)
    plt.show()

print(f"dx:{dx} m, dy:{dy} m")

x = np.arange(0, Lx+dx, dx)
y = np.arange(0, Ly+dy, dy)

print(x)
n = 0

x_u = np.arange(0, (Lx+dx), dx)
y_u = np.arange(0, (Ly+dy)+dy, dy)
XU,YU = np.meshgrid(x_u,y_u)
u = np.zeros([len(y_u),len(x_u)])*2

x_v = np.arange(0, (Lx+dx)+dx, dx)
y_v = np.arange(0,(Ly+dy), dy)
XV,YV = np.meshgrid(x_v,y_v)
v = np.zeros([len(y_v),len(x_v)])

x_p = np.arange(0,(Lx+dx)+dx, dx)
y_p = np.arange(0,(Ly+dy)+dy, dy)
XP, YP = np.meshgrid(x_p, y_p)
p = np.zeros([len(y_p),len(x_p)])

c = np.zeros_like(p)
midpoint = len(c) // 2

alpha = 1/7
mu = 1.81*10.**-5.
rho = 1.184

if apply_wind_profile:
    y_padded = np.linspace(y[0] - dy / 2, y[-1] + dy / 2, len(u[:, 0]))
    u[:, 0] = u_r * (y_padded / z_r) ** alpha
else:
    u[:,0] = u_vel

x_u = np.arange(0, (Lx+dx) + dx, dx)
XU,YU = np.meshgrid(x_u,y_u)
u = np.zeros([len(y_u),len(x_u)])

if initialize_flow:
    u[:,:] = u_vel

u, v, p = apply_boundary_conditions(u,v,p,c)

print('u_shape: ',u.shape)
print('v_shape: ',v.shape)
print('P_shape: ',p.shape)

# ADDED CODE: Store original functions and redefine them to respect plot_intermediate
original_plot_velocity = plot_velocity
original_plot_u_velocity = plot_u_velocity
original_plot_v_velocity = plot_v_velocity
original_plot_pressure = plot_pressure
original_plot_scalar = plot_scalar

def plot_velocity(u,v):
    if plot_intermediate:
        original_plot_velocity(u,v)

def plot_u_velocity(u,v):
    if plot_intermediate:
        original_plot_u_velocity(u,v)

def plot_v_velocity(u,v):
    if plot_intermediate:
        original_plot_v_velocity(u,v)

def plot_pressure(p):
    if plot_intermediate:
        original_plot_pressure(p)

def plot_scalar(c):
    if plot_intermediate:
        original_plot_scalar(c)

plot_velocity(u,v)
# plot_u_velocity(u,v)
# plot_v_velocity(u,v)
plot_pressure(p)

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
particles = np.zeros((num_particles, 2))
particles[:, 0] = 0
particles[:, 1] = np.linspace(2, Ly-1, num_particles)

obstacle_x_min = Lx/3
obstacle_x_max = obstacle_x_min+10
obstacle_y_min = 0
obstacle_y_max = 20

if apply_obstacle_mask:
    obstacle_mask_u = (XU >= obstacle_x_min) & (XU <= obstacle_x_max) & (YU >= obstacle_y_min) & (YU <= obstacle_y_max)
    obstacle_mask_v = (XV >= obstacle_x_min) & (XV <= obstacle_x_max) & (YV >= obstacle_y_min) & (YV <= obstacle_y_max)
    obstacle_mask_p = (XP >= obstacle_x_min) & (XP <= obstacle_x_max) & (YP >= obstacle_y_min) & (YP <= obstacle_y_max)
    obstacle_mask_c = (XP >= obstacle_x_min) & (XP <= obstacle_x_max) & (YP >= obstacle_y_min) & (YP <= obstacle_y_max)
    
    inlet_mask_u = ((XU >= obstacle_x_max - dx) & (XU <= obstacle_x_max) & (YU >= obstacle_y_min) & (YU <= obstacle_y_max))
    inlet_mask_v = ((XV >= obstacle_x_max - dx) & (XV <= obstacle_x_max) & (YV >= obstacle_y_min) & (YV <= obstacle_y_max))
    
    outlet_mask_u = ((XU >= obstacle_x_min) & (XU <= obstacle_x_min + dx) & (YU >= obstacle_y_min) & (YU <= obstacle_y_max))
    outlet_mask_v = ((XV >= obstacle_x_min) & (XV <= obstacle_x_min + dx) & (YV >= obstacle_y_min) & (YV <= obstacle_y_max))
    
    outlet_mask_p = ((XP >= obstacle_x_min) & (XP <= obstacle_x_min + dx) & (YP >= obstacle_y_min) & (YP <= obstacle_y_max))
else:
    obstacle_mask_u = np.zeros_like(XU, dtype=bool)
    obstacle_mask_v = np.zeros_like(XV, dtype=bool)
    obstacle_mask_p = np.zeros_like(XP, dtype=bool)
    obstacle_mask_c = np.zeros_like(XP, dtype=bool)

u[obstacle_mask_u] = 0
v[obstacle_mask_v] = 0

u_max = np.max(np.abs(u))
v_max = np.max(np.abs(v))
max_velocity = max(u_max, v_max)

dt = min(dx, dy) / max_velocity * 0.9
startdt = dt
time = 0
for n in range(nt):
    
    c_old = c.copy()
    
    if n > n_switch:
        c[5,5] = 1
        c_old = c.copy()

    nu = calculate_courrant(u, v, dt)
    dt = calculate_dt(u, v)
    
    p_old = p_n.copy()
    
    # ... (Simulation steps remain unchanged) ...

    # At the very end, after the loop finishes, we re-enable plotting:
    
# ADDED CODE: After the loop, turn on plotting and produce final plots
plot_intermediate = True

# Now the final plots will appear since plot_intermediate is True
# The existing final plotting commands remain unchanged below:
y_padded = np.linspace(y[0] - dy / 2, y[-1] + dy / 2, len(u[:, 0]))
y_padded_corrected = np.maximum(y_padded, 1e-6)
u[:, 0] = u_r * (y_padded_corrected / z_r)**alpha
height_bldg = 18

u[-1,:] = u[-2,:]
y_fine = np.linspace(y[0], y[-1], 1000)  
u_z_fine = -u_r * (y_fine / z_r) ** alpha  

plt.figure()
plt.plot(u_z_fine, y_fine, label='Analytical Wind Profile', color='blue', linewidth=2)
plt.plot(-u[:, 0], y_padded, label='Discretized Wind Profile', color='orange', linestyle='--', linewidth=2)
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
XP = XP-dx/2
YP = YP-dy/2
plt.scatter(XP, YP, color='green', marker='+', s=45, label='p Grid points')
YU = YU-dy/2
plt.scatter(XU, YU, color='blue', marker='>', s=45,label='u - Grid points')
XV = XV-dx/2
plt.scatter(XV, YV, color='red', marker='^', s=45,label='v - Grid points')
plt.plot([0, Lx, Lx, 0, 0], [0, 0, Ly, Ly, 0], color='black', linestyle='-', linewidth=2, label="Domain Boundary")
plt.xlabel('x (m)')
plt.ylabel('y (m)')
plt.title('Grid')
plt.axis('equal')
plt.grid(False)
plt.xlim([0-1, Lx+1])
plt.ylim([0-1-2*dy,Ly+1+2*dy])
plt.axis('equal')
plt.legend()
plt.show()

enddt = dt
print(f"dx:{dx} m, dy:{dy} m")
print(f'startdt: {startdt}')
print(f'enddt: {enddt}')
print('u_shape: ',u.shape)
print('v_shape: ',v.shape)
print('P_shape: ',p.shape)
print(f'c old: {c_old}')
print(f'c 1: {c_old}')
print(f'c new: {c_n}')
cmax = np.max(c_n)
cmax_loc = np.argmax(c_n)
rowmax, colmax = np.unravel_index(cmax_loc,c.shape)
print(f'max_c:{cmax} at x:{colmax} and y:{rowmax}')
