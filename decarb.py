#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 19 00:09:06 2024

@author: tejjolly
"""

import numpy as np
import matplotlib.pyplot as plt


""" vvv Number of iterations and wind profile vvv  and other parameters"""

nt = 1000 #timesteps
apply_wind_profile = True #apply power-law profile (true) or uniform profile)

u_vel = 2      # m/s (uniform wind profile speed)
u_r = 8.49      # m/s (reference wind speed)
z_r = 6        # feet (reference height of measurement)


Lx = 60    # Length of domain in x-direction (change to 1)
Ly = 20      # Length of domain in y-direction (change to 1)

dx = Lx/(180)
dy = Ly/(60)

""" ^^^ Number of iterations and wind profile^^^ """


def calculate_dt(u,v):
    u_max = np.max(np.abs(u))
    v_max = np.max(np.abs(v))
    vel_dlength_max = max(u_max/dx,v_max/dy)
    courrant_num_adjuster = 0.3
    dt = 1/vel_dlength_max
    print(f'n: {n}, calculated dt: {dt}')
    
    # Viscous/diffusive dt?
    dt_diff = dx**2/(mu/rho)
    if dt_diff < dt:
        print('\n\n\n\ndt_diff used!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n\n')
    return courrant_num_adjuster*min(dt,dt_diff)

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
        # v[1,:] = u_vel      # velocity at left-side boundary
    
    # Left boundary, inflow BC (continued)
    p[:,0] = p[:,1]     # dP/dx|x=0 = 0
    v[:,0] = 0         # v|x=0 = 0
    
    # Right boundary, outflow BC
    u[:,-1] = u[:,-2]   # u(M+1,1) = u(M,1)
    p[:,-2] = 0         # P|x=Lx = 0 (0 gauge outlet)
    p[:,-1] = -p[:,-2]  # P(M+1,1) = -P(M,1)
    v[:,-1] = v[:,-2]   # dv/dx = 0
    
    # Top boundary, wall
    u[-1,:] = -u[-2,:]  # no slip at top wall
    p[-1,:] = p[-2,:]   # top wall, horizontal wall dP/dy = 0
    v[-1,:] = 0         # impermeable wall, v = 0
     
    
    # Bottom boundary, wall
    u[0,:] = -u[1,:]    # no slip at bottom wall
    p[0,:] = p[1,:]     # bottom wall, horizontal wall dP/dy = 0 
    v[0,:] = 0          # bottom wall (ground), v = 0
    
    # Top boundary, wall
    # c[-1,:] = -c[-2,:]  # no slip at top wall
    # c[0,:] = -c[1,:]    # no slip at bottom wall
    # c[:,0] = c[:,1]     # dP/dx|x=0 = 0
    # c[:,-1] = c[:,-2]   # u(M+1,1) = u(M,1)
    
    return u, v, p
    
def plot_velocity(u,v):
    u_centered = 0.5 * (u[:-1, :-1] + u[1:, :-1])   # Average to cell centers in x-direction
    v_centered = 0.5 * (v[:, :-1] + v[:, 1:])       # Average to cell centers in y-direction
    velocity_magnitude = np.sqrt(u_centered**2 + v_centered**2)
    
    # # New meshgrid for plotting
    # x_plot = np.arange(u_centered.shape[1])
    # y_plot = np.arange(u_centered.shape[0])
    # X_plot, Y_plot = np.meshgrid(x_plot, y_plot)
    
    x_plot = x[:u_centered.shape[1]]
    y_plot = y[:u_centered.shape[0]]
    X_plot, Y_plot = np.meshgrid(x_plot, y_plot)

    # velocity contour
    plt.figure(figsize=(8, 6))
    # plt.figure()
    plt.contourf(X_plot, Y_plot, velocity_magnitude, levels=20, cmap='viridis')  # Contour plot
    plt.colorbar(label='Velocity Magnitude [m/s]')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title(f'Velocity Contour, n: {n}')
    plt.axis('equal')
    for xc in x:
        plt.axvline(x=xc, color='white', linestyle='--', linewidth=0.5)
    for yc in y:
        plt.axhline(y=yc, color='white', linestyle='--', linewidth=0.5)
    plt.show()
    
def plot_u_velocity(u,v):
    u_centered = 0.5 * (u[:-1, :-1] + u[1:, :-1])   # Average to cell centers in x-direction
    # v_centered = 0.5 * (v[:, :-1] + v[:, 1:])       # Average to cell centers in y-direction
    # velocity_magnitude = np.sqrt(u_centered**2 + v_centered**2)
    
    # New meshgrid for plotting
    # x_plot = np.arange(u_centered.shape[1])
    # y_plot = np.arange(u_centered.shape[0])
    X_plot, Y_plot = np.meshgrid(x, y)
    
    # u-velocity contour
    plt.figure(figsize=(8, 6))
    plt.contourf(X_plot, Y_plot, u_centered, levels=10, cmap='viridis')  # Contour plot
    plt.colorbar(label='Velocity Magnitude [m/s]')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title(f'u-Velocity Contour, n: {n}')
    plt.show()
    
def plot_v_velocity(u,v):
    # u_centered = 0.5 * (u[:-1, :-1] + u[1:, :-1])   # Average to cell centers in x-direction
    v_centered = 0.5 * (v[:, :-1] + v[:, 1:])       # Average to cell centers in y-direction
    # velocity_magnitude = np.sqrt(u_centered**2 + v_centered**2)
    
    # New meshgrid for plotting
    # x_plot = np.arange(u_centered.shape[1])
    # y_plot = np.arange(u_centered.shape[0])
    X_plot, Y_plot = np.meshgrid(x, y)
    
    # v-velocity contour
    plt.figure(figsize=(8, 6))
    plt.contourf(X_plot, Y_plot, v_centered, levels=10, cmap='viridis')  # Contour plot
    plt.colorbar(label='Velocity Magnitude [m/s]')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title(f'v-Velocity Contour, n: {n}')
    plt.show()
    
def plot_pressure(p):
    # pressure contour
    plt.figure(figsize=(8, 6))
    plt.contourf(XP, YP, p, levels=10, cmap='coolwarm')  # Contour plot
    plt.colorbar(label='Pressure [Pa]')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title(f'Pressure Contour, n: {n}')
    plt.show()
    
def plot_scalar(c):
    # pressure contour
    # plt.figure(figsize=(8, 6))
    # plt.contourf(XP, YP, p, levels=10, cmap='coolwarm')  # Contour plot
    # plt.colorbar(label='Pressure [Pa]')
    # plt.xlabel('x')
    # plt.ylabel('y')
    # plt.title(f'Pressure Contour, n: {n}')
    # plt.show()
    
    # scalar contour
    vmin = 0
    vmax = 2
    # c = np.clip(c, vmin, vmax)
    plt.figure(figsize=(8, 6))
    plt.contourf(XP, YP, c, levels=20, cmap='coolwarm') # Contour plot
    # plt.contourf(XP, YP, c, levels=np.linspace(vmin, vmax, 21), cmap='coolwarm')  # Contour plot
    plt.colorbar(label='Scalar concentration [c]')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title(f'Scalar Concentration Contour, n: {n}')
    plt.axis('equal')
    
    for xc in x:
        plt.axvline(x=xc, color='white', linestyle='--', linewidth=0.5)
    for yc in y:
        plt.axhline(y=yc, color='white', linestyle='--', linewidth=0.5)
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
#TODO
#c[4,4] = 10

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
    
dt = dx/u_vel * .9  #seconds
startdt = dt

""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
"""BCs"""

# Outflow BC at right boundary
# Adding ghost nodes in x_u at right boundary
x_u = np.arange(0, (Lx+dx) + dx, dx)  # Adding extra column of ghost nodes at right side boundary
XU,YU = np.meshgrid(x_u,y_u)        # Recreating meshgrid for plotting
u = np.zeros([len(y_u),len(x_u)])

## Adding ghost nodes in x_p at right boundary
# x_p = np.arange(0,(Lx+dx)+dx + dx, dx)   # Adding extra column of ghost nodes at right side boundary
# XP, YP = np.meshgrid(x_p, y_p)      # Recreating meshgrid for plotting
# p = np.zeros([len(y_p),len(x_p)])

u, v, p = apply_boundary_conditions(u,v,p,c)

print('u_shape: ',u.shape)
print('v_shape: ',v.shape)
print('P_shape: ',p.shape)

plot_velocity(u,v)
# plot_u_velocity(u,v)
# plot_v_velocity(u,v)
plot_pressure(p)


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


num_particles = 20
particles = np.zeros((num_particles, 2))  # Each row is [x, y]

# particles x-coords
particles[:, 0] = 0  # x-coordinate is fixed

# particle y-coords
particles[:, 1] = np.linspace(1, Ly-1, num_particles)


for n in range(nt):
        
    if n == 30:
        c[5,5] = 1
        c_old = c.copy()
        # print(f'start c at n:{n}: ')
        
    # if n == (nt//2 + 6):
    #     # c[4,4] = 10
    #     c_1 = c_n.copy()

    
    nu = calculate_courrant(u, v, dt)
    dt = calculate_dt(u, v)
    
    p_old = p_n.copy()
    
    # print(f'n:{n}')
    # print(f'u_max:{u_max}')
    # print(f'v_max:{v_max}')
    # print(f'nu_u:{nu_u}')
    # print(f'nu_v:{nu_v}')
    # print(f'u_dlength_max:{u_max/dx}')
    # print(f'v_dlength_max:{v_max/dy}')
    # print(f'new dt: {dt}')
    
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

        
    tolerance = 1e-4  # Convergence tolerance
    max_iterations = 100000  # Maximum number of iterations to prevent infinite looping

    for it in range(max_iterations):  # Maximum number of Poisson iterations per time step
        pn = p.copy()
        p[1:-1,1:-1] = (((pn[1:-1,2:] + pn[1:-1,0:-2]) * dy**2 +
                         (pn[2:,1:-1] + pn[0:-2,1:-1]) * dx**2) /
                        (2 * (dx**2 + dy**2)) -
                        dx**2 * dy**2 / (2 * (dx**2 + dy**2)) * b[1:-1,1:-1])

        # pressure boundary conditions
        p = apply_boundary_conditions(u_frac, v_frac, p, c)[2]
        max_diff = np.max(np.abs(p - pn))  # Maximum difference at any node
        if max_diff < tolerance:
            print(f"n: {n}, Pressure convergence, iter: {it} iterations")
            break
    else:
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
    
    """scalar transport"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    for j in range(1, len(y_v) - 1):
        for i in range(1, len(x_u) - 1):
            ducdx = (1/dx * (u_n[j,i] * (c[j,i+1] + c[j,i])/2)
                    - u_n[j,i-1] * (c[j,i] + c[j,i-1])/2)
            dvcdy = (1/dy * (v_n[j+1,i] * (c[j+1,i] + c[j,i])/2)
                    - v_n[j-1,i] * (c[j,i] + c[j-1,i])/2)
            c_n[j,i] = dt * -(ducdx + dvcdy) + c[j,i]#
            # c_n[j,i] = c[j,i] - dt * (ducdx + dvcdy)

    
    
    """"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    vel_tolerance = 1e-6
    p_tolerance = 1e-4
    
    u_max_diff = np.max(np.abs(u_n - u))
    v_max_diff = np.max(np.abs(v_n - v))
    p_max_diff = np.max(np.abs(p_n - p_old))
    
    max_vel_diff = max(u_max_diff, v_max_diff)  # combined velocity criterion
    
    if max_vel_diff < vel_tolerance and p_max_diff < p_tolerance and n!= 0:
        print(f"Steady state, n: {n}")
        break
    # print(f'c {n}: {c}')

    u = u_n.copy()
    v = v_n.copy()
    p = p_n.copy()
    c = c_n.copy()
    
    # if n % 1 == 0:
        # plot_scalar(u_n,v_n)
        # plot_scalar(c)
    
    for k in range(num_particles):
        particle_x, particle_y = particles[k]
    
        # ensures particles remain within domain
        i = min(max(int(particle_x / dx), 0), len(x_u) - 2) 
        j = min(max(int(particle_y / dy), 0), len(y_u) - 2) 
    
        # interp for u velocity
        u_interp = (u[j,i] * (1 - (particle_x % dx) / dx) * (1 - (particle_y % dy) / dy) +
                    u[j,i+1] * ((particle_x % dx) / dx) * (1 - (particle_y % dy) / dy) +
                    u[j+1,i] * (1 - (particle_x % dx) / dx) * ((particle_y % dy) / dy) +
                    u[j+1,i+1] * ((particle_x % dx) / dx) * ((particle_y % dy) / dy))
    
        # interp for v velocity
        v_interp = (v[j,i] * (1 - (particle_x % dx) / dx) * (1 - (particle_y % dy) / dy) +
                    v[j,i+1] * ((particle_x % dx) / dx) * (1 - (particle_y % dy) / dy) +
                    v[j+1,i] * (1 - (particle_x % dx) / dx) * ((particle_y % dy) / dy) +
                    v[j+1,i+1] * ((particle_x % dx) / dx) * ((particle_y % dy) / dy))
    
        # Update particle position
        particle_x += u_interp * dt
        particle_y += v_interp * dt
    
        # Reflect particle if it hits the boundary
        particle_x = max(0, min(Lx, particle_x))
        particle_y = max(0, min(Ly, particle_y))
    
        # Store updated position
        particles[k] = [particle_x, particle_y]

        # Extract x and y positions of all particles
    particle_x_positions = particles[:, 0]
    particle_y_positions = particles[:, 1]
    
    # Plot particle positions
    plt.figure(figsize=(8, 6))
    plt.scatter(particle_x_positions, particle_y_positions, color='red', label='Particles')  # Adjust y
    plt.xlim(0, Lx)
    plt.ylim(0, Ly)
    plt.xlabel('x [m]')
    plt.ylabel('y [m]')
    plt.title('Particle Trajectories')
    plt.gca().set_aspect('equal', adjustable='box')
    plt.legend()
    # plt.grid(True, which='both', linestyle='--', color='gray', linewidth=0.5)
    for xc in x:
        plt.axvline(x=xc, color='black', linestyle='--', linewidth=0.5)
    for yc in y:
        plt.axhline(y=yc, color='black', linestyle='--', linewidth=0.5)
    plt.show()



""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
"""Plotting"""

# wind profile
y_padded = np.linspace(y[0] - dy / 2, y[-1] + dy / 2, len(u[:, 0]))
# u_z = u_r * (y_padded/z_r)**alpha
height_bldg = 18 # meters

u[-1,:] = u[-2,:] # just for plotting
y_fine = np.linspace(y[0], y[-1], 1000)  
u_z_fine = -u_r * (y_fine / z_r) ** alpha  

height_bldg = 18  # meters

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