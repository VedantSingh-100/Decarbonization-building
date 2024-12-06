#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 19 00:09:06 2024

@author: tejjolly
"""

import numpy as np
import matplotlib.pyplot as plt


""" vvv Number of iterations and wind profile vvv """

nt = 1000
apply_wind_profile = True

""" ^^^ Number of iterations and wind profile^^^ """


def calculate_dt(u,v):
    u_max = np.max(np.abs(u))
    v_max = np.max(np.abs(v))
    vel_dlength_max = max(u_max/dx,v_max/dy)
    courrant_num_adjuster = 0.3
    dt = 1/vel_dlength_max * courrant_num_adjuster
    print(f'n: {n}, calculated dt: {dt}')
    
    # Viscous/diffusive dt?
    dt_diff = dx**2/(mu/rho)
    if dt_diff < dt:
        print('\n\n\n\ndt_diff used!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n\n')
    return min(dt,dt_diff)

def calculate_courrant(u,v,dt):
    u_max = np.max(np.abs(u))
    v_max = np.max(np.abs(v))
    nu_u = u_max*dt/dx
    nu_v = v_max*dt/dy
    nu = np.max([nu_u,nu_v])
    return nu

def apply_boundary_conditions(u, v, p):
    # Left boundary, inflow BC
    if apply_wind_profile:
        y_padded = np.linspace(y[0] - dy / 2, y[-1] + dy / 2, len(u[:, 0]))
        u[:, 0] = u_r * (y_padded / z_r) ** alpha
    else:
        u[:,0] = u_vel      # velocity at left-side boundary
    
    # Left boundary, inflow BC (continued)
    p[:,0] = p[:,1]     # dP/dx|x=0 = 0
    v[:,0] = 0          # v|x=0 = 0
    
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
    
    
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

"""Grid"""""
u_vel = 10 #m/s

Lx = 50    # Length of domain in x-direction (change to 1)
Ly = 25      # Length of domain in y-direction (change to 1)

dx = Lx/(20)
dy = Ly/(10)

print(f"dx:{dx} m, dy:{dy} m")

x = np.arange(0, Lx+dx, dx) # x-grid
y = np.arange(0, Ly+dy, dy) # y-grid

print(x)
n = 0

x_u = np.arange(0, (Lx+dx), dx)          # u-velocity in x-nodes
y_u = np.arange(0, (Ly+dy)+dy, dy)       # u-velocity in y-nodes
XU,YU = np.meshgrid(x_u,y_u)
u = np.zeros([len(y_u),len(x_u)])

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

""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
"""wind"""
u_r = 8.49      # m/s (reference wind speed)
z_r = 10        # meters (reference height of measurement)
alpha = 1/7
# Implemented
mu = 1.81*10.**-5.  #Pa-s
rho = 1.184         #kg/ms
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
x_p = np.arange(0,(Lx+dx)+dx + dx, dx)   # Adding extra column of ghost nodes at right side boundary
XP, YP = np.meshgrid(x_p, y_p)      # Recreating meshgrid for plotting
p = np.zeros([len(y_p),len(x_p)])

u, v, p = apply_boundary_conditions(u,v,p)

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

u_frac = u*0
v_frac = v*0

interior_x_u = range(1, len(x_u) - 1)
interior_y_u = range(1, len(y_u) - 1)

interior_x_v = range(1, len(x_v) - 1)
interior_y_v = range(1, len(y_v) - 1)
b = np.zeros_like(p)


for n in range(nt):
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

    # For-loop version for b calculation (pressure Poisson)
    for j in range(1, len(y_p) - 1):
        for i in range(1, len(x_p) - 1):
            b[j, i] = (1 / dt
                       * ((u_frac[j,i] - u_frac[j,i-1]) / dx + (v_frac[j,i] - v_frac[j-1,i]) / dy))

        
    tolerance = 1e-6  # Convergence tolerance
    max_iterations = 100000  # Maximum number of iterations to prevent infinite looping

    for it in range(max_iterations):  # Maximum number of Poisson iterations per time step
        pn = p.copy()
        p[1:-1,1:-1] = (((pn[1:-1,2:] + pn[1:-1,0:-2]) * dy**2 +
                         (pn[2:,1:-1] + pn[0:-2,1:-1]) * dx**2) /
                        (2 * (dx**2 + dy**2)) -
                        dx**2 * dy**2 / (2 * (dx**2 + dy**2)) * b[1:-1,1:-1])

        # pressure boundary conditions
        p = apply_boundary_conditions(u_frac, v_frac, p)[2]
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
    
    u_n, v_n, p_n = apply_boundary_conditions(u_n, v_n, p_n)
    
    
    vel_tolerance = 1e-6
    p_tolerance = 1e-4
    
    u_max_diff = np.max(np.abs(u_n - u))
    v_max_diff = np.max(np.abs(v_n - v))
    p_max_diff = np.max(np.abs(p_n - p_old))
    
    max_vel_diff = max(u_max_diff, v_max_diff)  # combined velocity criterion
    
    if max_vel_diff < vel_tolerance and p_max_diff < p_tolerance and n!= 0:
        print(f"Steady state, n: {n}")
        break
    
    u = u_n.copy()
    v = v_n.copy()
    p = p_n.copy()
    
    
    if n % 10 == 0:
        plot_velocity(u,v)
        print(f'n:{n}')
        print(f'udiff: {u_max_diff}')
        print(f'vdiff: {v_max_diff}')
        print(f'pdiff: {p_max_diff}')

    
    

""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
"""Plotting"""

# wind profile
y_padded = np.linspace(y[0] - dy / 2, y[-1] + dy / 2, len(u[:, 0]))
u_z = u_r * (y_padded/z_r)**alpha
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

# print(u[:])