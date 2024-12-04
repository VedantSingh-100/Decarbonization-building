#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 19 00:09:06 2024

@author: tejjolly
"""

import numpy as np
import matplotlib.pyplot as plt

def apply_boundary_conditions(u, v, p)


"""Grid from class"""""
Lx = 1 * 20  # Length of domain in x-direction (change to 1)
Ly = 1  * 20 # Length of domain in y-direction (change to 1)

dx = .25 * 20  #(change to 0.25)
dy = .5 * 20   #(change to 0.25)

x = np.arange(0, Lx+dx, dx) # x-grid
y = np.arange(0, Ly+dy, dy) # y-grid

M = len(x)  # num grid points, x
N = len(y)  # num grid points, y

# Or below????
# u-vel
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
u_z = -u_r * (y/z_r)**alpha

mu = 1
rho = 1
dt = 1

""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
"""BCs"""
# Outflow BC at right boundary
x_u = np.arange(0, (Lx+dx) + dx, dx)  # Adding extra column of ghost nodes at right side boundary
XU,YU = np.meshgrid(x_u,y_u)        # Recreating meshgrid for plotting
u = np.zeros([len(y_u),len(x_u)])
u[:,-1] = u[:,-2]                   # u(M+1,1) = u(M,1)

x_p = np.arange(0,(Lx+dx)+dx + dx, dx)   # Adding extra column of ghost nodes at right side boundary
y_p = np.arange(0,(Ly+dy)+dy + dy, dy)   # Adding extra row of ghost nodes at top side boundary
XP, YP = np.meshgrid(x_p, y_p)      # Recreating meshgrid for plotting
p = np.zeros([len(y_p),len(x_p)])
p[:,-2] = 0                         # P|x=Lx = 0 (0 gauge outlet)
p[:,-1] = -p[:,-2]                  # P(M+1,1) = -P(M,1)


# Inflow BC at left boundary
u[:,0] = .01    # velocity at left-side boundary = 0.1 m/s
p[:,0] = p[:,1] # dP/dx|x=0 = 0


# Top and Bottom boundaries
u[0,:] = -u[1,:]    # no slip at bottom wall
u[-1,:] = u[-2,:]   # free shear at top boundary
p[0,:] = p[1,:]     # bottom wall, horizontal wall dP/dy = 0 

p[:,-2] = 0         # P|y=Ly = 0 (0 gauge outlet)
p[-1,:] = -p[-2,:]  # P(N+1,1) = -P(N,1)

# TODO: BCs for v #############################################################
y_v = np.arange(0,(Ly+dy) + dy, dy) # Adding extra row of ghost nodes at top side boundary
XV,YV = np.meshgrid(x_v,y_v)
v = np.zeros([len(y_v),len(x_v)])

v[:,-1] = v[:,-2]   # Outflow at right boundary
v[:,0] = 0          # left side boundary, v = 0 at ghost nodes left of left-side boundary
v[0,:] = 0          # bottom wall (ground), v = 0
v[-1,:] = v[-2,:]   # top wall (free shear), dv/dy = 0


""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
"""Internal Nodes"""

"""Fractional step method, version 1"""

"""Step 1, ignore P in the momentum eqn."""
u_new = u*0
v_new = v*0
p_new = p*0

# x-momentum eqn. components
interior_x_u = range(1,len(x_u) - 1)
interior_y_u = range(1,len(y_u) - 1)

interior_x_v = range(1,len(x_u) - 1)
interior_y_v = range(1,len(y_u) - 1)

for i in range(len(x_u)):
    for j in range(len(y_u)):
        
        duudx = 1/dx * (((u[j,i]+u[j,i+1])/2)**2 - ((u[j,i]+u[j,i-1])/2)**2)
        
        dvudy = 1/dy * (v[j,i] + v[j,i+1])/2 * (u[j,i] + u[j+1,i])/2 \
                - (v[j-1,i] + v[j-1,i+1])/2 * (u[j-1,i] + u[j,i])/2
                
        d2udx2 = (u[j,i+1] - 2*u[j,i] + u[j,i-1])/dx**2
        
        d2udy2 = u[j+1,i] - 2*u[j,i] + u[j-1,i]/dy**2
        
        u_new[j,i] = dt * (mu/rho * (d2udx2 + d2udy2) - (duudx + dvudy)) + u[j,i]
    
# y-momentum eqn. components
for i in range(len(x_v)):
    for j in range(len(y_v)):
        
        duvdx = 1/dx * (u[j,i] + u[j+1,i])/2 * (v[j,i+1] + v[j,i])/2 \
                - (u[j,i-1] + u[j+1,i-1])/2 * (v[j,i] + v[j,i-1])/2
                
        dvvdy = 1/dy * ((v[j,i] + v[j+1,i])/2)**2 - ((v[j,i] + v[j-1,i])/2)**2
        
        d2vdx2 = (v[j,i+1] - 2*v[j,i] + v[j,i-1])/dx**2
        
        d2vdy2 = (v[j+1,i] - 2*v[j,i] + v[j-1,i])/dy**2
        
        v_new[j,i] = dt * (mu/rho * (d2vdx2 + d2vdy2) - (duvdx + dvvdy)) + v[j,i]
        
u[j,i] = u_new[j,i]
v[j,i] = v_new[j,i]


"""Step 2, solve Pressure Poisson eq."""
interior_x_p = range(1,len(x_p) - 1)
interior_y_p = range(1,len(y_p) - 1)

for i in interior_x_p:
    for j in interior_y_p:
            b = 1/dt * ((u[i,j] - u[i-1,j]) / dx + \
                          (v[i,j] - v[i,j-1]) / dy)
                
    for it in range(50):  # Number of Poisson iterations per time step
        pn = p.copy()
        p[1:-1,1:-1] = (((pn[1:-1,2:] + pn[1:-1,0:-2]) * dy**2 +
                         (pn[2:,1:-1] + pn[0:-2,1:-1]) * dx**2) /
                        (2 * (dx**2 + dy**2)) -
                        dx**2 * dy**2 / (2 * (dx**2 + dy**2)) * b[1:-1,1:-1])

        # Apply pressure boundary conditions
        p = apply_boundary_conditions(u, v, p)[2]


### Probably need to combine the above into one loop?^^^^

""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
"""Plotting"""
plt.figure()
plt.plot(u_z,y)
plt.title("Wind profile")
plt.xlabel("Velocity [m/s]")
plt.ylabel("height [m]")
plt.gca().yaxis.set_label_position("right")
plt.gca().yaxis.tick_right()
# plt.axhline(y=height_bldg * ft2m, color='gray', linestyle='--', linewidth=1, label=f'bldg. height')  # Add horizontal line
plt.axis('equal')
plt.gca().set_aspect(.25, adjustable='box')
plt.xlim(right=0)  
# plt.legend()


plt.figure(figsize=(12, 8))

XP = XP-dx/2
YP = YP-dy/2
plt.scatter(XP, YP, color='green', marker = '+', s=45,label='p Grid points')

YU = YU-dy/2
plt.scatter(XU, YU, color='blue', marker='>', s=45,label='u - Grid points')

XV = XV-dx/2
plt.scatter(XV, YV, color='red', marker='^', s=45,label='v - Grid points')


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

print('u_shape: ',u.shape)
print('v_shape: ',v.shape)
print('P_shape: ',p.shape)


