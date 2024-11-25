#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 19 00:09:06 2024

@author: tejjolly
"""

import numpy as np
import matplotlib.pyplot as plt

"""Grid from class"""""
Lx = 4  # Length of domain in x-direction (change to 1)
Ly = 4  # Length of domain in y-direction (change to 1)

dx = 1  #(change to 0.25)
dy = .5  #(change to 0.25)

x = np.arange(0, Lx+dx, dx) # x-grid
y = np.arange(0, Ly+dy, dy) # y-grid

M = len(x)  # num grid points, x
N = len(y)  # num grid points, y

# Or below????
# u-vel
x_u = np.arange(0, (Lx+dx), dx)          # u-velocity in x-nodes
y_u = np.arange(0, (Ly+dy)+dy, dy)       # u-velocity in y-nodes
XU,YU = np.meshgrid(x_u,y_u)
u = np.zeros([len(x_u),len(y_u)])

# v-vel
x_v = np.arange(0, (Lx+dx)+dx, dx)       # v-velocity, x-nodes
y_v = np.arange(0,(Ly+dy), dy)           # v-velocity, y-nodes
XV,YV = np.meshgrid(x_v,y_v)
v = np.zeros([len(x_v),len(y_v)])


# pressure
x_p = np.arange(0,(Lx+dx)+dx, dx)        # Pressure x-nodes
y_p = np.arange(0,(Ly+dy)+dy, dy)        # Pressure y-nodes
XP, YP = np.meshgrid(x_p, y_p)
p = np.zeros([len(x_p),len(y_p)])

""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
"""wind"""
u_r = 8.49      # m/s (reference wind speed)
z_r = 10        # meters (reference height of measurement)
alpha = 1/7
u_z = -u_r * (y/z_r)**alpha

""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
"""BCs"""
# Outflow BC at left boundary
x_u = np.arange(0, (Lx+dx) + dx, dx)  # Adding extra column of ghost nodes at right side boundary
XU,YU = np.meshgrid(x_u,y_u)        # Recreating meshgrid for plotting
u = np.zeros([len(x_u),len(y_u)])   # Recreating array for u velocity
u[-1,:] = u[-2,:]                   # u_M+1,1 = u_M,1

x_p = np.arange(0,(Lx+dx)+dx + dx, dx)   # Adding extra column of ghost nodes at right side boundary
XP, YP = np.meshgrid(x_p, y_p)      # Recreating meshgrid for plotting
p = np.zeros([len(x_p),len(y_p)])   # Recreating array for Pressure
p[-1,:] = -p[-2,:]                  # P_M+1,1 = -P_M,1


# Inflow BC at left boundary
u[0,:] = .01    # velocity at left-side boundary = 0.1 m/s
v[0,:] = 0      # v = 0 at ghost nodes left of left-side boundary
p[0,:] = p[1,:] #dP/dx|x=0 = 0


# Top and Bottom boundaries
u[:,0] = -u[:,1]    # no slip at bottom wall
u[:,-1] = u[:,-2]   # free shear at top boundary

# TODO: BCs for v #############################################################

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
plt.legend()


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
plt.ylim([0-1,Ly+1])
plt.legend()
plt.show()
