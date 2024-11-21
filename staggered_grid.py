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
dy = 1  #(change to 0.25)

x = np.arange(0, Lx+dx, dx) # x-grid
y = np.arange(0, Ly+dy, dy) # y-grid

M = len(x)  # num grid points, x
N = len(y)  # num grid points, y

# Or below????
# u-vel
u_x = np.arange(0, (Lx+dx), dx)          # u-velocity in x-nodes
u_y = np.arange(0, (Ly+dy)+dy, dy)       # u-velocity in y-nodes

# v-vel
v_x = np.arange(0, (Lx+dx)+dx, dx)       # v-velocity, x-nodes
v_y = np.arange(0,(Ly+dy), dy)           # v-velocity, y-nodes

# pressure
p_x = np.arange(0,(Lx+dx)+dx, dx)        # Pressure x-nodes
p_y = np.arange(0,(Ly+dy)+dy, dy)        # Pressure y-nodes

print("\nClass grid")
print('u_x:', len(u_x))
print('u_y:', len(u_y))
print('v_x:', len(v_x))
print('v_y:', len(v_y))
print('p_x:', len(p_x))
print('p_y:', len(p_y))




""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
"""Grid for project"""""
height_bldg = 65    # ft
width_bldg = 65     # ft
width_domain = (2*width_bldg) + (4 * width_bldg)    # (upstream length) + (downstream length)
height_domain = (3*height_bldg)                     # length

ft2m = 0.3048 # Conversion, feet to meters

Lx = width_domain * ft2m   # Length of domain in x-direction (change to 1)
Ly = height_domain * ft2m    # Length of domain in y-direction (change to 1)

dy = 0.03
dx = dy*5

x = np.arange(0, Lx+dx, dx) # x-grid
y = np.arange(0, Ly+dy, dy) # y-grid

M = len(x)  # num grid points, x
N = len(y)  # num grid points, y

# Or below????
# u-vel
u_x = np.arange(0, (Lx+dx), dx)          # u-velocity in x-nodes
u_y = np.arange(0, (Ly+dy)+dy, dy)       # u-velocity in y-nodes

# v-vel
v_x = np.arange(0, (Lx+dx)+dx, dx)       # v-velocity, x-nodes
v_y = np.arange(0,(Ly+dy), dy)           # v-velocity, y-nodes

# pressure
p_x = np.arange(0,(Lx+dx)+dx, dx)        # Pressure x-nodes
p_y = np.arange(0,(Ly+dy)+dy, dy)        # Pressure y-nodes

print('\nProject grid')
print('u_x:', len(u_x))
print('u_y:', len(u_y))
print('v_x:', len(v_x))
print('v_y:', len(v_y))
print('p_x:', len(p_x))
print('p_y:', len(p_y))


""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
"""wind"""
u_r = 8.49      # m/s (reference wind speed)
z_r = 10        # meters (reference height of measurement)
alpha = 1/7
u_z = -u_r * (y/z_r)**alpha







""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
"""Plotting"""
plt.figure()
plt.plot(u_z,y)
plt.title("Wind profile")
plt.xlabel("Velocity [m/s]")
plt.ylabel("height [m]")
plt.gca().yaxis.set_label_position("right")
plt.gca().yaxis.tick_right()
plt.axhline(y=height_bldg * ft2m, color='gray', linestyle='--', linewidth=1, label=f'bldg. height')  # Add horizontal line
# plt.axis('equal')
plt.gca().set_aspect(.25, adjustable='box')
plt.xlim(right=0)  
plt.legend()


X, Y = np.meshgrid(x, y)
plt.figure(figsize=(12, 8))
plt.plot(X[:, ::100], Y[:, ::100], 'k-', linewidth=0.5)
plt.plot(X.T[:, ::100], Y.T[:, ::100], 'k-', linewidth=0.5)
plt.scatter(X[::100, ::100], Y[::100, ::100], color='red', s=10,label='100th Grid points')

plt.plot([0, Lx, Lx, 0, 0], [0, 0, Ly, Ly, 0], color='blue', linestyle='--', linewidth=2, label="Domain Boundary")

plt.xlabel('x (m)')
plt.ylabel('y (m)')
plt.title('Improved Grid Visualization')
plt.axis('equal')  # Ensure equal aspect ratio
plt.grid(False)  # Turn off default matplotlib grid
plt.xlim([0-10, Lx+10])
plt.ylim([0-10,Ly+10])
plt.legend()
plt.show()
