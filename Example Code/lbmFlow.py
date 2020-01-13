#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
@author: Dr. Gabor Zavodszky
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

###### Flow definition #########################################################
maxIter = 200000  # Total number of time iterations.
Re = 2000  # Reynolds number.
nx = 520
ny = 180
ly = ny-1.0
q = 9  # Lattice dimensions and populations.
cx = nx/4
cy = ny/2
r = ny/9          # Coordinates of the cylinder.
uLB = 0.04                       # Velocity in lattice units.
nulb = uLB*r/Re
omega = 1.0 / (3.*nulb+0.5)  # Relaxation parameter.

###### Lattice Constants #######################################################
# Lattice velocities.
c = np.array([(x, y) for x in [-1, 0, 1] for y in [-1, 0, 1]])
t = np.array([
    1/36, 1/9, 1/36,
    1/9, 4/9, 1/9,
    1/36, 1/9, 1/36
])
noslip = np.arange(9)[::-1]
i1 = np.arange(3)
i2 = np.arange(3, 6)
i3 = np.arange(6, 9)
# exit()
###### Function Definitions ####################################################


# Helper function for density computation.
def sumpop(fin): return np.sum(fin, axis=0)


def equilibrium(rho, u):              # Equilibrium distribution function.
    cu = 3.0 * np.dot(c, u.transpose(1, 0, 2))
    usqr = 3./2.*(u[0]**2+u[1]**2)
    feq = np.zeros((q, nx, ny))
    for i in np.arange(q):
        feq[i, :, :] = rho*t[i]*(1.+cu[i]+0.5*cu[i]**2-usqr)
    return feq


###### Setup: cylindrical obstacle and velocity inlet with perturbation ########
obstacle = np.fromfunction(lambda x, y: (x-cx)**2+(y-cy)**2 < r**2, (nx, ny))
vel = np.fromfunction(lambda d, x, y: (1-d)*uLB *
                      (1.0+1e-4*np.sin(y/ly*2*np.pi)), (2, nx, ny))
feq = equilibrium(1.0, vel)
fin = feq.copy()

plt.ion()

###### Main time loop ##########################################################
for time in np.arange(maxIter):
    fin[i1, -1, :] = fin[i1, -2, :]  # Right wall: outflow condition.
    rho = sumpop(fin)           # Calculate macroscopic density and velocity.
    u = np.dot(c.transpose(), fin.transpose((1, 0, 2)))/rho

    # Left wall: compute density from known populations.
    u[:, 0, :] = vel[:, 0, :]
    rho[0, :] = 1./(1.-u[0, 0, :]) * \
        (sumpop(fin[i2, 0, :])+2.*sumpop(fin[i1, 0, :]))

    feq = equilibrium(rho, u)  # Left wall: Zou/He boundary condition.
    fin[i3, 0, :] = fin[i1, 0, :] + feq[i3, 0, :] - fin[i1, 0, :]
    fout = fin - omega * (fin - feq)  # Collision step.
    for i in np.arange(q):
        fout[i, obstacle] = fin[noslip[i], obstacle]
    for i in np.arange(q):  # Streaming step.
        fin[i, :, :] = np.roll(
            np.roll(fout[i, :, :], c[i, 0], axis=0), c[i, 1], axis=1)

    if (time % 100 == 0):  # Visualization
        print("Iteration:", time)
        plt.imshow(np.sqrt(u[0]**2+u[1]**2).transpose(), cmap=cm.jet)
        plt.pause(0.0001)
