#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
@author: Dr. Gabor Zavodszky
"""

from PIL import Image
import cv2
import numpy as np
import matplotlib.pyplot as plt

###### Flow definition #########################################################
maxIter = 500  # Total number of time iterations.
Re = 50  # Reynolds number.
# nx, ny = 520, 180
nx, ny = 520, 180
ly = ny-1.0
q = 9  # Lattice dimensions and populations.
u_LB = 0.04                       # Velocity in lattice units.
nu_lb = u_LB*1/Re
Omega = 1.0 / (3.*nu_lb+0.5)  # Relaxation parameter.

###### Lattice Constants #######################################################
# Lattice velocities.
c = np.array([(x, y) for x in [0, -1, 1] for y in [0, -1, 1]])
w = 1./36. * np.ones(q)                                   # Lattice weights.
w[np.asarray([np.linalg.norm(ci) < 1.1 for ci in c])] = 1./9.
w[0] = 4./9.
noslip = [c.tolist().index((-c[i]).tolist()) for i in range(q)]
# Unknown on right wall.
i1 = np.arange(q)[np.asarray([ci[0] < 0 for ci in c])]
i2 = np.arange(q)[np.asarray([ci[0] == 0 for ci in c])]  # Vertical middle.
i3 = np.arange(q)[np.asarray([ci[0] > 0 for ci in c])]  # Unknown on left wall.
###### Function Definitions ####################################################


def sumpop(fin):
    """
    Helper function for density computation.
    """
    return np.sum(fin, axis=0)


def equilibrium(rho, u):
    """
    Equilibrium distribution function.
    """
    cu = 3.0 * np.dot(c, u.transpose(1, 0, 2))
    usqr = 3./2.*(u[0]**2+u[1]**2)
    feq = np.zeros((q, nx, ny))
    for i in np.arange(q):
        feq[i, :, :] = rho*w[i]*(1+cu[i]+0.5*cu[i]**2-usqr)
    return feq


######## Setup: Create Geometry ##########
img = Image.open('geometry.png')
img = img.convert('1')
geometry = 1 - np.asarray(img, dtype=bool)
##### Setup: velocity inlet with perturbation ########
vel = np.fromfunction(lambda d, x, y: (1-d)*u_LB *
                      (1.0+1e-4*np.sin(y/ly*2*np.pi)), (2, nx, ny))
feq = equilibrium(1, vel)
fin = feq.copy()

##### Main time loop ##########################################################
for t in np.arange(maxIter):
    fin[i1, -1, :] = fin[i1, -2, :]  # Right wall: outflow condition.'
    # Calculate macroscopic density and velocity.
    rho = sumpop(fin)
    u = np.dot(c.transpose(), fin.transpose((1, 0, 2)))/rho

    #########################################################################
    # Left wall: compute density from known populations.
    u[:, 0, :] = vel[:, 0, :]
    rho[0, :] = 1/(1-u[0, 0, :]) * \
        (sumpop(fin[i2, 0, :])+2.*sumpop(fin[i1, 0, :]))
    feq = equilibrium(rho, u)  # Left wall: Zou/He boundary condition.
    fin[i3, 0, :] = fin[i1, 0, :] + feq[i3, 0, :] - fin[i1, 0, :]
    #########################################################################

    # Collision
    fout = fin - Omega * (fin - feq)
    for i in np.arange(q):
        fout[i, geometry] = fin[noslip[i], geometry]

    # Streaming
    for i in np.arange(q):
        fin[i, :, :] = np.roll(
            np.roll(fout[i, :, :], c[i, 0], axis=0), c[i, 1], axis=1)
    if (t % 10 == 0):  # Visualization
        print("Iteration:", t)
        image = np.sqrt(u[0]**2+u[1]**2).transpose()
        plt.imshow(image, cmap=plt.get_cmap('Reds'))
        plt.pause(0.0001)
