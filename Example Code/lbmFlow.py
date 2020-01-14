#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
@author: Dr. Gabor Zavodszky
"""
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image

###### Flow definition #########################################################
# Total number of time iterations.
maxIter = 200000
# Reynolds number.
Re = 50
# Lattice dimensions and populations.
nx = 520
ny = 180
ly = ny-1.0
q = 9
# Coordinates of the cylinder.
cx = nx/4
cy = ny/2
r = ny/9
# Velocity in lattice units.
uLB = 0.04
nulb = uLB*r/Re
# Relaxation parameter.
omega = 1.0 / (3.*nulb+0.5)

###### Lattice Constants #######################################################
# Lattice velocities.
c = np.array([(x, y) for x in [0, -1, 1] for y in [0, -1, 1]])
t = 1./36. * np.ones(q)                                   # Lattice weights.
t[np.asarray([np.linalg.norm(ci) < 1.1 for ci in c])] = 1./9.
t[0] = 4./9.
print(t)
exit()
noslip = [c.tolist().index((-c[i]).tolist()) for i in range(q)]
# Unknown on right wall.
i1 = np.arange(q)[np.asarray([ci[0] < 0 for ci in c])]
# Vertical middle.
i2 = np.arange(q)[np.asarray([ci[0] == 0 for ci in c])]
# Unknown on left wall.
i3 = np.arange(q)[np.asarray([ci[0] > 0 for ci in c])]

###### Function Definitions ####################################################


# Helper function for density computation.
def sumpop(fin): return np.sum(fin, axis=0)


def equilibrium(rho, u):              # Equilibrium distribution function.
    cu = 3.0 * np.dot(c, u.transpose(1, 0, 2))
    usqr = 3./2.*(u[0]**2+u[1]**2)
    feq = np.zeros((q, nx, ny))
    for i in range(q):
        feq[i, :, :] = rho*t[i]*(1.+cu[i]+0.5*cu[i]**2-usqr)
    return feq


###### Setup: cylindrical obstacle and velocity inlet with perturbation ########
# obstacle = np.fromfunction(lambda x, y: (x-cx)**2+(y-cy)**2 < r**2, (nx, ny))
img = Image.open('geometry.png')
img = img.convert('1')
obstacle = (1 - np.asarray(img, dtype=bool)).transpose()
# print(obstacle.shape)
# print(obstacle)
# plt.imshow(obstacle)
# plt.show()
vel = np.fromfunction(lambda d, x, y: (1-d)*uLB *
                      (1.0+1e-4*np.sin(y/ly*2*np.pi)), (2, nx, ny))
feq = equilibrium(1.0, vel)
fin = feq.copy()

plt.ion()

###### Main time loop ##########################################################
for time in range(maxIter):
    # Right wall: outflow condition.
    fin[i1, -1, :] = fin[i1, -2, :]
    # Calculate macroscopic density and velocity.
    rho = sumpop(fin)
    u = np.dot(c.transpose(), fin.transpose((1, 0, 2)))/rho

    # Left wall: compute density from known populations.
    u[:, 0, :] = vel[:, 0, :]
    rho[0, :] = 1./(1.-u[0, 0, :]) * \
        (sumpop(fin[i2, 0, :])+2.*sumpop(fin[i1, 0, :]))
    # Left wall: Zou/He boundary condition.
    feq = equilibrium(rho, u)
    fin[i3, 0, :] = fin[i1, 0, :] + feq[i3, 0, :] - fin[i1, 0, :]
    # Collision step.
    fout = fin - omega * (fin - feq)
    for i in range(q):
        fout[i, obstacle] = fin[noslip[i], obstacle]
    # Streaming step.
    for i in range(q):
        fin[i, :, :] = np.roll(
            np.roll(fout[i, :, :], c[i, 0], axis=0), c[i, 1], axis=1)
    # Visualization
    if (time % 1 == 0):
        print("Iteration:", time)
        plt.imshow(np.sqrt(u[0]**2+u[1]**2).transpose(), cmap="Reds")
        plt.pause(0.01)
