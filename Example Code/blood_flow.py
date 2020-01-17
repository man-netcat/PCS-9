from __future__ import division, print_function

import sys
import time

import matplotlib.animation
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image


def sumpop(fin): return np.sum(fin, axis=0)


def equilibrium(rho, u):
    cu = 3.0 * np.dot(c, u.transpose(1, 0, 2))
    usqr = 3./2.*(u[0]**2+u[1]**2)
    feq = np.zeros((9, width, height))
    for i in range(9):
        feq[i, :, :] = rho*t[i]*(1.+cu[i]+0.5*cu[i]**2-usqr)
    return feq


def mag(u_x, u_y):
    return np.sqrt(u_x**2+u_y**2)


def stream():
    global c, fin, fout, feq
    for i in range(9):
        fin[i, :, :] = np.roll(
            np.roll(fout[i, :, :], c[i, 0], axis=0), c[i, 1], axis=1)


# Collide particles within each cell to redistribute velocities
# (can be optimized a little more):
def collide():
    global c, fin, fout, feq
    fout = fin - omega * (fin - feq)
    for i in range(9):
        fout[i, geometry.transpose()] = fin[noslip[i], geometry.transpose()]


# Impose boundary conditions
def boundary():
    global c, fin, fout, feq
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


def nextFrame(arg):
    global startTime

    fin[i1, -1, :] = fin[i1, -2, :]

    rho = sumpop(fin)
    u = np.dot(c.transpose(), fin.transpose((1, 0, 2)))/rho

    # Progress our simulation
    for step in range(1):  # adjust number of steps for smooth animation
        boundary()
        collide()
        stream()

    fluidImage.set_array(vis[0](u[0], u[1]).transpose())
    return (fluidImage, wallImage)  # return the figure elements to redraw


if __name__ == '__main__':
    geometry = np.asarray(Image.open(sys.argv[1]).convert('1')) == 0
    height, width = geometry.shape

    viscosity = 0.02                    # viscosity
    omega = 1 / (3*viscosity + 0.5)     # parameter for "relaxation"
    u0 = 0.1                            # initial and in-flow speed

    # Lattice Constants
    c = np.array([(x, y) for x in [0, -1, 1] for y in [0, -1, 1]])
    t = 1./36. * np.ones(9)
    t[np.asarray([np.linalg.norm(ci) < 1.1 for ci in c])] = 1./9.
    t[0] = 4./9.
    noslip = [c.tolist().index((-c[i]).tolist()) for i in range(9)]
    i1 = np.arange(9)[np.asarray([ci[0] < 0 for ci in c])]
    i2 = np.arange(9)[np.asarray([ci[0] == 0 for ci in c])]
    i3 = np.arange(9)[np.asarray([ci[0] > 0 for ci in c])]

    # Calculate macroscopic density and velocity.
    vel = np.array([np.full((width, height), u0), np.full((width, height), 0)])
    feq = equilibrium(1.0, vel)
    fin = feq.copy()
    rho = sumpop(fin)

    # Graphics helper stuff here
    theFig = plt.figure(figsize=(8, 3))
    vis = [mag, 0.2]
    # vis = [curl, 0.02]
    fluidImage = plt.imshow(
        vis[0](vel[0], vel[1]).transpose(),
        origin='lower',
        norm=plt.Normalize(-vis[1], vis[1]),
        interpolation='none'
    )
    wImageArray = np.zeros((height, width, 4), np.uint8)  # an RGBA image
    # set alpha=255 wall sites only
    print(geometry.shape)
    wImageArray[geometry, 3] = 255
    wallImage = plt.imshow(wImageArray, origin='lower', interpolation='none')

    startTime = time.clock()    # Start clock for performance counting

    animate = matplotlib.animation.FuncAnimation(
        theFig, nextFrame, interval=10, blit=True)
    plt.show()
