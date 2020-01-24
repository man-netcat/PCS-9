from __future__ import division, print_function

import sys

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


def update(frame):
    # Progress our simulation
    for _ in range(1):  # adjust number of steps for smooth animation
        fin[i1, -1, :] = fin[i1, -2, :]  # Right wall: outflow condition.'
        # Calculate macroscopic density and velocity.
        rho = sumpop(fin)
        u = np.dot(c.transpose(), fin.transpose((1, 0, 2)))/rho

        # Left wall: compute density from known populations.
        u[:, 0, :] = vel[:, 0, :]
        rho[0, :] = 1/(1-u[0, 0, :]) * \
            (sumpop(fin[i2, 0, :])+2.*sumpop(fin[i1, 0, :]))
        feq = equilibrium(rho, u)  # Left wall: Zou/He boundary condition.
        fin[i3, 0, :] = fin[i1, 0, :] + feq[i3, 0, :] - fin[i1, 0, :]

        # Collision
        fout = fin - Omega * (fin - feq)
        for i in np.arange(9):
            fout[i, geometry.transpose()] = fin[noslip[i],
                                                geometry.transpose()]

        # Streaming
        for i in np.arange(9):
            fin[i, :, :] = np.roll(
                np.roll(fout[i, :, :], c[i, 0], axis=0), c[i, 1], axis=1)
    fluidImage.set_array(
        # Velocity
        # vis[0](u[0], u[1]).transpose()
        # Density
        rho.transpose()
    )
    printProgressBar(frame + 1, frames, prefix='Progress:',
                     suffix='Complete', length=50)
    return (fluidImage, wallImage)  # return the figure elements to redraw


def printProgressBar(iteration, total, prefix='', suffix='', decimals=1,
                     length=100, fill='█', printEnd="\r"):
    percent = ("{0:." + str(decimals) + "f}"
               ).format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print('\r%s |%s| %s%% %s' % (prefix, bar, percent, suffix), end=printEnd)
    # Print New Line on Complete
    if iteration == total:
        print()


if __name__ == '__main__':
    geometry = np.asarray(Image.open(sys.argv[1]).convert('1'))
    height, width = geometry.shape

    Re = 50  # Reynolds number.
    u_0 = 0.1  # Velocity in lattice units.
    viscosity = u_0/Re
    Omega = 1.0 / (3.*viscosity+0.5)  # Relaxation parameter.

    frames = 400
    fps = 60

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
    # vel = np.fromfunction(lambda d, x, y: (1-d)*u_0 *
    #                       (1.0+1e-4*np.sin(y/(height-1)*2*np.pi)),
    #                       (2, width, height))
    vel = np.array([np.full((width, height), u_0),
                    np.full((width, height), 0)])
    feq = equilibrium(1, vel)
    fin = feq.copy()
    rho = sumpop(fin)
    u = np.dot(c.transpose(), fin.transpose((1, 0, 2)))/rho

    # Graphics helper stuff here
    fig = plt.figure(figsize=(8, 3))
    vis = [mag, 0.2]
    # Velocity
    fluidImage = plt.imshow(
        vis[0](vel[0], vel[1]).transpose(),
        origin='lower',
        norm=plt.Normalize(-vis[1], vis[1]),
        interpolation='none',
        cmap=plt.get_cmap('Reds')
    )
    # Density
    # fluidImage = plt.imshow(
    #     rho.transpose(),
    #     origin='lower',
    #     norm=plt.Normalize(1, 1.1),
    #     interpolation='none',
    #     cmap=plt.get_cmap('jet')
    # )
    wImageArray = np.zeros((height, width, 4), np.uint8)  # an RGBA image
    wImageArray[geometry, 3] = 255
    wallImage = plt.imshow(wImageArray, origin='lower', interpolation='none')

    animate = matplotlib.animation.FuncAnimation(
        fig, update, interval=1, blit=True, frames=frames)
    # plt.show()
    Writer = matplotlib.animation.writers['ffmpeg']
    writer = Writer(fps=fps, metadata=dict(artist='Me'), bitrate=1800)
    animate.save('im.mp4', writer=writer)
