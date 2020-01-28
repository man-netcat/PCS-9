from __future__ import division, print_function

import argparse

import matplotlib.animation
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image

# Add Command-Line Arguments
parser = argparse.ArgumentParser(description='Simulate Blood Flow')
parser.add_argument('geometry', help='Path to geometry for simulation')
parser.add_argument('output', help='Output Method: plot or video')
parser.add_argument('-o', '--out', default='out.mp4', help='Video Name')
parser.add_argument('-f', '--fps', default=30, help='Frames Per Second')
parser.add_argument('-l', '--length', default=30,
                    help='Video Length in seconds')
parser.add_argument('-R', '--Reynolds', default=10, help='Reynolds Number')
parser.add_argument('-U', '--velocity', default=0.1, help='Initial Velocity')
parser.add_argument('--method', default='velocity',
                    help='Display Method: velocity or density')
parser.add_argument('-g', '--graph', default='graph.png',
                    help='Graph File Name')
args = parser.parse_args()


def sum_populations(fin):
    '''Helper function for adding up distributions'''
    return np.sum(fin, axis=0)


def equilibrium(rho, u):
    '''
    Calculates the equilibrium distribution from a given density and
    velocity
    '''
    cu = 3.0 * np.dot(c, u.transpose(1, 0, 2))
    usqr = 3./2.*(u[0]**2+u[1]**2)
    feq = np.zeros((9, width, height))
    for i in range(9):
        feq[i, :, :] = rho*w[i]*(1.+cu[i]+0.5*cu[i]**2-usqr)
    return feq


def calculate_magnitude(u_x, u_y):
    '''
    Helper function that calculates the magnitude of individual x and y
    components of the velocity
    '''
    return np.sqrt(u_x**2+u_y**2)


# Initialise Geometry and its dimensions
geometry = np.asarray(Image.open(args.geometry).convert('1'))
height, width = geometry.shape

# Initialise simulation variables
viscosity = float(args.velocity)/int(args.Reynolds)
# Relaxation parameter.
Omega = 1.0 / (3.*viscosity+0.5)

# Initialise amount of frames for video
frames = int(args.fps)*int(args.length)

# Lattice Constants
c = np.array([
    [0,  0],
    [1,  0],
    [0,  1],
    [-1,  0],
    [0, -1],
    [1,  1],
    [-1,  1],
    [-1, -1],
    [1, -1]
])
w = np.array([4 / 9, 1 / 9, 1 / 9, 1 / 9, 1 /
              9, 1 / 36, 1 / 36, 1 / 36, 1 / 36])
noslip = np.array([0, 3, 4, 1, 2, 7, 8, 5, 6])
x_neg = np.array([3, 6, 7])
x_neu = np.array([0, 2, 4])
x_pos = np.array([1, 5, 8])
y_neg = np.array([4, 7, 8])
y_neu = np.array([0, 1, 3])
y_pos = np.array([2, 5, 6])

# Calculate macroscopic density and velocity.
vel = np.array([np.full((width, height), args.velocity),
                np.full((width, height), 0)])
feq = equilibrium(1, vel)
fin = feq.copy()
rho = sum_populations(fin)
u = np.dot(c.transpose(), fin.transpose((1, 0, 2)))/rho


def update(frame):
    '''
    Main step function. This function is called every iteration and updates
    the density and velocity in the simulation.
    '''
    # Right wall: outflow condition.
    fin[x_neg, -1, :] = fin[x_neg, -2, :]

    # Calculate macroscopic density and velocity.
    rho = sum_populations(fin)
    u = np.dot(c.transpose(), fin.transpose((1, 0, 2)))/rho

    # Left wall: compute density from known populations.
    u[:, 0, :] = vel[:, 0, :]
    rho[0, :] = 1/(1-u[0, 0, :]) * \
        (sum_populations(fin[x_neu, 0, :]) +
         2.*sum_populations(fin[x_neg, 0, :]))
    feq = equilibrium(rho, u)

    # Left wall: Zou/He boundary condition.
    fin[x_pos, 0, :] = feq[x_pos, 0, :]

    # Collision
    fout = fin - Omega * (fin - feq)
    for i in np.arange(9):
        fout[i, geometry.transpose()] = fin[noslip[i], geometry.transpose()]

    # Streaming
    for i in np.arange(9):
        fin[i, :, :] = np.roll(
            np.roll(fout[i, :, :], c[i, 0], axis=0), c[i, 1], axis=1)

    # Update Data
    for desc in keypoints.keys():
        if args.method == 'velocity':
            data[desc] = np.append(data[desc], vis[0](
                u[0], u[1]).transpose()[keypoints[desc]])
        elif args.method == 'density':
            data[desc] = np.append(
                data[desc], rho.transpose()[keypoints[desc]])

    # Update Array
    if args.method == 'velocity':
        fluidImage.set_array(vis[0](u[0], u[1]).transpose())
    elif args.method == 'density':
        fluidImage.set_array(rho.transpose())

    if args.output == 'video':
        printProgressBar(frame + 1, frames, prefix='Progress:',
                         suffix='Complete', length=50)

    return (fluidImage, wallImage)


def printProgressBar(iteration, total, prefix='', suffix='', decimals=1,
                     length=100, fill='█', printEnd="\r"):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent
                                  complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
        printEnd    - Optional  : end character (e.g. "\r", "\r\n") (Str)
    (Source:)
    https://stackoverflow.com/questions/3173320/text-progress-bar-in-the-console
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 *
                                                     (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print('\r%s |%s| %s%% %s' % (prefix, bar, percent, suffix), end=printEnd)
    # Print New Line on Complete
    if iteration == total:
        print()


if __name__ == '__main__':
    # Dictionary for data
    with open("./out/probe.txt", "r") as f:
        keypoints = {}
        data = {}
        for line in f:
            desc, x, y = line.split(',')
            keypoints[desc] = (round(float(y)), round(float(x)))
            data[desc] = np.array([])

    # Initialise Figure
    fig = plt.figure(figsize=(8, 3))
    vis = [calculate_magnitude, 0.2]
    if args.method == 'velocity':
        fluidImage = plt.imshow(
            vis[0](vel[0], vel[1]).transpose(),
            origin='lower',
            norm=plt.Normalize(-vis[1], vis[1]),
            interpolation='none',
            cmap=plt.get_cmap('jet')
        )
    elif args.method == 'density':
        fluidImage = plt.imshow(
            rho.transpose(),
            origin='lower',
            norm=plt.Normalize(1, 1.1),
            interpolation='none',
            cmap=plt.get_cmap('jet')
        )
    else:
        raise ValueError("Invalid Method Specified")
    wImageArray = np.zeros((height, width, 4), np.uint8)

    # Set the alpha value of the geometry to 255 (not translucent)
    wImageArray[geometry, 3] = 255
    wallImage = plt.imshow(wImageArray, origin='lower', interpolation='none')

    # Initialise Animation
    animate = matplotlib.animation.FuncAnimation(
        fig, update, interval=1, blit=True, frames=frames)
    if args.output == 'plot':
        plt.show()
    elif args.output == 'video':
        Writer = matplotlib.animation.writers['ffmpeg']
        writer = Writer(fps=int(args.fps), metadata=dict(
            artist='Me'), bitrate=1800)
        animate.save(args.out, writer=writer)
    else:
        raise ValueError("Invalid Output Method Specified")

    # Plot Graph
    fig2, ax = plt.subplots()

    for desc in keypoints.keys():
        plt.plot(np.arange(frames+1), data[desc], label=desc)
    ax.set_xlabel('timesteps')
    ax.set_ylabel(args.method)
    plt.legend()
    plt.savefig(args.graph)
