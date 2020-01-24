"""
Lattice-Boltzmann method for fluid simulation
D2Q9 grid
BGK collision operator
Karman-vortices
"""

from __future__ import division
from __future__ import print_function
import matplotlib.pyplot as plt
import matplotlib.animation
from PIL import Image
import numpy as np
import sys
import os
import cv2

# Fixed parameters:
width = 400
height = 160
viscosity = 0.05
u0 = 0.25                           # initial and in-flow speed
omega = 1 / (3*viscosity + 0.5)     # parameter for "relaxation"

# Custom parameters
geometry = "./out/geometry2.png"
videoname = "video.avi"
fps = 20
frames = 300

if len(sys.argv) > 1:
    geometry = sys.argv[1]
if len(sys.argv) > 2:
    videoname = sys.argv[2]
if len(sys.argv) > 3:
    fps = sys.argv[3]
if len(sys.argv) > 4:
    frames = sys.argv[4]

# Lattice weight constants for D2Q9
f_n = 4.0/9.0
o_n = 1.0/9.0
o_36 = 1.0/36.0

# Initialize arrays for a steady rightward flow
ux0 = u0
uy0 = 0.0
rho0 = np.ones((height, width))
ux02 = ux0 * ux0               # pre-compute terms used repeatedly...
uy02 = uy0 * uy0
u02 = ux02 + uy02
omu0215 = 1 - 1.5*u02         # "one minus u2 times 1.5"
uxuy0 = ux0 * uy0

# Set populations to f^eq(rho, ux0, uy0)
n0 = f_n * rho0 * omu0215
nN = o_n * rho0 * (omu0215 + 3*uy0 + 4.5*uy02)
nS = o_n * rho0 * (omu0215 - 3*uy0 + 4.5*uy02)
nE = o_n * rho0 * (omu0215 + 3*ux0 + 4.5*ux02)
nW = o_n * rho0 * (omu0215 - 3*ux0 + 4.5*ux02)
nNE = o_36 * rho0 * (omu0215 + 3*(ux0+uy0) + 4.5*(u02+2*uxuy0))
nNW = o_36 * rho0 * (omu0215 + 3*(-ux0+uy0) + 4.5*(u02-2*uxuy0))
nSE = o_36 * rho0 * (omu0215 + 3*(ux0-uy0) + 4.5*(u02-2*uxuy0))
nSW = o_36 * rho0 * (omu0215 + 3*(-ux0-uy0) + 4.5*(u02+2*uxuy0))

# Get macroscopic values from initial populations
rho = n0 + nN + nS + nE + nW + nNE + nSE + nNW + nSW
ux = (nE + nNE + nSE - nW - nNW - nSW) / rho
uy = (nN + nNE + nNW - nS - nSE - nSW) / rho

# Initialize wall locations:
wall = np.asarray(Image.open(geometry).convert('1')) == 1
geometry_input = [x for x in range(len(wall.transpose()[0])) if not wall.transpose()[0][x]]
geometry_outputs = [x for x in range(len(wall.transpose()[0])) if not wall.transpose()[-1][x]]

index = 0
for value in geometry_outputs:
    if index == 0:
        oldValue = value
        index += 1
        continue
    if oldValue - value < -1:
        break
    oldValue = value
    index += 1

geometry_output1 = geometry_outputs[:index]
geometry_output2 = geometry_outputs[index:]


n0_input_eq = n0.transpose()[0][geometry_input]
nN_input_eq = nN.transpose()[0][geometry_input]
nS_input_eq = nS.transpose()[0][geometry_input]
nE_input_eq = nE.transpose()[0][geometry_input]
nW_input_eq = nW.transpose()[0][geometry_input]
nNE_input_eq = nNE.transpose()[0][geometry_input]
nNW_input_eq = nNW.transpose()[0][geometry_input]
nSE_input_eq = nSE.transpose()[0][geometry_input]
nSW_input_eq = nSW.transpose()[0][geometry_input]

# Set up indices for fast evaluation of wall neighbors
barrierN = np.roll(wall,  1, axis=0)
barrierS = np.roll(wall, -1, axis=0)
barrierE = np.roll(wall,  1, axis=1)
barrierW = np.roll(wall, -1, axis=1)
barrierNE = np.roll(barrierN,  1, axis=1)
barrierNW = np.roll(barrierN, -1, axis=1)
barrierSE = np.roll(barrierS,  1, axis=1)
barrierSW = np.roll(barrierS, -1, axis=1)

# Move all particles by one step along their directions of motion (periodic boundary):
def stream():
    global rho, ux, uy, n0, nN, nS, nE, nW, nNE, nNW, nSE, nSW
    # axis 0 is north-south; + direction is north
    nN = np.roll(nN,   1, axis=0)
    nNE = np.roll(nNE,  1, axis=0)
    nNW = np.roll(nNW,  1, axis=0)
    nS = np.roll(nS,  -1, axis=0)
    nSE = np.roll(nSE, -1, axis=0)
    nSW = np.roll(nSW, -1, axis=0)
    # axis 1 is east-west; + direction is east
    nE = np.roll(nE,   1, axis=1)
    nNE = np.roll(nNE,  1, axis=1)
    nSE = np.roll(nSE,  1, axis=1)
    nW = np.roll(nW,  -1, axis=1)
    nNW = np.roll(nNW, -1, axis=1)
    nSW = np.roll(nSW, -1, axis=1)
    # Update left input stream
    np.put(n0.transpose(), geometry_input, n0_input_eq)
    np.put(nN.transpose(), geometry_input, nN_input_eq)
    np.put(nS.transpose(), geometry_input, nS_input_eq)
    np.put(nE.transpose(), geometry_input, nE_input_eq)
    np.put(nW.transpose(), geometry_input, nW_input_eq)
    np.put(nNE.transpose(), geometry_input, nNE_input_eq)
    np.put(nNW.transpose(), geometry_input, nNW_input_eq)
    np.put(nSE.transpose(), geometry_input, nSE_input_eq)
    np.put(nSW.transpose(), geometry_input, nSW_input_eq)



# Collide particles within each cell to redistribute velocities (can be optimized a little more):
def collide():
    global rho, ux, uy, n0, nN, nS, nE, nW, nNE, nNW, nSE, nSW

    rho = n0 + nN + nS + nE + nW + nNE + nSE + nNW + nSW
    ux = (nE + nNE + nSE - nW - nNW - nSW) / rho
    uy = (nN + nNE + nNW - nS - nSE - nSW) / rho
    ux2 = ux * ux               # pre-compute terms used repeatedly...
    uy2 = uy * uy
    u2 = ux2 + uy2
    omu215 = 1 - 1.5*u2         # "one minus u2 times 1.5"
    uxuy = ux * uy
    n0 = (1-omega)*n0 + omega * f_n * rho * omu215
    nN = (1-omega)*nN + omega * o_n * rho * (omu215 + 3*uy + 4.5*uy2)
    nS = (1-omega)*nS + omega * o_n * rho * (omu215 - 3*uy + 4.5*uy2)
    nE = (1-omega)*nE + omega * o_n * rho * (omu215 + 3*ux + 4.5*ux2)
    nW = (1-omega)*nW + omega * o_n * rho * (omu215 - 3*ux + 4.5*ux2)
    nNE = (1-omega)*nNE + omega * o_36 * rho * \
        (omu215 + 3*(ux+uy) + 4.5*(u2+2*uxuy))
    nNW = (1-omega)*nNW + omega * o_36 * rho * \
        (omu215 + 3*(-ux+uy) + 4.5*(u2-2*uxuy))
    nSE = (1-omega)*nSE + omega * o_36 * rho * \
        (omu215 + 3*(ux-uy) + 4.5*(u2-2*uxuy))
    nSW = (1-omega)*nSW + omega * o_36 * rho * \
        (omu215 + 3*(-ux-uy) + 4.5*(u2+2*uxuy))


# Impose boundary conditions
def boundary():
    global rho, ux, uy, n0, nN, nS, nE, nW, nNE, nNW, nSE, nSW

    # Bounce-back
    nNc = nN.copy()
    nNWc = nNW.copy()
    nWc = nW.copy()
    nSWc = nSW.copy()

    nN[wall] = nS[wall]
    nNW[wall] = nSE[wall]
    nW[wall] = nE[wall]
    nSW[wall] = nNE[wall]
    nS[wall] = nNc[wall]
    nE[wall] = nWc[wall]
    nNE[wall] = nSWc[wall]
    nSE[wall] = nNWc[wall]

    # Force steady rightward flow at left side (no need to set 0, N, and S components):
    nE[:, 0] = o_n * (1 + 3*u0 + 4.5*u0**2 - 1.5*u0**2)
    nW[:, 0] = o_n * (1 - 3*u0 + 4.5*u0**2 - 1.5*u0**2)
    nNE[:, 0] = o_36 * (1 + 3*u0 + 4.5*u0**2 - 1.5*u0**2)
    nSE[:, 0] = o_36 * (1 + 3*u0 + 4.5*u0**2 - 1.5*u0**2)
    nNW[:, 0] = o_36 * (1 - 3*u0 + 4.5*u0**2 - 1.5*u0**2)
    nSW[:, 0] = o_36 * (1 - 3*u0 + 4.5*u0**2 - 1.5*u0**2)

# Helper functions
def curl(ux, uy):
    return np.roll(uy, -1, axis=1) - np.roll(uy, 1, axis=1) - np.roll(ux, -1, axis=0) + np.roll(ux, 1, axis=0)


def mag(ux, uy):
    return np.sqrt(ux**2+uy**2)

#### Graphics helper stuff here
fig = plt.figure(figsize=(10,4))
vis = [mag, 0.2]
# vis = [curl, 0.02]
# fluidImage = plt.imshow(vis[0](ux, uy), origin='lower', norm=plt.Normalize(-vis[1], vis[1]), cmap=plt.get_cmap('jet'), interpolation='none')
fluidImage = plt.imshow(rho, origin='lower', cmap=plt.get_cmap('Reds'), interpolation='none')  
wImageArray = np.zeros((height, width, 4), np.uint8)  # an RGBA image
wImageArray[wall, 3] = 255                            # set alpha=255 wall sites only
wallImage = plt.imshow(wImageArray, origin='lower', interpolation='none')

# Print iterations progress
def printProgressBar (iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '█', printEnd = "\r"):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
        printEnd    - Optional  : end character (e.g. "\r", "\r\n") (Str)
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print('\r%s |%s| %s%% %s' % (prefix, bar, percent, suffix), end = printEnd)
    # Print New Line on Complete
    if iteration == total: 
        print()

plt.ion()
# Function called to update plot -> also progresses the simulation
def nextFrame(frame):
    plt.close
    
    # Progress our simulation
    for step in range(1): # adjust number of steps for smooth animation
        stream()
        collide()
        boundary()

    if not os.path.exists("pics"):
        os.mkdir("pics")

    plt.savefig("./pics/" + str(frame) + ".png")
    plt.close()
    
    # Show the speed in the image
    # fluidImage = plt.imshow(vis[0](ux, uy), origin='lower', norm=plt.Normalize(-vis[1], vis[1]), 
    # cmap=plt.get_cmap('jet'), interpolation='none')

    # Show the pressure in the image
    fluidImage = plt.imshow(rho, origin='lower', cmap=plt.get_cmap('Reds'), interpolation='none')
    
    # Add the walls to the image
    plt.imshow(wImageArray, origin='lower', interpolation='none')
    
    # Show progress of simulation
    printProgressBar(frame + 1, frames, prefix = 'Progress:', suffix = 'Complete', length = 50)
    return (fluidImage, wallImage)

for i in range(frames):
    nextFrame(i)

file_array = []
for filename in os.listdir("./pics"):
    file_array.append(int(filename.split('.')[0]))

file_array.sort()
file_array = [str(frame_nr) + ".png" for frame_nr in file_array]
frame = cv2.imread(os.path.join("./pics", file_array[1]))
height, width, layers = frame.shape
 
video = cv2.VideoWriter(videoname, 0, fps, (width,height))

for image in file_array:
    video.write(cv2.imread(os.path.join("./pics", image)))

cv2.destroyAllWindows()
video.release()