import argparse
import os
import numpy as np
import matplotlib.pyplot as plt
from PIL import Image

import geometry.bifurcation as bifurcation

# parse arguments
parser = argparse.ArgumentParser(description='Generate multiple bifurcation geom files')
parser.add_argument('--width', '-w', type=int, default=1000, help='Width of image')
parser.add_argument('--height', '-h', type=int, default=500, help='Height of image')
parser.add_argument('--num', '-n', type=int, default=10, help='Number of geometries generated')
parser.add_argument('--min', '-m', type=float, default=0.0, help='Min width of the narrowing of one of the vein bifurcations')
parser.add_argument('--max', '-M', type=float, default=0.9, help='Max width of the narrowing of one of the vein bifurcations')
parser.add_argument('--out', '-o', type=str, default='out', help='The name of the output folder')
parser.add_argument('--clear', '-c', action='store_true', help='Clear the output folder')
args = parser.parse_args()

# make folder
if not args.out:
    print('No out path')
    exit(1)
folder = args.out
if not os.path.exists(folder):
    os.mkdir(folder)
elif args.clear:
    os.rmdir(folder)
    os.mkdir(folder)

# generate
widths = np.linspace(args.min, args.max, num=args.num)
split_vein, v0, v1, v2 = None, None, None, None
for width in widths:
    split_vein, v0, v1, v2 = bifurcation.build(args.width, args.height, args.width/5, angle=20)
    v1.add_narrowing(loc=0.5, width=0.4, height=width)
    image = split_vein.get_image()
    filename = f'bifurcation_{width:9.4f}.png'
    new_im = Image.fromarray(image)
    new_im.save(filename)

if not v0:
    print('Nothing generated')
    exit(1)

# display
image = split_vein.get_image()

probe_in = v0.get_probe_point(0.3)
probe_normal = v2.get_probe_point(0.8)
probe_nar_before = v1.get_probe_point(0.2)
probe_nar_after = v1.get_probe_point(0.8)

print(probe_in, probe_normal, probe_nar_before, probe_nar_after)
fig = plt.figure()
ax = fig.add_subplot(111)
ax.imshow(image, cmap=plt.cm.gray)
ax.scatter(*zip(*[probe_in, probe_normal, probe_nar_before, probe_nar_after]), s=20)
plt.show()
