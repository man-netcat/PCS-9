import argparse
import numpy as np
import matplotlib.pyplot as plt

import geometry.bifurcation as bifurcation

parser = argparse.ArgumentParser(description='Generate multiple bifurcation geom files')
parser.add_argument('--num', '-n', type=int, default=10, help='Number of geometries generated')
parser.add_argument('--min', '-m', type=float, default=0.0, help='Min width of the narrowing of one of the vein bifurcations')
parser.add_argument('--max', '-M', type=float, default=0.9, help='Max width of the narrowing of one of the vein bifurcations')
args = parser.parse_args()
print(args.min, args.max, args.num)

widths = np.linspace(args.min, args.max, num=args.num)
split_vein, v0, v1, v2 = None, None, None, None
for width in widths:
    split_vein, v0, v1, v2 = bifurcation.build(1000, 500, 200, angle=20)
    v1.add_narrowing(loc=0.5, width=0.4, height=width)
    image = split_vein.get_image()

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
