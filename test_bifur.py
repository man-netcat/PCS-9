
import matplotlib.pyplot as plt
import numpy as np
import geometry.bifurcatie as bifur
from geometry.veins import Veins

split_vein, v0, v1, v2 = bifur.build(1000, 500, 200, angle=10)

v1.add_narrowing(0.5, 0.4, 0.5)


image = split_vein.get_image()

probe1 = v1.get_probe_point(.3)

fig = plt.figure()
ax = fig.add_subplot(111)
ax.imshow(image, cmap=plt.cm.gray)
ax.scatter(*probe1, s=20)
plt.show()
