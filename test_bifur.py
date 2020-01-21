
import matplotlib.pyplot as plt
import numpy as np
import geometry.bifurcatie as bifur
from geometry.veins import Veins

split_vein = bifur.build(1000, 500, 50, angle=30)
image = split_vein.get_image()
fig = plt.figure()
ax = fig.add_subplot(111)
ax.imshow(image, cmap=plt.cm.gray)
plt.show()
