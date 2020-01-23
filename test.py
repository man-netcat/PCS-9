
import matplotlib.pyplot as plt
import numpy as np

from geometry.veins import Veins

split_vein = Veins(1000, 500)
v0 = split_vein.add_vein((0, 250), (200, 450), angle_from=0, width=50)
a1 = v0.add_end(40)
a2 = v0.add_end(-40)
v0.append_vein((1000, 400), end_i=a1, width=50)
v0.append_vein((1000, 100), end_i=a2, width=20)
print('test')
# h = veins(res=10)
# v0 = h.new_vein(...)
# v1 = v0.appendVein(length, angle, width, bevel_start=True, bevel_end=False)
# v2 = v0.appendVeinAt(pos, width, ...)
# h.addInlet(flow)
# h.addOutlet()


image = split_vein.get_image()
fig = plt.figure()
ax = fig.add_subplot(111)
ax.imshow(image, cmap=plt.cm.gray)
plt.show()
# print(image)