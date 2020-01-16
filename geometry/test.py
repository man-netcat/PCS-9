from geometry.generator import Veins
import matplotlib.pyplot as plt
import numpy as np

split_vein = Veins(10)
v0 = split_vein.add_vein((0, 250), (200, 250), width=30)
v0.append_vein((1000, 100), width=30)
v0.append_vein((1000, 400), width=20)
# h = veins(res=10)
# v0 = h.new_vein(...)
# v1 = v0.appendVein(length, angle, width, bevel_start=True, bevel_end=False)
# v2 = v0.appendVeinAt(pos, width, ...)
# h.addInlet(flow)
# h.addOutlet()

# h.size()
# h.get_inlets() --> {postitions: [], flow:number}[]
# h.get_outlets() --> {postitions: []}

image = split_vein.get_image()
fig = plt.figure()
ax = fig.add_subplot(111)
ax.imshow(image, cmap=plt.cm.gray)
plt.show()
print(image)