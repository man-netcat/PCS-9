import source.geometry.aorta as aorta
import matplotlib.pyplot as plt

abdominal, veins = aorta.build_abdominal(20)
print(type(veins))
image = abdominal.get_image()
fig = plt.figure()
ax = fig.add_subplot(111)
ax.imshow(image, cmap=plt.cm.gray)
plt.show()
