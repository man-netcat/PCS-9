import geometry.aorta as aorta
import matplotlib.pyplot as plt

abdominal = aorta.build_abdominal(20)

image = abdominal.get_image()
fig = plt.figure()
ax = fig.add_subplot(111)
ax.imshow(image, cmap=plt.cm.gray)
plt.show()
