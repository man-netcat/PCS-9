import cv2
import numpy as np
import os
import sys

fps = 10
name = "video"
if len(sys.argv) > 1:
    fps = int(sys.argv[1])
if len(sys.argv) > 2:
    name = sys.argv[2]


file_array = []
for filename in os.listdir("./pics"):
    file_array.append(int(filename.split('.')[0]))

file_array.sort()
file_array = [str(frame_nr) + ".jpg" for frame_nr in file_array]
frame = cv2.imread(os.path.join("./pics", file_array[0]))
height, width, layers = frame.shape
 
video = cv2.VideoWriter(name + ".avi", 0, fps, (width,height))

for image in file_array:
    video.write(cv2.imread(os.path.join("./pics", image)))

cv2.destroyAllWindows()
video.release()