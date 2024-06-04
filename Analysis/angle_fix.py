import numpy as np

frame_angle = []
frame_angle_save = []

with open("analysis_files/a3a-a3b_light_apo_angle.dat","r") as fin:
     next(fin) # title read skip
     frame_angle = fin.readlines()

angle = [0] * len(frame_angle)
angle_save = [0] * len(frame_angle)
frame = [0] * len(frame_angle)

print("frame" + "\t" + "angle")
for i in range(len(frame_angle)):
    frame_angle_save = frame_angle[i]
    angle[i] = float(frame_angle_save.split()[1])
    frame[i] = int(frame_angle_save.split()[0])
    
    if angle[i] > 90.0:
        print(str(frame[i]) + "\t" + "{:.6}".format(180.0 - angle[i]))
    else:
        print(str(frame[i]) + "\t" + "{:.6}".format(angle[i]))
