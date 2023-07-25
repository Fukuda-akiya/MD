import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from math import pi

# load trajctory and parmtop
traj = md.load("/Users/pride/kenkyu/OaPAC/MD/dark/4YUS/4YUS_dark_equil_strip.nc",
		top="/Users/pride/kenkyu/OaPAC/MD/dark/4YUS/4YUS_dark_equil_strip.prmtop")

# load flame of each time in traj
time = traj.time

psi_indices = [108,110,111,112]
#psi_indices = [756,759,761,763]
#psi_indices = [5962,5964,5965,5966]
#psi_indices = [6610,6613,6615,6617]
angles = md.compute_dihedrals(traj, [psi_indices])

# output file
fo = open("test.dat","w")

psi = [0] * len(time)

for i in range(len(time)):
	psi[i] = (angles[i]*180)/pi
	fo.write(str(time[i]*0.1 - 531) + "\t" + str(*psi[i]) + "\n")

fo.close()

# plot
plt.plot(traj.time*0.1 - 531, psi, label="Dihedrals")
plt.xlabel("time [frame]")
plt.xlim(0,9999)
plt.ylabel("angle [degree]")
plt.show()
