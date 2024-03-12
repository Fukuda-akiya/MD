import matplotlib.pyplot as plt
import numpy as np
import seaborn as sb
import pandas as pd

chi_number = 1
# chain A
name = "BLUF/dihedral_light_holo_Phe47"
dataA = np.genfromtxt( name + ".dat" )
dataA_save = np.zeros((len(dataA), chi_number))

# chain B
name = "BLUF/dihedral_light_holo_Phe413"
dataB = np.genfromtxt( name + ".dat" )
dataB_save = np.zeros((len(dataB), chi_number))

# reference load

name = "/home/users/dzn/Kenkyu/OaPAC/Analysis/Mutation/dark/holo/dihedral/BLUF/dihedral_dark_holo_Phe47_ref"
refA = np.genfromtxt( name + ".dat" )
name = "/home/users/dzn/Kenkyu/OaPAC/Analysis/Mutation/dark/holo/dihedral/BLUF/dihedral_dark_holo_Phe413_ref"
refB = np.genfromtxt( name + ".dat" )

for i in range(len(dataA)):
    for j in range(0,chi_number):
        if dataA[i][j+1] < -180:
            dataA_save[i][j] = dataA[i][j+1] + 360
        elif dataA[i][j+1] > 180:
            dataA_save[i][j] = dataA[i][j+1] - 360
        else:
            dataA_save[i][j] = dataA[i][j+1]

for i in range(len(dataB)):
    for j in range(0,chi_number):
        if dataB[i][j+1] < -180:
            dataB_save[i][j] = dataB[i][j+1] + 360
        elif dataB[i][j+1] > 180:
            dataB_save[i][j] = dataB[i][j+1] - 360
        else:
            dataB_save[i][j] = dataB[i][j+1]

plt.rcParams["figure.figsize"] = (8, 6)
"""
for i in range(0,chi_number):
    sb.kdeplot(dataA_save[:,i], lw = 1.0, label="chi" + str(i+1), alpha = 1.0)
"""
sb.kdeplot(dataA_save[:,0], lw = 1.0, label="chi" + str(1+1), alpha = 1.0, color = "blue")
sb.kdeplot(dataB_save[:,0], lw = 1.0, label="chi" + str(1+1), alpha = 1.0, color = "darkorange")

plt.axvline(refA[1], lw = 0.8, color = "blue", label = "Initial Angle2")
plt.axvline(refB[1], lw = 0.8, color = "darkorange", label = "Initial Angle2")
plt.title("Phe47", fontsize=16, fontweight='bold') 
plt.ylabel("ratio", fontsize=16, fontweight='bold')
plt.xlabel("angle [degree]", fontsize=16, fontweight='bold')
plt.xticks(fontsize=11, fontweight='bold')
plt.yticks(fontsize=11, fontweight='bold')
plt.xlim(-180,180)
plt.ylim(0,0.03)
outputname = "dihedral_light_holo_Phe47"
plt.savefig( "BLUF/png/" + outputname + "_distribution.png")
