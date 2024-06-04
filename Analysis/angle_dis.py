import matplotlib.pyplot as plt
import numpy as np
import seaborn as sb

png_name = "angle_distribution_all.png"

# Fix Angle Function
def fix_angle(data):
    angle = np.zeros(len(data))
    for i in range(len(data)):
        if data[i,1] > 90.0:
            angle[i] = 180.0 - data[i,1]
        else:
            angle[i] = data[i,1]

    return angle

# WT dark
WT_dark = "a3a-a3b_WT_dark_angle.dat"
data01 = np.genfromtxt(WT_dark)
angle01 = fix_angle(data01)

# WT light
WT_light = "a3a-a3b_WT_light_angle.dat"
data02 = np.genfromtxt(WT_light)
angle02 = fix_angle(data02)

# L111A/L115A
Mut_light = "a3a-a3b_L111A-L115A_light_angle.dat"
data03 = np.genfromtxt(Mut_light)
angle03 = fix_angle(data03)

# Trp90/Met92
switch = "a3a-a3b_Trp90-Met92_angle.dat"
data04 = np.genfromtxt(switch) 
angle04 = fix_angle(data04)

sb.kdeplot(angle01[:], fill=True, label = "WT dark")
sb.kdeplot(angle02[:], fill=True, label = "WT light")
sb.kdeplot(angle03[:], fill=True, label = "L111A/L115A")
sb.kdeplot(angle04[:], fill=True, label = "Trp90/Met92")
plt.xticks(fontweight='bold')
plt.yticks(fontweight='bold')
plt.xlabel("Angle [degree]",fontweight = "bold")
plt.ylabel("Density", fontweight = "bold")
plt.legend(bbox_to_anchor=(1,1), loc="upper right", borderaxespad=0, fontsize=12) 
plt.savefig(png_name)
