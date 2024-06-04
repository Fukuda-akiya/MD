import numpy as np
import mdtraj as md
import matplotlib.pyplot as plt
from contact_map import ContactFrequency, AtomMismatchedContactDifference
import pandas as pd

# resd trajctory and topology
traj1 = md.load("/Users/fukuda/Kenkyu/OaPAC/dark_4YUS/Run/4YUS_dark_run_strip_comp10.nc",
		top = "/Users/fukuda/Kenkyu/OaPAC/dark_4YUS/Run/4YUS_dark_run_strip.prmtop")
traj2 = md.load("/Users/fukuda/Kenkyu/OaPAC/dark_4YUT/4YUT_dark_run_strip_all.nc",
		top = "/Users/fukuda/Kenkyu/OaPAC/dark_4YUT/4YUT_dark_run_strip.prmtop")
topology1 = traj1.topology	
topology2 = traj2.topology

# calculate contact
resid1 = topology1.select("name CA")
protein1 = topology1.select("resSeq 1 to 350 or resSeq 367 to 717 and symbol != 'H'")
protein2 = topology2.select("resSeq 1 to 350 or resSeq 368 to 718 and symbol != 'H'")
traj_contacts1 = ContactFrequency(traj1, query=protein1, haystack=protein1)
traj_contacts2 = ContactFrequency(traj2, query=protein2, haystack=protein2)
diff = AtomMismatchedContactDifference(traj_contacts1,traj_contacts2)

f = open("diff_contact_frequency.dat","w")
f.write(str(diff.residue_contacts.most_common()[:15]) + "\n")
f.write("\n")
f.write(str(list(reversed(diff.residue_contacts.most_common()))[:15]) + "\n")
f.close()

# plot
fig, ax = diff.residue_contacts.plot(cmap='seismic', vmin=-1, vmax=1)
colorbar_obj = fig.axes[1]
colorbar_obj.yaxis.set_ticks(np.arange(-1.0,1.1,0.25))
colorbar_obj.set_yticklabels(np.arange(-1.0,1.1,0.25), weight='bold',  size=13, fontname='Arial')
ax.set_aspect('equal', adjustable='box')
plt.xlabel("Residue number", fontsize=18, fontname='Arial', fontweight='bold')
plt.ylabel("Residue number", fontsize=18, fontname='Arial', fontweight='bold')
plt.xticks(fontsize=13, fontname='Arial', fontweight='bold')
plt.yticks(fontsize=13, fontname='Arial', fontweight='bold')
plt.savefig('Png/Diff_contact_map.png', bbox_inches='tight', pad_inches = 0, dpi=300)
plt.show()
