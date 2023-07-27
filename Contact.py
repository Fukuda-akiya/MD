import numpy as np
import mdtraj as md
import matplotlib.pyplot as plt
from contact_map import ContactFrequency

traj = md.load("/Users/fukuda/Kenkyu/OaPAC/dark_4YUS/Run/4YUS_dark_run_all.nc",
		top = "/Users/fukuda/Kenkyu/OaPAC/dark_4YUS/Run/4YUS_dark_run_BB.prmtop")
topology = traj.topology

traj_contacts = ContactFrequency(traj)
fig, ax = traj_contacts.residue_contacts.plot(cmap='seismic', vmin=-1, vmax=1);
plt.show()
