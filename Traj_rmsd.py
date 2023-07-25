import numpy as np
import MDAnalysis as mda
import matplotlib.pyplot as plt
from MDAnalysis.analysis import align,rms

# Make Universe object from traj and parm
md_uni = mda.Universe("../4YUS_dark_equil_strip.prmtop",
		      "../4YUS_dark_equil_strip.nc")

# set traj frame 0
md_uni.trajectory[0]

# alignment trajectory
alignment = align.AlignTraj(mobile = md_uni,
		   	    reference = md_uni,
			    select = "name CA",
			    in_memory = True)
alignment.run()

# RMSD
RMSD_analysis = rms.RMSD(md_uni,
			 md_uni,
			 select = "resid 104 to 125",
			 ref_frame = 0
			 )
RMSD_analysis.run()

time = RMSD_analysis.results.rmsd[1:, 1]
RMSD = RMSD_analysis.results.rmsd[1:, 2]

plt.plot(time, RMSD)

plt.title("RMSD 3 alpha herix")
plt.xlabel("Time (ps)")
plt.ylabel("RMSD [$\AA$]")
plt.show()
