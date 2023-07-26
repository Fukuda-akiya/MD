import numpy as np
import MDAnalysis as mda
import matplotlib.pyplot as plt
from MDAnalysis.analysis import align,rms
import seaborn as sb

# Make Universe object from traj and parm
md_uni01 = mda.Universe("../4YUS_dark_run_strip.prmtop",
	  	        "../01/4YUS_dark_run_strip01.nc")

md_uni02 = mda.Universe("../4YUS_dark_run_strip.prmtop",
	  	        "../02/4YUS_dark_run_auto02.nc")

md_uni03 = mda.Universe("../4YUS_dark_run_strip.prmtop",
	  	        "../03/4YUS_dark_run_auto03.nc")

# load average structure pdb
ave_str01 = mda.Universe("../01/4YUS_dark_run_average01.pdb")
ave_str02 = mda.Universe("../02/4YUS_dark_run_average02.pdb")
ave_str03 = mda.Universe("../03/4YUS_dark_run_average03.pdb")

# set traj frame 0
md_uni01.trajectory[0]
md_uni02.trajectory[0]
md_uni03.trajectory[0]

# alignment trajectory1
alignment = align.AlignTraj(mobile = md_uni01,
		   	    reference = ave_str01,
			    select = "name CA",
			    in_memory = True)
alignment.run()

# RMSF
C_alpha01 = md_uni01.select_atoms("name CA")
RMSF_analysis01 = rms.RMSF(C_alpha01)
RMSF_analysis01.run()


# alignment trajectory2
alignment = align.AlignTraj(mobile = md_uni02,
		   	    reference = ave_str02,
			    select = "name CA",
			    in_memory = True
			   )
alignment.run()

# RMSF
C_alpha02 = md_uni02.select_atoms("name CA")
RMSF_analysis02 = rms.RMSF(C_alpha02)
RMSF_analysis02.run()

# alignment trajectory3
alignment = align.AlignTraj(mobile = md_uni03,
		   	    reference = ave_str03,
			    select = "name CA",
			    in_memory = True
			   )
alignment.run()

# RMSF
C_alpha03 = md_uni03.select_atoms("name CA")
RMSF_analysis03 = rms.RMSF(C_alpha03)
RMSF_analysis03.run()

# load resid
resid = C_alpha01.resnums

# plot
plt.plot(resid, RMSF_analysis01.results.rmsf, lw = 0.5, label ="traj1")
plt.plot(resid, RMSF_analysis02.results.rmsf, lw = 0.5, label ="traj2")
plt.plot(resid, RMSF_analysis03.results.rmsf, lw = 0.5, label ="traj3")

plt.title("RMSF C alpha")
plt.xlabel("residue number")
plt.ylabel("RMSF [$\AA$]")
plt.legend(bbox_to_anchor=(1, 1), loc='upper right', borderaxespad=1, fontsize=12)

plt.show()
