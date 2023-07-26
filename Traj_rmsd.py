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

# RMSD
RMSD_analysis01 = rms.RMSD(md_uni01,
			   md_uni01,
			   select = "protein",
			   ref_frame = 0
			  )
RMSD_analysis01.run()


# alignment trajectory2
alignment = align.AlignTraj(mobile = md_uni02,
		   	    reference = ave_str02,
			    select = "name CA",
			    in_memory = True
			   )
alignment.run()

# RMSD
RMSD_analysis02 = rms.RMSD(md_uni02,
			   md_uni02,
			   select = "protein",
			   ref_frame = 0
			  )
RMSD_analysis02.run()

# alignment trajectory3
alignment = align.AlignTraj(mobile = md_uni03,
		   	    reference = ave_str03,
			    select = "name CA",
			    in_memory = True
			   )
alignment.run()

# RMSD
RMSD_analysis03 = rms.RMSD(md_uni03,
			   md_uni03,
			   select = "protein",
			   ref_frame = 0
			  )
RMSD_analysis03.run()

# plot
"""
time = RMSD_analysis01.results.rmsd[0:, 1]
RMSD01 = RMSD_analysis01.results.rmsd[0:, 2]
RMSD02 = RMSD_analysis02.results.rmsd[0:, 2]
RMSD03 = RMSD_analysis03.results.rmsd[0:, 2]

plt.plot(time, RMSD01, lw = 0.3, label ="traj1")
plt.plot(time, RMSD02, lw = 0.3, label ="traj2")
plt.plot(time, RMSD03, lw = 0.3, label ="traj3")

plt.title("RMSD AC domain")
plt.xlabel("Time [microsecond]")
plt.ylabel("RMSD [$\AA$]")
plt.legend(bbox_to_anchor=(1, 0), loc='lower right', borderaxespad=1, fontsize=12)

"""
sb.kdeplot(RMSD_analysis01.results.rmsd[0:, 2], shade=True, label = "traj1")
sb.kdeplot(RMSD_analysis02.results.rmsd[0:, 2], shade=True, label = "traj2")
sb.kdeplot(RMSD_analysis03.results.rmsd[0:, 2], shade=True, label = "traj3")
plt.legend(bbox_to_anchor=(1, 1), loc='upper right', borderaxespad=0, fontsize=12)
plt.title("distribution of the RMSD")
plt.xlabel("RMSD [$\AA$]")

plt.show()
