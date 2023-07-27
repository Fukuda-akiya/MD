import numpy as np
import MDAnalysis as mda
import matplotlib.pyplot as plt
from MDAnalysis.analysis import align,rms
import seaborn as sb

# Make Universe object from traj and parm
md_uni01 = mda.Universe("/Users/fukuda/Kenkyu/OaPAC/dark_4YUS/Run/4YUS_dark_run_BB.prmtop",
	  	        "/Users/fukuda/Kenkyu/OaPAC/dark_4YUS/Run/01/4YUS_dark_run_BB.nc")

md_uni02 = mda.Universe("/Users/fukuda/Kenkyu/OaPAC/dark_4YUS/Run/4YUS_dark_run_BB.prmtop",
	  	        "/Users/fukuda/Kenkyu/OaPAC/dark_4YUS/Run/02/4YUS_dark_run_auto.nc")

md_uni03 = mda.Universe("/Users/fukuda/Kenkyu/OaPAC/dark_4YUS/Run/4YUS_dark_run_BB.prmtop",
	  	        "/Users/fukuda/Kenkyu/OaPAC/dark_4YUS/Run/03/4YUS_dark_run_auto.nc")

# load average structure pdb
ave_str01 = mda.Universe("/Users/fukuda/Kenkyu/OaPAC/dark_4YUS/Run/01/4YUS_dark_run_rep01.pdb")
ave_str02 = mda.Universe("/Users/fukuda/Kenkyu/OaPAC/dark_4YUS/Run/02/4YUS_dark_run_rep02.pdb")
ave_str03 = mda.Universe("/Users/fukuda/Kenkyu/OaPAC/dark_4YUS/Run/03/4YUS_dark_run_rep03.pdb")

# set traj frame 0
md_uni01.trajectory[0]
md_uni02.trajectory[0]
md_uni03.trajectory[0]

# alignment trajectory1
alignment = align.AlignTraj(mobile = md_uni01,
		   	    reference = ave_str01,
			    select = "backbone",
			    in_memory = True)
alignment.run()

# RMSD
RMSD_analysis01 = rms.RMSD(md_uni01,
			   ave_str01,
			   select = "resid 151 to 279 518 to 646",
			   ref_frame = 0
			  )
RMSD_analysis01.run()

#print(np.min(RMSD_analysis01.results.rmsd,axis=0))


# alignment trajectory2
alignment = align.AlignTraj(mobile = md_uni02,
		   	    reference = ave_str02,
			    select = "backbone",
			    in_memory = True
			   )
alignment.run()

# RMSD
RMSD_analysis02 = rms.RMSD(md_uni02,
			   ave_str02,
			   select = "resid 151 to 279 518 to 646",
			   ref_frame = 0
			  )
RMSD_analysis02.run()

#print(np.min(RMSD_analysis02.results.rmsd,axis=0))


# alignment trajectory3
alignment = align.AlignTraj(mobile = md_uni03,
		   	    reference = ave_str03,
			    select = "backbone",
			    in_memory = True
			   )
alignment.run()

# RMSD
RMSD_analysis03 = rms.RMSD(md_uni03,
			   ave_str03,
			   select = "resid 151 to 279 518 to 646",
			   ref_frame = 0
			  )
RMSD_analysis03.run()

#print(np.min(RMSD_analysis03.results.rmsd,axis=0))

# plot
"""
f = open("RMSD03.dat","w")
for i in range(len(RMSD_analysis03.results.rmsd)):
	f.write(str(RMSD_analysis03.results.rmsd[i]) + "\n")

f.close()
"""
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
