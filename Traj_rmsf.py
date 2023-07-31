import numpy as np
import MDAnalysis as mda
import matplotlib.pyplot as plt
from MDAnalysis.analysis import align,rms
import seaborn as sb
import pandas as pd

# Make Universe object from traj and parm
md_uni = mda.Universe("/Users/fukuda/Kenkyu/OaPAC/dark_4YUS/Run/4YUS_dark_run_BB_C_strip.prmtop",
	  	        "/Users/fukuda/Kenkyu/OaPAC/dark_4YUS/Run/4YUS_dark_run_BB_C_strip_all.nc")

# load average structure pdb
ave_str = mda.Universe("/Users/fukuda/Kenkyu/OaPAC/dark_4YUS/Run/4YUS_dark_run_rep_all.pdb")

# set traj frame 0
md_uni.trajectory[0]

# alignment trajectory1
alignment = align.AlignTraj(mobile = md_uni,
		   	    reference = ave_str,
			    select = "backbone",
			    in_memory = True)
alignment.run()

# RMSF
C_alpha01 = md_uni.select_atoms("name CA")
RMSF_analysis01 = rms.RMSF(C_alpha01)
RMSF_analysis01.run()

# load resid
resid = C_alpha01.resnums

# output
f = open("RMSF.dat","w")
for i in range(len(RMSF_analysis01.results.rmsf)):
	f.write(str(resid[i]) + "\t" + str(RMSF_analysis01.results.rmsf[i]) + "\n")

f.close()

# Make dataframe
df = pd.DataFrame(
	data={"Resid":resid,
	      "RMSF":RMSF_analysis01.results.rmsf}
	      )

df_s = df.sort_values('RMSF')

pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)

f = open("RMSF_data.dat","w")
f.write(str(df_s))
f.close()

"""
# plot
plt.plot(resid, RMSF_analysis01.results.rmsf, lw = 0.5)

plt.title("RMSF C alpha")
plt.xlabel("residue number")

plt.show()
"""
