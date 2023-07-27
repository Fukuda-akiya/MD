import numpy as np
import MDAnalysis as mda
import matplotlib.pyplot as plt
from MDAnalysis.analysis import align
import MDAnalysis.analysis.pca as pca
import pandas as pd
from sklearn.cluster import KMeans

# Make Universe object from traj and parm
md_uni = mda.Universe("/Users/fukuda/Kenkyu/OaPAC/dark_4YUS/Run/4YUS_dark_run_BB.prmtop",
	  	      "/Users/fukuda/Kenkyu/OaPAC/dark_4YUS/Run/4YUS_dark_run_all.nc")

# load average structure pdb
ave_str = mda.Universe("/Users/fukuda/Kenkyu/OaPAC/dark_4YUS/Run/4YUS_dark_run_rep_all.pdb")

# set traj frame 0
md_uni.trajectory[0]

# load time
frame_time = md_uni.trajectory.dt

# alignment trajectory1
alignment = align.AlignTraj(mobile = md_uni,
		   	    reference = ave_str,
			    select = "backbone",
			    in_memory = True)
alignment.run()

# PCA of CA
Ca_pca = pca.PCA(md_uni,
		 select = ' name CA '
		)
Ca_pca.run()

# decrease dimension for CA
transformed = Ca_pca.transform(md_uni.select_atoms("name CA"), n_components=3)
print(transformed)

# Store to dataframe
df = pd.DataFrame(transformed, columns=['PC1', 'PC2', 'PC3'])
df['Time (ps)'] = df.index * frame_time
df.head()

# k-means methods
kmeans = KMeans(n_clusters=3, max_iter=30, init="random")
cluster = kmeans.fit_predict(transformed[:,0:1])

df["cluster"] = cluster

# plot

df.plot(kind="scatter", x=0,y=1,c="cluster", cmap="winter")
#plt.scatter(x, y)
plt.title("Clustering in PCA of CA")
plt.xlabel("1st principal component[$\AA$]")
plt.ylabel("2nd principal component[$\AA$]")
plt.show()
