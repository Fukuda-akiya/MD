import numpy as np
import MDAnalysis as mda
import matplotlib.pyplot as plt
from MDAnalysis.analysis import align
import MDAnalysis.analysis.pca as pca
import pandas as pd
from sklearn.cluster import KMeans
import seaborn as sns

# Make Universe object from traj and parm
md_uni = mda.Universe("/Users/fukuda/Kenkyu/OaPAC/dark_4YUS/Run/4YUS_dark_run_BB_C_strip.prmtop",
	  	      "/Users/fukuda/Kenkyu/OaPAC/dark_4YUS/Run/4YUS_dark_run_BB_C_strip_all.nc")

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

# Store to dataframe
#df = pd.DataFrame(transformed, columns=['PC1', 'PC2', 'PC3'])
#df['Time (ps)'] = df.index * frame_time

# k-means methods
kmeans = KMeans(n_clusters=3, max_iter=30, init="random")
cluster = kmeans.fit_predict(transformed[:,0:1])
df = pd.DataFrame(transformed, columns=['PC1', 'PC2', 'PC3'])
df["cluster"] = cluster
centers = kmeans.cluster_centers_
#df_cluster_centers = pd.DataFrame(cluster_centers, columns=["x1", "x2"])
print(centers[:,0] + centers[:,1])
"""
# plot clustering
sns.scatterplot(data=df, x='PC1',y='PC2',hue="cluster")
# plot cluster center
sns.scatterplot(data=df_cluster_centers, x="x1", y="x2",  s=200, marker='*', color='gold', linewidth=0.5)

plt.title("Clustering in PCA of CA")
#plt.xlabel("1st principal component[$\AA$]")
#plt.ylabel("2nd principal component[$\AA$]")
plt.show()
"""
