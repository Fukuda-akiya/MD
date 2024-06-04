import numpy as np
import MDAnalysis as mda
import matplotlib.pyplot as plt
from MDAnalysis.analysis import align
import MDAnalysis.analysis.pca as pca
import pandas as pd
from sklearn.cluster import KMeans
import seaborn as sns
import nglview as nv

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
		 select = 'backbone '
		)
Ca_pca.run()

# decrease dimension for CA
transformed = Ca_pca.transform(md_uni.select_atoms("backbone"), n_components=3)

# PC1
pc_1 = Ca_pca.results.p_components[:,2]

weight_1 = transformed[:,2]

mean_coord = Ca_pca.mean

pro_1 = np.outer(weight_1, pc_1) + mean_coord.flatten()

coordinates = pro_1.reshape(len(weight_1),-1,3)
print(coordinates.shape)

# Make universe 
pro_uni = mda.Merge(md_uni.select_atoms("backbone"))
pro_uni.load_new(coordinates, order = "fac")
pro_select = pro_uni.select_atoms("backbone")

# traj output
pro_select.write('PC3_backbone_all.nc', frames="all")

# View NGL
#view = nv.show_mdanalysis(pro_uni.atoms) 
#view

"""
# k-means methods
kmeans = KMeans(n_clusters=3, max_iter=30, init="random")
cluster = kmeans.fit_predict(transformed[:,0:1:2])
df = pd.DataFrame(transformed, columns=['PC1', 'PC2', 'PC3'])
df["cluster"] = cluster
#centers = kmeans.cluster_centers_
#df_cluster_centers = pd.DataFrame(centers, columns=["x", "y"])


# plot clustering
sns.scatterplot(data=df, x='PC1',y='PC2',hue="cluster")
# plot cluster center
#sns.scatterplot(data=df_cluster_centers, x="x", y="y", s=200, marker='*', color='gold', linewidth=0.5)
plt.scatter(kmeans.cluster_centers_[:,0],kmeans.cluster_centers_[:,1], s=200, marker='*', color='gold', linewidth=0.5)

plt.title("Clustering in PCA of CA")
#plt.xlabel("1st principal component[$\AA$]")
#plt.ylabel("2nd principal component[$\AA$]")
plt.show()
"""
