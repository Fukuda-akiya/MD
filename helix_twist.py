import numpy as np
import math
import MDAnalysis as mda
from MDAnalysis.analysis import align
import matplotlib.pyplot as plt

np.set_printoptions(threshold=np.inf)

#########################################
###############Function##################
#########################################

# Universeから 座標をndarrayに取り出す関数
def extract(name):
    a3helix = mda_uni.select_atoms(name)
    coord_from_mda = np.array([a3helix.positions for frame in mda_uni.trajectory])

    return coord_from_mda

# 原子数の取得
def atom_number(name):
    a3helix = mda_uni.select_atoms(name)

    return a3helix

# 重心を計算する関数
def center(coord,atom):
    # set ndarray for center
    center = np.zeros([len(coord),3])
    for i in range(len(coord)):
        for j in range(len(atom)):
            for k in range(3):
                center[i][k] += coord[i,j,k] / len(atom)

    return center

# 法線ベクトルを作成する関数
def normal(theta,coord,R):
    rot = np.matrix([
                    [np.cos(theta), -np.sin(theta), 0],
                    [np.sin(theta),  np.cos(theta), 0],
                    [0, 0, 1],
                      ])
    n_vector = np.zeros([len(coord),3])
    for i in range(len(coord)):
        for j in range(3):
            n_vector[i][j] = rot[j] @ R[i]
    
    # 規格化
    n_norm = np.zeros([len(coord),3])
    for i in range(len(coord)):
        norm = np.linalg.norm(n_vector[i])
        for j in range(3):
            n_norm[i][j] = n_vector[i][j] / norm

    return n_norm

########################################
########################################
########################################

# make universe
# parmtopとtrajectoryの読み込み
mda_uni = mda.Universe("../parm_files/4YUS_dark_run_for_secstruct.prmtop",
                       "../traj_files/4YUS_dark_run_for_secstruct.nc")

# Referenceの読み込み
ave_str = mda.Universe("../pdb_files/ref_for_secstruct.pdb")

# trajectoryを初期frameにセット
mda_uni.trajectory[0]

# Align trajectory
alignment = align.AlignTraj(mobile = mda_uni,
                            reference = ave_str,
                            select = "resid 146 to 366 512 to 733 ",
                            in_memory = True)
alignment.run()

# MDAnalysisのUniverseオブジェクトから座標をnumpyのarrayに取り出す
coord_from_mdaA_bottom = extract("resid 103:107")
coord_from_mdaA_top = extract("resid 121:125")
coord_from_mdaB_bottom = extract("resid 469:473")
coord_from_mdaB_top = extract("resid 487:491")

atom_number_bottom = atom_number("resid 103:107")
atom_number_top = atom_number("resid 121:125")

# 重心の計算
centerA_bottom = np.zeros(len(coord_from_mdaA_bottom))
centerA_bottom = center(coord_from_mdaA_bottom,atom_number_bottom)
centerA_top = np.zeros(len(coord_from_mdaA_top))
centerA_top = center(coord_from_mdaA_top,atom_number_top)

centerB_bottom = np.zeros(len(coord_from_mdaB_bottom))
centerB_bottom = center(coord_from_mdaB_bottom,atom_number_bottom)
centerB_top = np.zeros(len(coord_from_mdaB_top))
centerB_top = center(coord_from_mdaB_top,atom_number_top)

# 重心間のベクトルの計算
R_bottom = np.zeros([len(coord_from_mdaA_bottom),3])
R_top = np.zeros([len(coord_from_mdaA_top),3])

for i in range(len(coord_from_mdaA_bottom)):
    for j in range(3):
        R_bottom[i][j] = centerA_bottom[i][j] - centerB_bottom[i][j]
        R_top[i][j] = centerA_top[i][j] - centerB_top[i][j]

# 方線ベクトルの生成（回転行列の作用）
theta = np.pi * 0.5
n_vector_bottom = normal(theta,coord_from_mdaA_bottom,R_bottom)
n_vector_top = normal(theta,coord_from_mdaA_top,R_top)

#for i in range(len(coord_from_mdaA_bottom)):
#    print(str(i+1) + "\t" + str(np.linalg.norm(n_vector_bottom[i])))

twist = np.zeros(len(coord_from_mdaA_bottom))
frame = np.zeros(len(coord_from_mdaA_bottom))
for i in range(len(coord_from_mdaA_bottom)):
    twist[i] = np.arcsin( np.linalg.norm( np.cross( n_vector_bottom[i],n_vector_top[i] ) ) ) * 180 / np.pi
    frame[i] = i + 1

# plot
plt.rcParams["figure.figsize"] = (8, 6)
plt.title("WT Dark", fontweight = "bold")
plt.plot(frame, twist, lw = 1.25, color = "blue")
plt.xlabel("Frame", fontweight = "bold", fontsize = 12)
plt.ylabel("Twist [Degree]", fontweight = "bold", fontsize = 12)
plt.xticks(fontweight = "bold", fontsize = 11)
plt.yticks(fontweight = "bold", fontsize = 11)
plt.xlim(0,9000)
plt.ylim(35,80)
plt.savefig("Twist_WT_dark_holo_3us.png")
