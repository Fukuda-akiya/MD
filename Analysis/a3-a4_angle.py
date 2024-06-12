import numpy as np
import math
import MDAnalysis as mda
import matplotlib.pyplot as p

np.set_printoptions(threshold=np.inf)

# 次元数(x,y,z etc.)
dim = 3
# 読み込むparmtop,trajectory
# 8qfg
traj_name = "/home/users/dzn/Kenkyu/OaPAC/Analysis/8QFG/traj/8QFG_run01-05.nc"
parm_name = "/home/users/dzn/Kenkyu/OaPAC/Analysis/8QFG/parmtop/updated_bb.prmtop"
# WT dark
#traj_name = "/home/users/dzn/Kenkyu/OaPAC/Analysis/dark/4YUS_holo/traj_files/4YUS_dark_run01-05_3us.nc"
#parm_name = "/home/users/dzn/Kenkyu/OaPAC/Analysis/dark/4YUS_holo/parm_files/4YUS_dark_holo_updated_3us.prmtop"
# select resid
resid_bottom_a3 = "resid 103:107 453:457"
resid_top_a3 = "resid 121:125 471:475"
resid_bottom_a4A = "resid 127:129" 
resid_top_a4A = "resid 134:136"
resid_bottom_a4B = "resid 477:479" 
resid_top_a4B = "resid 484:486"

#########################################
###############Function##################
#########################################

# Universeから座標をndarrayに取り出す関数
def extract(name):
    atom = mda_uni.select_atoms(name)
    coord_from_mda = np.array([atom.positions for frame in mda_uni.trajectory])

    return coord_from_mda

# 原子数の取得
def atom_number(name):
    atom = mda_uni.select_atoms(name)

    return atom

# 重心を計算する関数
def center(coord,atom):
    # set ndarray for center
    center = np.zeros([len(coord),dim])
    for i in range(len(coord)):
        for j in range(len(atom)):
            for k in range(dim):
                center[i][k] += coord[i,j,k] / len(atom)

    return center

# 重心間のベクトルを計算する関数
def diff(coord,centerA,centerB):
    R = np.zeros([len(coord),dim])
    for i in range(len(coord)):
        for j in range(dim):
            R[i][j] = centerA[i][j] - centerB[i][j]

    return R

#　角度の計算（内積）
def inner(coord, R_a, R_b):
    angle = np.zeros(len(coord))
    norm_a = np.zeros(len(coord))
    norm_b = np.zeros(len(coord))
    for i in range(len(coord)):
        # 内積の計算
        dot_product = np.dot(R_a.T[:,i], R_b[i,:])
        # ベクトルの大きさを計算
        norm_a[i] = np.linalg.norm(R_a[i,:])
        norm_b[i] = np.linalg.norm(R_b[i,:])
        # 角度
        angle[i] = np.arccos(dot_product / (norm_a[i] * norm_b[i])) * 180 / np.pi

    return angle

########################################
########################################
########################################


# parmtopとtrajectoryの読み込み
mda_uni = mda.Universe(parm_name,traj_name)

# trajectoryを初期frameにセット
mda_uni.trajectory[0]

# MDAnalysisのUniverseオブジェクトから座標をnumpyのarrayに取り出す
# a3 helix
coord_from_mda_bottom_a3 = extract(resid_bottom_a3)
coord_from_mda_top_a3 = extract(resid_top_a3)
# a4 helix
coord_from_mda_bottom_a4A = extract(resid_bottom_a4A)
coord_from_mda_top_a4A = extract(resid_top_a4A)
coord_from_mda_bottom_a4B = extract(resid_bottom_a4B)
coord_from_mda_top_a4B = extract(resid_top_a4B)

# 原子数の取得
# a3 helix
atom_number_bottom_a3 = atom_number(resid_bottom_a3)
atom_number_top_a3 = atom_number(resid_top_a3)
# a4 helix
atom_number_bottom_a4A = atom_number(resid_bottom_a4A)
atom_number_top_a4A = atom_number(resid_top_a4A)
atom_number_bottom_a4B = atom_number(resid_bottom_a4B)
atom_number_top_a4B = atom_number(resid_top_a4B)

# 重心の計算
# 配列の初期化
center_bottom_a3 = np.zeros(len(coord_from_mda_bottom_a3))
center_top_a3 = np.zeros(len(coord_from_mda_top_a3))
center_bottom_a4A = np.zeros(len(coord_from_mda_bottom_a4A))
center_top_a4A = np.zeros(len(coord_from_mda_top_a4A))
center_bottom_a4B = np.zeros(len(coord_from_mda_bottom_a4B))
center_top_a4B = np.zeros(len(coord_from_mda_top_a4B))
# a3 helix
center_bottom_a3 = center(coord_from_mda_bottom_a3,atom_number_bottom_a3)
center_top_a3 = center(coord_from_mda_top_a3,atom_number_top_a3)
# a4 helix
center_bottom_a4A = center(coord_from_mda_bottom_a4A,atom_number_bottom_a4A)
center_top_a4A = center(coord_from_mda_top_a4A,atom_number_top_a4A)
center_bottom_a4B = center(coord_from_mda_bottom_a4B,atom_number_bottom_a4B)
center_top_a4B = center(coord_from_mda_top_a4B,atom_number_top_a4B)

# 重心間のベクトルの計算
# a3 helix
R_a3 = diff(coord_from_mda_bottom_a3,center_top_a3,center_bottom_a3)
R_a4A = diff(coord_from_mda_bottom_a4A,center_top_a4A,center_bottom_a4A)
R_a4B = diff(coord_from_mda_bottom_a4B,center_top_a4B,center_bottom_a4B)

# 角度の計算
angle_a3_a4A = inner(coord_from_mda_bottom_a3,R_a3,R_a4A)
angle_a3_a4B = inner(coord_from_mda_bottom_a3,R_a3,R_a4B)

print("Frame" + "\t" + "a3-a4a" + "\t" + "a3-a4b")
for i in range(len(angle_a3_a4A)):
    print(str(i) + "\t" + str(angle_a3_a4A[i]) + "\t" + str(angle_a3_a4B[i]))
