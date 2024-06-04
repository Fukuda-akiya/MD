import numpy as np
import math
import MDAnalysis as mda
from MDAnalysis.analysis import align
import matplotlib.pyplot as plt

np.set_printoptions(threshold=np.inf)

##### 定数や読み込むファイル #####

# 次元数(x,y,z etc.)
dim = 3
# Runの本数
Run = 3.0
# window幅
width = 100
# 読み込むparmtop,trajectory,reference
parm_name = "../parm_files/4YUS_light_holo_updated_for_secstruct.prmtop"
traj_name = "../traj_files/4YUS_light_run_for_secstruct.nc"
ref_name  = "/home/users/dzn/Kenkyu/OaPAC/Analysis/dark/4YUS_holo/pdb_files/ref_for_secstruct.pdb"
# fitting
fitting_area = "resid 146 to 366 512 to 733" 
# select resid
residA_bottom = "resid 103:107"
residA_top = "resid 121:125"
residB_bottom = "resid 469:473"
residB_top = "resid 487:491"
# About fig
fig_title = "WT Light"
fig_name = "Twist_WT_light_holo_3us.png"
# 表示領域
ymin = 40
ymax = 70
# plot color
twist_color = "orange"

##### 基本的にここから下はいじらない #####


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

# 法線ベクトルを作成する関数
def normal(theta,coord,R):
    rot = np.matrix([
                    [np.cos(theta), -np.sin(theta), 0],
                    [np.sin(theta),  np.cos(theta), 0],
                    [0, 0, 1],
                    ])
    n_vector = np.zeros([len(coord),dim])
    for i in range(len(coord)):
        for j in range(dim):
            n_vector[i][j] = rot[j] @ R[i]
    
    # 規格化
    n_norm = np.zeros([len(coord),dim])
    for i in range(len(coord)):
        norm = np.linalg.norm(n_vector[i])
        for j in range(dim):
            n_norm[i][j] = n_vector[i][j] / norm

    return n_norm

# ねじれ角の計算する関数
def twist_cal(coord,n_bottom,n_top):
    twist = np.zeros(len(coord))
    for i in range(len(coord)):
        twist[i] = np.arcsin( np.linalg.norm( np.cross( n_bottom[i],n_top[i] ) ) ) * 180 / np.pi
    
    # Run平均
    twist_mean = np.zeros(int(len(coord)/Run))
    for i in range(int(len(coord)/Run)):
        twist_mean[i] = ( twist[i] + twist[3000+i] + twist[6000+i] ) / Run
    
    return twist_mean

# Time windowをかける関数
def window(width,coord,twist,frame):
    for i in range(len(frame)):

        imin = i - width * 0.5
        if imin < 1:
            imin = 1

        imax = i + width * 0.5
        if imax > int(len(frame)):
            imax = int(len(frame))

        twist_window[i] = sum(twist[int(imin):int(imax)]) / (imax - imin + 1)

    return twist_window

########################################
########################################
########################################

# make universe
# parmtopとtrajectoryの読み込み
mda_uni = mda.Universe(parm_name,
                       traj_name)

# Referenceの読み込み
ave_str = mda.Universe(ref_name)

# trajectoryを初期frameにセット
mda_uni.trajectory[0]

# Align trajectory
alignment = align.AlignTraj(mobile = mda_uni,
                            reference = ave_str,
                            select = fitting_area,
                            in_memory = True)
alignment.run()

# MDAnalysisのUniverseオブジェクトから座標をnumpyのarrayに取り出す
coord_from_mdaA_bottom = extract(residA_bottom)
coord_from_mdaA_top = extract(residA_top)
coord_from_mdaB_bottom = extract(residB_bottom)
coord_from_mdaB_top = extract(residB_top)

atom_number_bottom = atom_number(residA_bottom)
atom_number_top = atom_number(residA_top)

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
R_bottom = diff(coord_from_mdaA_bottom,centerA_bottom,centerB_bottom)
R_top = diff(coord_from_mdaA_top,centerA_top,centerB_top)

# 方線ベクトルの生成（回転行列の作用）
theta = np.pi * 0.5
n_vector_bottom = normal(theta,coord_from_mdaA_bottom,R_bottom)
n_vector_top = normal(theta,coord_from_mdaA_top,R_top)

# ねじれ角の計算
twist = np.zeros(int(len(coord_from_mdaA_bottom)/Run))
frame = np.zeros(int(len(coord_from_mdaA_bottom)/Run))
twist = twist_cal(coord_from_mdaA_bottom,n_vector_bottom,n_vector_top)

for i in range(int(len(coord_from_mdaA_bottom)/Run)):
    frame[i] = i + 1

# Time windowをかける 
twist_window = np.zeros(int(len(coord_from_mdaA_bottom)/Run))
twist_window = window(width,coord_from_mdaA_bottom,twist,frame)

# plot
plt.rcParams["figure.figsize"] = (8, 6)
plt.title(fig_title, fontweight = "bold", fontsize = 12)
plt.plot(frame, twist, lw = 1.25, color = twist_color)
plt.plot(frame, twist_window, lw = 1.25, color = "black")
plt.xlabel("Frame", fontweight = "bold", fontsize = 12)
plt.ylabel("Twist [Degree]", fontweight = "bold", fontsize = 12)
plt.xticks(fontweight = "bold", fontsize = 11)
plt.yticks(fontweight = "bold", fontsize = 11)
plt.xlim(0,int(len(frame)))
plt.ylim(ymin,ymax)
plt.savefig(fig_name)
