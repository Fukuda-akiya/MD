import MDAnalysis as mda
from sklearn.decomposition import PCA
import numpy as np

np.set_printoptions(threshold=np.inf)

# トポロジーとトラジェクトリーのファイルパス
topology = '/home/users/dzn/Kenkyu/OaPAC/Analysis/8QFG/parmtop/updated_bb.prmtop'
trajectory = '/home/users/dzn/Kenkyu/OaPAC/Analysis/8QFG/traj/8QFG_run01-05.nc'

#########################################
###############Function##################
#########################################

# Universeから座標をndarrayに取り出す関数
def extract(name):
    atom = u.select_atoms(name)
    coord_from_mda = np.array([atom.positions for frame in u.trajectory])

    return coord_from_mda

# 二つのヘリックス間の重心ベクトルを計算する関数
def center(data1,data2):
    vector = np.zeros(data1.shape)
    for i in range(data1.shape[0]):
        for j in range(data1.shape[1]):
            vector[i, j] = ( data2[i,j] - data1[i,j] ) * 0.5

    return vector

# 鈍角を鋭角に修正する関数
def fix_angle(data):
    angle = np.zeros(len(data))
    for i in range(len(data)):
        if data[i] > 90.0:
           angle[i] = 180.0 - data[i]
        else:
           angle[i] = data[i]

    return angle

########################################
########################################
########################################

# Universeを作成
u = mda.Universe(topology, trajectory)

# 交差するα-ヘリックスの選択（例：残基範囲を指定）
helix1 = 'resid 103:126 and name CA'
helix2 = 'resid 453:476 and name CA'
a4a = 'resid 127:136 and name CA'
a4b = 'resid 477:486 and name CA'


# 座標をNumPy配列に変換
positions_helix1 = extract(helix1)
positions_helix2 = extract(helix2)
positions_a4a = extract(a4a)
positions_a4b = extract(a4b)

# 交差する2つのヘリックスの重心の計算
center_vec = center(positions_helix1, positions_helix2)

# PCAの結果を保存するリスト
principal_vectors_helix1 = []
#principal_vectors_helix2 = []
principal_vectors_a4a = []
principal_vectors_a4b = []

for i in range(len(u.trajectory)):
        # a3 ヘリックスの主軸
        pca_helix1 = PCA(n_components=1)
        pca_helix1.fit(center_vec[i])
        principal_vector_helix1 = pca_helix1.components_[0]
        principal_vectors_helix1.append(principal_vector_helix1)
        # a4a ヘリックスの主軸
        pca_a4a = PCA(n_components=1)
        pca_a4a.fit(positions_a4a[i])
        principal_vector_a4a = pca_a4a.components_[0]
        principal_vectors_a4a.append(principal_vector_a4a)
        # a4b ヘリックスの主軸
        pca_a4b = PCA(n_components=1)
        pca_a4b.fit(positions_a4b[i])
        principal_vector_a4b = pca_a4b.components_[0]
        principal_vectors_a4b.append(principal_vector_a4b)

# NumPy配列に変換
principal_vectors_helix1 = np.array(principal_vectors_helix1)
#principal_vectors_helix2 = np.array(principal_vectors_helix2)
principal_vectors_a4a = np.array(principal_vectors_a4a)
principal_vectors_a4b = np.array(principal_vectors_a4b)

angle_degrees = np.zeros(len(u.trajectory))
angle_degrees_fix = np.zeros(len(u.trajectory))

#ベクトル間の角度を計算
for i in range(len(u.trajectory)):
    dot_product = np.dot(principal_vectors_helix1[i], principal_vectors_a4b[i])
    norm_helix1 = np.linalg.norm(principal_vectors_helix1[i])
    norm_helix2 = np.linalg.norm(principal_vectors_a4b[i])
    cos_theta = dot_product / (norm_helix1 * norm_helix2)
    angle = np.arccos(cos_theta)  # ラジアンでの角度
    angle_degrees[i] = np.degrees(angle)  # 度単位の角度
    
    angle_degrees_fix = fix_angle(angle_degrees)
    print(str(i) + "\t" + f"{angle_degrees_fix[i]:.4f}")
