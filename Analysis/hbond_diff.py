"""
This scripts perform hydrogen analysis with mdtraj and create the heatmap for vizualization
"""

# Information
# File: hbond.py
# Author: Fukuda Akiya
# Date: 2025-03-21
# Description: Producing the hydrogen anlysis

import sys
import os
import mdtraj as md
import numpy as np
np.set_printoptions(threshold=np.inf)
from collections import Counter
import matplotlib.pyplot as plt

trajin1 = "/home/users/dzn/Kenkyu/OaPAC/Analysis/compare/holo/Angle/cluster/8qfg/average/closedA_closedB/closedA_closedB.nc"
parm1 = "/home/users/dzn/Kenkyu/OaPAC/Analysis/8QFG/parmtop/updated_for_contact_wofmn.prmtop"

trajin2 = "/home/users/dzn/Kenkyu/OaPAC/Analysis/compare/holo/Angle/cluster/8qfg/average/openA_openB/openA_openB.nc" 
parm2 = "/home/users/dzn/Kenkyu/OaPAC/Analysis/8QFG/parmtop/updated_for_contact_wofmn.prmtop"

select1 = 'resSeq 1 to 700'

def hbond(trajin, parm, select):
    # load trajectory and topology
    traj = md.load(trajin, top=parm)
    top = traj.topology
    selected_atoms = traj.topology.select(select)
    selected_residues = set(top.atom(idx).residue.index for idx in selected_atoms)
    print(f"traj information:{traj}")
    
    # get the information of the hydrogen bond
    hbond_count = Counter()
    for i in range(len(traj)):
        hbonds = md.baker_hubbard(traj[i], freq=0.0, distance_cutoff=0.3, angle_cutoff=135)
        
        seen_res_pairs = set()
        for j in range(len(hbonds)):
            donor_residue = top.atom(hbonds[j][0]).residue
            acceptor_residue = top.atom(hbonds[j][2]).residue
            res_pair = tuple([donor_residue, acceptor_residue])  
            if donor_residue.index in selected_residues and acceptor_residue.index in selected_residues:
                if res_pair not in seen_res_pairs:
                    hbond_count[res_pair] += 1
                    seen_res_pairs.add(res_pair)

    for res_pair in hbond_count:
            hbond_count[res_pair] /= len(traj)

    residues = list(set([residue for pair in hbond_count.keys() for residue in pair]))
    # mapping of residue name and index
    residue_to_idx = {residue.index: residue.index for residue in residues}

    # Initialize the matrix
    contact_matrix = np.zeros((700, 700))
    for i in range(700):
        for j in range(700):
            if i == j:
                contact_matrix[i, j] = 1
            else:
                contact_matrix[i, j] = 0

    # Substitute occupancy from hbond_count to matrix
    for res_pair, count in hbond_count.items():
        if count > 0.00:
            i = residue_to_idx[res_pair[0].index]
            j = residue_to_idx[res_pair[1].index]
            contact_matrix[i, j] = count
            contact_matrix[j, i] = count

    return contact_matrix, hbond_count

def diff(hbond_count1, hbond_count2,contact_matrix1, contact_matrix2):
    all_res_pairs = set(hbond_count1.keys()).union(set(hbond_count2.keys()))

    hbond_diff = Counter()
    for res_pair in all_res_pairs:
            hbond_diff[res_pair] = hbond_count1.get(res_pair, 0) - hbond_count2.get(res_pair, 0)
                    
    diff_matrix = np.subtract(contact_matrix1, contact_matrix2)

    return hbond_diff,diff_matrix

def plot_heatmap(selecter,matrix):
    plt.figure(figsize=(8, 8), facecolor='white')
    if selecter == "contact":
        plt.imshow(matrix, cmap="Purples", origin="lower", vmin=0, vmax=1)
    elif selecter == "diff":
        plt.imshow(matrix, cmap="seismic", origin="lower", vmin=-1, vmax=1)
    else:
        exit(1)

    plt.colorbar()
    plt.xlabel("Residue Number", fontsize=18, fontweight='bold')
    plt.ylabel("Residue Number", fontsize=18, fontweight='bold')
    plt.xticks(fontsize=13, fontweight='bold')
    plt.yticks(fontsize=13, fontweight='bold')
    plt.title("Hydrogen Bond Contact Map", fontsize=18, fontweight='bold')
    plt.xlim(0,143)
    plt.ylim(0,143)
    plt.savefig("diff_switch_wtdark_AvsA.png", bbox_inches='tight', pad_inches = 0, dpi=300)

def file_write(selecter, file_name, hbond_count):
    if selecter == "diff":
        with open(file_name, "w") as fo:
            for res_pair, count in hbond_count.most_common():
                if count >= 0.1 or count <= -0.1:
                    fo.write(f"{res_pair}: {count:.2f}\n")
    elif selecter == "contact":
        with open(file_name, "w") as fo:
            for res_pair, count in hbond_count.most_common():
                if count >= 0.1:
                    fo.write(f"{res_pair}: {count:.2f}\n")
    else:
        exit(1)
     
def main():
    contact_matrix1, hbond_count1  = hbond(trajin1, parm1, select1)
    contact_matrix2, hbond_count2  = hbond(trajin2, parm2, select1)
    hbond_diff, diff_matrix = diff(hbond_count1, hbond_count2, contact_matrix1, contact_matrix2)
    #file_write("diff", "diff.dat", hbond_diff)
    plot_heatmap("diff",diff_matrix)

if __name__ == "__main__":
    main()
