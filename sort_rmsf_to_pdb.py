import numpy as np
import pandas as pd
import itertools as it

resnum = 701
datanum = 2804
rmsf = [0] * (resnum+1)
rmsf_split = [] * 2804
rmsf_save = 0.0
rmsf_save2 = [0] * (resnum + 1)

with open("4YUS_holo/RMSF_data.dat","r") as frmsf:
     rmsf = frmsf.readlines()

for i in range(0,resnum):
    rmsf_save = rmsf[i]
    rmsf_save2[i] = rmsf_save.split()[1]
    for j in range(0,4):
        rmsf_split.append(rmsf_save2[i])

pdb_data = []
pdb_HET = [0] * datanum
pdb_atom_name = [0] * datanum
pdb_atom_index = [0] * datanum
pdb_resname = [0] * datanum
pdb_chain = [0] * datanum
pdb_resid = [0] * datanum
pdb_x = [0] * datanum
pdb_y = [0] * datanum
pdb_z = [0] * datanum
pdb_occupancy = [0] * datanum
pdb_atom = [0] * datanum
pdb_rmsf = []
pdb_data_save = 0.0

with open("4YUS_holo/4YUS_dark_run_strip_allstep_reference.pdb","r") as fpdb:
     next(fpdb) # title read skip
     pdb_data = fpdb.readlines()

for i in range(datanum):
    pdb_data_save = pdb_data[i]
    pdb_HET[i] = pdb_data_save.split()[0] 
    pdb_atom_index[i] = pdb_data_save.split()[1] 
    pdb_atom_name[i] = pdb_data_save.split()[2] 
    pdb_resname[i] = pdb_data_save.split()[3] 
    pdb_chain[i] = pdb_data_save.split()[4] 
    pdb_resid[i] = pdb_data_save.split()[5] 
    pdb_x[i] = pdb_data_save.split()[6] 
    pdb_y[i] = pdb_data_save.split()[7] 
    pdb_z[i] = pdb_data_save.split()[8] 
    pdb_occupancy[i] = pdb_data_save.split()[9]
    pdb_atom[i] = pdb_data_save.split()[11]

    if len(pdb_atom_index[i]) == 1:
        if len(pdb_atom_name[i]) == 2:
            if len(pdb_resid[i]) == 1:
                print(pdb_HET[i] + "      " + pdb_atom_index[i] + "  " + pdb_atom_name[i] + "  " + 
                      pdb_resname[i] + " " + pdb_chain[i] + "   " + pdb_resid[i] + "      " +
                      pdb_x[i] + "  " + pdb_y[i] + "  " + pdb_z[i] + "  " +
                      pdb_occupancy[i] + "  " + "{:.6}".format(rmsf_split[i]) + "           " + pdb_atom[i])
            elif len(pdb_resid[i]) == 2:
                print(pdb_HET[i] + "      " + pdb_atom_index[i] + "  " + pdb_atom_name[i] + "  " + 
                      pdb_resname[i] + " " + pdb_chain[i] + "  " + pdb_resid[i] + "      " +
                      pdb_x[i] + "  " + pdb_y[i] + "  " + pdb_z[i] + "  " +
                      pdb_occupancy[i] + "  " + "{:.6}".format(rmsf_split[i]) + "           " + pdb_atom[i])
            elif len(pdb_resid[i]) == 3:
                print(pdb_HET[i] + "      " + pdb_atom_index[i] + "  " + pdb_atom_name[i] + "  " + 
                      pdb_resname[i] + " " + pdb_chain[i] + " " + pdb_resid[i] + "      " +
                      pdb_x[i] + "  " + pdb_y[i] + "  " + pdb_z[i] + "  " +
                      pdb_occupancy[i] + "  " + "{:.6}".format(rmsf_split[i]) + "           " + pdb_atom[i])
        elif len(pdb_atom_name[i]) == 1:
            if len(pdb_resid[i]) == 1:
                print(pdb_HET[i] + "      " + pdb_atom_index[i] + "  " + pdb_atom_name[i] + "   " + 
                      pdb_resname[i] + " " + pdb_chain[i] + "   " + pdb_resid[i] + "      " +
                      pdb_x[i] + "  " + pdb_y[i] + "  " + pdb_z[i] + "  " +
                      pdb_occupancy[i] + "  " + "{:.6}".format(rmsf_split[i]) + "           " + pdb_atom[i])
            elif len(pdb_resid[i]) == 2:
                print(pdb_HET[i] + "      " + pdb_atom_index[i] + "  " + pdb_atom_name[i] + "   " + 
                      pdb_resname[i] + " " + pdb_chain[i] + "  " + pdb_resid[i] + "      " +
                      pdb_x[i] + "  " + pdb_y[i] + "  " + pdb_z[i] + "  " +
                      pdb_occupancy[i] + "  " + "{:.6}".format(rmsf_split[i]) + "           " + pdb_atom[i])
            elif len(pdb_resid[i]) == 3:
                print(pdb_HET[i] + "      " + pdb_atom_index[i] + "  " + pdb_atom_name[i] + "   " + 
                      pdb_resname[i] + " " + pdb_chain[i] + " " + pdb_resid[i] + "      " +
                      pdb_x[i] + "  " + pdb_y[i] + "  " + pdb_z[i] + "  " +
                      pdb_occupancy[i] + "  " + "{:.6}".format(rmsf_split[i]) + "           " + pdb_atom[i])
    elif len(pdb_atom_index[i]) == 2:
        if len(pdb_atom_name[i]) == 2:
            if len(pdb_resid[i]) == 1:
                print(pdb_HET[i] + "     " + pdb_atom_index[i] + "  " + pdb_atom_name[i] + "  " + 
                      pdb_resname[i] + " " + pdb_chain[i] + "   " + pdb_resid[i] + "      " +
                      pdb_x[i] + "  " + pdb_y[i] + "  " + pdb_z[i] + "  " +
                      pdb_occupancy[i] + "  " + "{:.6}".format(rmsf_split[i]) + "           " + pdb_atom[i])
            elif len(pdb_resid[i]) == 2:
                print(pdb_HET[i] + "     " + pdb_atom_index[i] + "  " + pdb_atom_name[i] + "  " + 
                      pdb_resname[i] + " " + pdb_chain[i] + "  " + pdb_resid[i] + "      " +
                      pdb_x[i] + "  " + pdb_y[i] + "  " + pdb_z[i] + "  " +
                      pdb_occupancy[i] + "  " + "{:.6}".format(rmsf_split[i]) + "           " + pdb_atom[i])
            elif len(pdb_resid[i]) == 3:
                print(pdb_HET[i] + "     " + pdb_atom_index[i] + "  " + pdb_atom_name[i] + "  " + 
                      pdb_resname[i] + " " + pdb_chain[i] + " " + pdb_resid[i] + "      " +
                      pdb_x[i] + "  " + pdb_y[i] + "  " + pdb_z[i] + "  " +
                      pdb_occupancy[i] + "  " + "{:.6}".format(rmsf_split[i]) + "           " + pdb_atom[i])
        elif len(pdb_atom_name[i]) == 1:
            if len(pdb_resid[i]) == 1:
                print(pdb_HET[i] + "     " + pdb_atom_index[i] + "  " + pdb_atom_name[i] + "   " + 
                      pdb_resname[i] + " " + pdb_chain[i] + "   " + pdb_resid[i] + "      " +
                      pdb_x[i] + "  " + pdb_y[i] + "  " + pdb_z[i] + "  " +
                      pdb_occupancy[i] + "  " + "{:.6}".format(rmsf_split[i]) + "           " + pdb_atom[i])
            elif len(pdb_resid[i]) == 2:
                print(pdb_HET[i] + "     " + pdb_atom_index[i] + "  " + pdb_atom_name[i] + "   " + 
                      pdb_resname[i] + " " + pdb_chain[i] + "  " + pdb_resid[i] + "      " +
                      pdb_x[i] + "  " + pdb_y[i] + "  " + pdb_z[i] + "  " +
                      pdb_occupancy[i] + "  " + "{:.6}".format(rmsf_split[i]) + "           " + pdb_atom[i])
            elif len(pdb_resid[i]) == 3:
                print(pdb_HET[i] + "     " + pdb_atom_index[i] + "  " + pdb_atom_name[i] + "   " + 
                      pdb_resname[i] + " " + pdb_chain[i] + " " + pdb_resid[i] + "      " +
                      pdb_x[i] + "  " + pdb_y[i] + "  " + pdb_z[i] + "  " +
                      pdb_occupancy[i] + "  " + "{:.6}".format(rmsf_split[i]) + "           " + pdb_atom[i])
    elif len(pdb_atom_index[i]) == 3:
        if len(pdb_atom_name[i]) == 2:
            if len(pdb_resid[i]) == 1:
                print(pdb_HET[i] + "    " + pdb_atom_index[i] + "  " + pdb_atom_name[i] + "  " + 
                      pdb_resname[i] + " " + pdb_chain[i] + "   " + pdb_resid[i] + "      " +
                      pdb_x[i] + "  " + pdb_y[i] + "  " + pdb_z[i] + "  " +
                      pdb_occupancy[i] + "  " + "{:.6}".format(rmsf_split[i]) + "           " + pdb_atom[i])
            elif len(pdb_resid[i]) == 2:
                print(pdb_HET[i] + "    " + pdb_atom_index[i] + "  " + pdb_atom_name[i] + "  " + 
                      pdb_resname[i] + " " + pdb_chain[i] + "  " + pdb_resid[i] + "      " +
                      pdb_x[i] + "  " + pdb_y[i] + "  " + pdb_z[i] + "  " +
                      pdb_occupancy[i] + "  " + "{:.6}".format(rmsf_split[i]) + "           " + pdb_atom[i])
            elif len(pdb_resid[i]) == 3:
                print(pdb_HET[i] + "    " + pdb_atom_index[i] + "  " + pdb_atom_name[i] + "  " + 
                      pdb_resname[i] + " " + pdb_chain[i] + " " + pdb_resid[i] + "      " +
                      pdb_x[i] + "  " + pdb_y[i] + "  " + pdb_z[i] + "  " +
                      pdb_occupancy[i] + "  " + "{:.6}".format(rmsf_split[i]) + "           " + pdb_atom[i])
        elif len(pdb_atom_name[i]) == 1:
            if len(pdb_resid[i]) == 1:
                print(pdb_HET[i] + "    " + pdb_atom_index[i] + "  " + pdb_atom_name[i] + "   " + 
                      pdb_resname[i] + " " + pdb_chain[i] + "   " + pdb_resid[i] + "      " +
                      pdb_x[i] + "  " + pdb_y[i] + "  " + pdb_z[i] + "  " +
                      pdb_occupancy[i] + "  " + "{:.6}".format(rmsf_split[i]) + "           " + pdb_atom[i])
            elif len(pdb_resid[i]) == 2:
                print(pdb_HET[i] + "    " + pdb_atom_index[i] + "  " + pdb_atom_name[i] + "   " + 
                      pdb_resname[i] + " " + pdb_chain[i] + "  " + pdb_resid[i] + "      " +
                      pdb_x[i] + "  " + pdb_y[i] + "  " + pdb_z[i] + "  " +
                      pdb_occupancy[i] + "  " + "{:.6}".format(rmsf_split[i]) + "           " + pdb_atom[i])
            elif len(pdb_resid[i]) == 3:
                print(pdb_HET[i] + "    " + pdb_atom_index[i] + "  " + pdb_atom_name[i] + "   " + 
                      pdb_resname[i] + " " + pdb_chain[i] + " " + pdb_resid[i] + "      " +
                      pdb_x[i] + "  " + pdb_y[i] + "  " + pdb_z[i] + "  " +
                      pdb_occupancy[i] + "  " + "{:.6}".format(rmsf_split[i]) + "           " + pdb_atom[i])
    elif len(pdb_atom_index[i]) == 4:
        if len(pdb_atom_name[i]) == 2:
            if len(pdb_resid[i]) == 1:
                print(pdb_HET[i] + "   " + pdb_atom_index[i] + "  " + pdb_atom_name[i] + "  " + 
                      pdb_resname[i] + " " + pdb_chain[i] + "   " + pdb_resid[i] + "      " +
                      pdb_x[i] + "  " + pdb_y[i] + "  " + pdb_z[i] + "  " +
                      pdb_occupancy[i] + "  " + "{:.6}".format(rmsf_split[i]) + "           " + pdb_atom[i])
            elif len(pdb_resid[i]) == 2:
                print(pdb_HET[i] + "   " + pdb_atom_index[i] + "  " + pdb_atom_name[i] + "  " + 
                      pdb_resname[i] + " " + pdb_chain[i] + "  " + pdb_resid[i] + "      " +
                      pdb_x[i] + "  " + pdb_y[i] + "  " + pdb_z[i] + "  " +
                      pdb_occupancy[i] + "  " + "{:.6}".format(rmsf_split[i]) + "           " + pdb_atom[i])
            elif len(pdb_resid[i]) == 3:
                print(pdb_HET[i] + "   " + pdb_atom_index[i] + "  " + pdb_atom_name[i] + "  " + 
                      pdb_resname[i] + " " + pdb_chain[i] + " " + pdb_resid[i] + "      " +
                      pdb_x[i] + "  " + pdb_y[i] + "  " + pdb_z[i] + "  " +
                      pdb_occupancy[i] + "  " + "{:.6}".format(rmsf_split[i]) + "           " + pdb_atom[i])
        elif len(pdb_atom_name[i]) == 1:
            if len(pdb_resid[i]) == 1:
                print(pdb_HET[i] + "   " + pdb_atom_index[i] + "  " + pdb_atom_name[i] + "   " + 
                      pdb_resname[i] + " " + pdb_chain[i] + "   " + pdb_resid[i] + "      " +
                      pdb_x[i] + "  " + pdb_y[i] + "  " + pdb_z[i] + "  " +
                      pdb_occupancy[i] + "  " + "{:.6}".format(rmsf_split[i]) + "           " + pdb_atom[i])
            elif len(pdb_resid[i]) == 2:
                print(pdb_HET[i] + "   " + pdb_atom_index[i] + "  " + pdb_atom_name[i] + "   " + 
                      pdb_resname[i] + " " + pdb_chain[i] + "  " + pdb_resid[i] + "      " +
                      pdb_x[i] + "  " + pdb_y[i] + "  " + pdb_z[i] + "  " +
                      pdb_occupancy[i] + "  " + "{:.6}".format(rmsf_split[i]) + "           " + pdb_atom[i])
            elif len(pdb_resid[i]) == 3:
                print(pdb_HET[i] + "   " + pdb_atom_index[i] + "  " + pdb_atom_name[i] + "   " + 
                      pdb_resname[i] + " " + pdb_chain[i] + " " + pdb_resid[i] + "      " +
                      pdb_x[i] + "  " + pdb_y[i] + "  " + pdb_z[i] + "  " +
                      pdb_occupancy[i] + "  " + "{:.6}".format(rmsf_split[i]) + "           " + pdb_atom[i])


