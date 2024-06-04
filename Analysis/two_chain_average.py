import numpy as np
import matplotlib.pyplot as plt

name = "Phe47A_Ala111B_distance"
dataA = np.genfromtxt( name + ".out", skip_header = 1)
name = "Phe47B_Ala111A_distance"
dataB = np.genfromtxt( name + ".out" )

mean = np.zeros(len(dataA))
for i in range(len(dataA)):
    mean[i] = 0.5 * ( dataA[i, 1] + dataB[i, 1] )
    if i == 0:
        print( "Frame" + "\t" + "Distance" )

    print( str(dataA[i, 0]) + "\t" + str(mean[i]) )
