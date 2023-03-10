import numpy as np
import matplotlib.pyplot as plt
import os, re

for f in os.listdir("liste_score/"):
    
    if f.startswith('score'):
        # read in data
        data = np.genfromtxt(f) 
        #close previous figure, if one exists
        plt.close()
        #create new figure and do plotting
        fig = plt.figure()
        ax = plt.subplot(111)
        ax.plot(data)
        #save figure
        plt.savefig("plot_score/"+f[:-4]+".png")