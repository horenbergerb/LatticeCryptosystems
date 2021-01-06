import numpy as np
import secrets
from math import isclose
from copy import copy
from matplotlib import pyplot as plt
from Utilities import hadamard_ratio, gcd, ext_gcd, rand_GLNZ, babai
from ZModule import ZModule

###############
#ZMODULE CLASS#
###############

class ZModule2D(ZModule):
    #initialized with a list of n-dimensional numpy vectors
    def __init__(self, vectors):
        self.basis = vectors
        self.dim = vectors[0].shape[0]
        if self.dim != 2:
            raise Exception("ZModule2D initialized with vectors not in R^2")

    def plot_points(self, ITER=3, plot_name = "lattice.png"):

        points = self.get_points(ITER)

        fig = plt.figure()
        pts = []
        for x in zip(*points):
            pts.append(x)
        plt.scatter(pts[0],pts[1], marker=".", s=.5*plt.rcParams['lines.markersize']**2)
        plt.savefig(plot_name)

#################
#DEBUG FUNCTIONS#
#################    
    
def debug_zmod_plot():
    #tests three examples from quasicrystals textbook, pg 36
    test_mod1 = ZModule2D([np.array([1,0]),np.array([np.cos(2.0*3.14159/5.0),np.sin(3.14159*2.0/5.0)]), np.array([(1.0/5.0)*(3+4*np.cos(2.0*3.14159/5.0)),4*np.sqrt(2.0*3.14159/5.0)])])
    test_mod1.plot_points(plot_name = "debuglattice1.png")

    test_mod2 = ZModule2D([np.array([1,0]),np.array([np.cos(2.0*3.14159/5.0),np.sin(3.14159*2.0/5.0)]), np.array([np.cos(4.0*3.14159/5.0),np.sin(4.0*3.14159/5.0)])])
    test_mod2.plot_points(plot_name = "debuglattice2.png")
    
    test_mod3 = ZModule2D([np.array([1,0]),np.array([np.cos(2.0*3.14159/5.0),np.sin(3.14159*2.0/5.0)]), np.array([np.sqrt(2.0),np.sqrt(3.0)])])
    test_mod3.plot_points(plot_name = "debuglattice3.png")

debug_zmod_plot()

    
