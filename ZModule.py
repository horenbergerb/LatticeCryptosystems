import numpy as np
import secrets
from math import isclose
from copy import copy
from matplotlib import pyplot as plt
from Utilities import hadamard_ratio, gcd, ext_gcd, rand_GLNZ, babai

###############
#ZMODULE CLASS#
###############

class ZModule():
    #initialized with a list of n-dimensional numpy vectors
    def __init__(self, vectors):
        self.basis = vectors
        self.dim = vectors[0].shape[0]

    #renders an image of 2D ZModules
    #needs to account for commutativity and inverses, i.e. not doing redundant paths to one lattice pt
    #  so i'm just looking for unique integer summations of basis vectors
    def get_points(self, ITER=3):

        points = [np.array([0]*self.dim)]
        for cur in range(0, len(self.basis)):
            cur_points = []
            for vec in points:
                for coeff in range(1,ITER+1):
                    cur_points.append(vec + coeff*self.basis[cur])
                    cur_points.append(vec - coeff*self.basis[cur])
            points.extend(cur_points)
        return points

    #gets the r_star configurations
    #needs an efficient algorithm... add linear combos of vectors rather than checking points
    #use the get_points func
    #limit search; don't want to search too near to edges
    #i'm thinking default to searching up to ITER/2
    #maybe i should make a more general get_pts utility... that takes lower and upper bounds. then i could create simpler sets for this and voronoi calc
    def get_r_stars(self, ITER=6, limit = None, save_plot = False, plot_name = "r_stars.png"):
        if limit == None:
            limit = ITER/2
        pass
                
#################
#DEBUG FUNCTIONS#
#################    
    
def debug_zmod_get_pts():
    #tests three examples from quasicrystals textbook, pg 36
    test_mod1 = ZModule([np.array([1,0]),np.array([np.cos(2.0*3.14159/5.0),np.sin(3.14159*2.0/5.0)]), np.array([(1.0/5.0)*(3+4*np.cos(2.0*3.14159/5.0)),4*np.sqrt(2.0*3.14159/5.0)])])
    print(test_mod1.get_points())

    test_mod2 = ZModule([np.array([1,0]),np.array([np.cos(2.0*3.14159/5.0),np.sin(3.14159*2.0/5.0)]), np.array([np.cos(4.0*3.14159/5.0),np.sin(4.0*3.14159/5.0)])])
    print(test_mod2.get_points())
    
    test_mod3 = ZModule([np.array([1,0]),np.array([np.cos(2.0*3.14159/5.0),np.sin(3.14159*2.0/5.0)]), np.array([np.sqrt(2.0),np.sqrt(3.0)])])
    print(test_mod3.get_points())

debug_zmod_get_pts()
