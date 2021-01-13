import numpy as np
import secrets
from math import isclose
from copy import copy
from matplotlib import pyplot as plt
import sympy
from Utilities import hadamard_ratio, gcd, ext_gcd, rand_GLNZ, babai


#################
#SAMPLE ZMODULES#
#################
SAMPLE_ZMODULES = dict()
SAMPLE_ZMODULES['textbook1'] = [np.array([1,0]),np.array([np.cos(2.0*3.14159/5.0),np.sin(3.14159*2.0/5.0)]), np.array([(1.0/5.0)*(3+4*np.cos(2.0*3.14159/5.0)),4*np.sqrt(2.0*3.14159/5.0)])]
SAMPLE_ZMODULES['textbook2'] = [np.array([1,0]),np.array([np.cos(2.0*3.14159/5.0),np.sin(3.14159*2.0/5.0)]), np.array([np.cos(4.0*3.14159/5.0),np.sin(4.0*3.14159/5.0)])]
SAMPLE_ZMODULES['textbook3'] = [np.array([1,0]),np.array([np.cos(2.0*3.14159/5.0),np.sin(3.14159*2.0/5.0)]), np.array([np.sqrt(2.0),np.sqrt(3.0)])]
                                                                                                    

###############
#ZMODULE CLASS#
###############

class ZModule():
    #initialized with a list of n-dimensional numpy vectors
    def __init__(self, vectors):
        self.basis = vectors
        self.dim = vectors[0].shape[0]

    #gets the orbit of the lattice
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

    #gets one possible lattice which projects to this Z-module
    def get_projected_lattice(self):
        '''
        1)Let L={b_1,...,b_n,b_{n+1},...,b_k} be vectors generating Z-module in n-space
          Let L be organized s.t. B={b_1,...,b_n} spans the space
        2)Embed all the vectors in k-space
        3)Find orthonormal basis A=a_1,...a_k s.t. a_{n+1},...,a_k is orthogonal complement of B
        4)"Lift" b_{n+1},...,b_k to lin. ind. vectors in k-space by adding multiples of a_j to b_j
        '''
        
        #1)
        matrix = sympy.Matrix(np.array(self.basis))
        _, indices = matrix.T.rref()
        #2)
        extended_basis = sympy.Matrix(np.array([np.append(self.basis[x], [0]*(len(self.basis)-len(indices)), 0) for x in range(len(self.basis))]))
        print("Original basis")
        print(matrix)
        print("Extended basis")
        print(extended_basis)
        #3)
        print("Bang")
        B, B_indices = extended_basis.T[indices,:].rref()
        print(B)
        print(B_indices)
        print(set(x for x in range(0, B.shape[1])))
        print(set(B_indices))
        A_indices = tuple(set(x for x in range(0, B.shape[1])).difference(set(B_indices)))
        print(A_indices)
        A = B[A_indices,:]
        print(A)
        A = sympy.Matrix(sympy.BlockMatrix([[A],[sympy.eye(len(self.basis)-A.shape[0])]]))
        print(A)
        
        
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
