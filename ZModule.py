import numpy as np
import secrets
from math import isclose
from copy import copy
from matplotlib import pyplot as plt

#Notes:

#Converting this to a general class for Z-Modules
#TO DO:
#add r_star and voronoi code
#create more general zmodule utilities
#think about file structure. maybe move utility funcs to another file, keep zmodule class in a separate file
#really diggin the new emacs setup; i should use emacs properly more often

###################
#UTILITY FUNCTIONS#
###################

def hadamard_ratio(matrix):
    det = abs(np.linalg.det(matrix))
    denom = 1.0
    for row in matrix:
        denom *= np.linalg.norm(row)
    return (det/denom)**(1.0/matrix.shape[0])

def gcd(a,b):
    s = 0
    old_s = 1
    r = b
    old_r = a

    while r != 0:
        quotient = old_r//r
        cur = r
        r = old_r-(quotient*r)
        old_r = cur

        cur = s
        s = old_s-(quotient*s)
        old_s = cur

    if b != 0:
        bez_t = (old_r-(old_s*a))//b
    else:
        bez_t = 0

    return old_s*a + bez_t*b
    
    
#extended euclidean algorithm
#courtesy of Wikipedia, pretty much
def ext_gcd(a, b):
    s = 0
    old_s = 1
    r = b
    old_r = a

    while r != 0:
        quotient = old_r//r
        cur = r
        r = old_r-(quotient*r)
        old_r = cur

        cur = s
        s = old_s-(quotient*s)
        old_s = cur

    if b != 0:
        bez_t = (old_r-(old_s*a))//b
    else:
        bez_t = 0

    return old_s, bez_t
        
#creates a random nxn matrix with determinant 1 or -1
#does this by using random elementary operations on the rows of the identity matrix
#inspiration from here: https://math.stackexchange.com/questions/2358484/integer-matrices-with-determinant-equal-to-1
def rand_GLNZ(n, max_int = 2**3):
    rands = secrets.SystemRandom()
    #the row vectors
    vecs = []
    for x in range(0, n):
        cur_vec = np.zeros(n)
        cur_vec[x] = 1.0
        vecs.append(cur_vec)


    #don't really like how this code is working
    for iterations in range(n*2 + secrets.randbelow(n)):
        row1 = secrets.randbelow(n)
        row2 = secrets.randbelow(n)
        if row1 == row2:
            continue
        rand_scalar = 0
        while rand_scalar == 0:
            rand_scalar = rands.randint(-max_int, max_int)
        vecs[row2] = np.add(vecs[row2], vecs[row1]*rand_scalar)

    for row in range(0, n):
        for row2 in range(0,n):
            if secrets.randbelow(2) == 1:
                temp = vecs[row]
                vecs[row] = vecs[row2]
                vecs[row2] = temp

    
    result = np.array(vecs)
    return result

#Babai's algorithm
def babai(ciphertext, basis):
    #first put the ciphertext in terms of the good basis
    good_inv = np.linalg.inv(basis)
    good_coeffs = np.matmul(ciphertext, good_inv)

    #rounding to nearest integers
    good_coeffs = np.rint(good_coeffs)

    good_vector = np.matmul(good_coeffs, basis)
    return good_vector

###############
#ZMODULE CLASS#
###############
#what features should we have?
# basic data: number of vectors, base space
# relative density of Deloine set
# compute Veronoi diagram (why?)
# compute r-stars, r-atlases
# some way to check symmetries?

#features for fun
# plot r-stars
# plot point set and Veronoi diagram for 2D case

class ZModule():
    #initialized with a list of n-dimensional numpy vectors
    def __init__(self, vectors):
        self.basis = vectors
        self.dim = vectors[0].shape[0]

    #renders an image of 2D ZModules
    #needs to account for commutativity and inverses, i.e. not doing redundant paths to one lattice pt
    #  so i'm just looking for unique integer summations of basis vectors
    def get_points(self, ITER=3, save_plot = False, plot_name = "lattice.png"):
        if self.dim != 2:
            raise Exception("Plot is only compatible with ZModules in R^2")

        points = [np.array([0]*self.dim)]
        for cur in range(0, len(self.basis)):
            cur_points = []
            for vec in points:
                for coeff in range(1,ITER+1):
                    cur_points.append(vec + coeff*self.basis[cur])
                    cur_points.append(vec - coeff*self.basis[cur])
            points.extend(cur_points)
        if save_plot:
            fig = plt.figure()
            pts = []
            for x in zip(*points):
                pts.append(x)
            plt.scatter(pts[0],pts[1], marker=".", s=.5*plt.rcParams['lines.markersize']**2)
            plt.savefig(plot_name)
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
        

##############
#CRYPTOSYSTEM#
##############
#number of characters allowed
class GGHSystem():
    def __init__(self, chars=1, r_rad = 0.5):
        self.dims = (8*chars)
        vecs = []
        for cur_dim in range(0,self.dims):
            cur_vec = [0.0]*self.dims
            #sets the current dimension to 1.0 or -1.0
            cur_vec[cur_dim]=float((-1)**(secrets.randbelow(2)))
            vecs.append(cur_vec)
            
        self.good_basis = np.array(vecs)
        print("Hadamard ratio of good basis: {}".format(hadamard_ratio(self.good_basis)))
        #matrix from good basis to bad basis
        self.in_between = rand_GLNZ(self.dims)
        self.bad_basis = np.matmul(self.in_between, self.good_basis)
        print("Hadamard ratio of bad basis: {}".format(hadamard_ratio(self.bad_basis)))
        
        self.r_rad = r_rad
        self.noise_generator = secrets.SystemRandom()
        
    #by default it assumes the message is a string
    def encrypt(self, message):
        payload = None
        
        #message is encoded as a binary vector
        #string form of binary first
        payload = ''.join(bin(ord(x)).replace('b','') for x in message)
        #now a list of floats
        payload = list(map(float, payload))
        print("Raw payload:")
        print(payload)

        #padding and checking length conditions
        if len(payload) > self.dims:
            raise Exception("Message to be encoded by encrypt() is longer than dimensions allow")
        while len(payload) < self.dims:
            payload.append(1.0)
        
        print("Padded payload:")
        print(payload)
            
        #linear combination of bad basis vectors
        scrambled = np.matmul(payload, self.bad_basis)

        #adding noise
        for bit in range(0,self.dims):
            scrambled[bit] += self.noise_generator.uniform(-self.r_rad, self.r_rad)
            
        return scrambled

    def decrypt(self, ciphertext):
        good_ciphertext = babai(ciphertext, self.good_basis)
        #we round here since we expect values of 1 or 0, then convert to an integer array
        raw_plaintext = np.rint(np.matmul(good_ciphertext, np.linalg.inv(self.bad_basis))).astype(np.int)
        raw_plaintext = raw_plaintext.tolist()
        raw_plaintext = map(str, raw_plaintext)
        raw_plaintext = ''.join(raw_plaintext)

        
        #converting the binary back into a string
        #not super graceful
        plaintext = ""
        for cur in range(0,(len(raw_plaintext)//8)):
            plaintext += chr(int(raw_plaintext[8*cur: 8*(cur+1)], 2))
        
        return plaintext
        
    
#################
#DEBUG FUNCTIONS#
#################

def debug_ext_gcd(quantity=100, print_out=False, max_int=2**10):
    for x in range(0, quantity):
        a = secrets.randbelow(max_int)
        #to prevent a==0
        a = a+1
        
        b = 0
        while (b == 0) or (gcd(b,a) != 1):
            b = secrets.randbelow(max_int)
        d,c = ext_gcd(a,b)

        if print_out:
            print("a: {}, b:{}, c:{}, d:{}".format(a,b,c,d))
            print("ad+bc={}".format(a*d+b*c))
        if a*d+b*c != 1:
            raise Exception("ext_gcd failed to calculate coefficients")

def debug_rand_GLNZ(n, quantity=100, print_out=False):
    for x in range(0, quantity):
        cur = rand_GLNZ(n)
        det = np.linalg.det(cur)

        if print_out:
            print(cur)
            print("Det: {}".format(det))

        #these rounding errors are gnarly
        if not isclose(abs(det), 1.0, rel_tol = 1e-3):
            raise Exception("rand_GL2N generated matrix has determinant != 1 or -1")

def debug_babai():
    #solution should be 53159, 81818 (from textbook)
    sol = babai(np.array([53172.0, 81743.0]), np.array([[137.0, 312.0],[215.0,-187.0]]))
    print("Should have: 53159, 81818")
    print("Calculated: {}".format(sol))

#simple "attack" using babai's algorithm directly
def attack(ciphertext, bad_basis):
    try:
        good_ciphertext = babai(ciphertext, bad_basis)

        raw_plaintext = np.rint(np.matmul(good_ciphertext, np.linalg.inv(bad_basis))).astype(np.int)
        raw_plaintext = raw_plaintext.tolist()
        raw_plaintext = map(str, raw_plaintext)
        raw_plaintext = ''.join(raw_plaintext)

        #converting the binary back into a string
        #not super graceful
        plaintext = ""
        for cur in range(0,(len(raw_plaintext)//8)):
            plaintext += chr(int(raw_plaintext[8*cur: 8*(cur+1)]))


        return plaintext
    except ValueError:
        print("Something went wrong decoding; failed to unscramble message")
    
def debug_zmod_plot():
    #tests three examples from quasicrystals textbook, pg 36
    test_mod1 = ZModule([np.array([1,0]),np.array([np.cos(2.0*3.14159/5.0),np.sin(3.14159*2.0/5.0)]), np.array([(1.0/5.0)*(3+4*np.cos(2.0*3.14159/5.0)),4*np.sqrt(2.0*3.14159/5.0)])])
    test_mod1.get_points(save_plot=True, plot_name = "debuglattice1.png")

    test_mod2 = ZModule([np.array([1,0]),np.array([np.cos(2.0*3.14159/5.0),np.sin(3.14159*2.0/5.0)]), np.array([np.cos(4.0*3.14159/5.0),np.sin(4.0*3.14159/5.0)])])
    test_mod2.get_points(save_plot=True, plot_name = "debuglattice2.png")
    
    test_mod3 = ZModule([np.array([1,0]),np.array([np.cos(2.0*3.14159/5.0),np.sin(3.14159*2.0/5.0)]), np.array([np.sqrt(2.0),np.sqrt(3.0)])])
    test_mod3.get_points(save_plot=True, plot_name = "debuglattice3.png")

debug_zmod_plot()
