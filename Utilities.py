import numpy as np
import secrets
from math import isclose
from copy import copy
from matplotlib import pyplot as plt

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

