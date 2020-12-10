import numpy as np
from fractions import gcd
import secrets
from math import isclose

#Notes:
#the noise, r_int, seems extremely touchy. i'm concerned if the bad basis values are too big that the noise is getting washed out. since the "attack" with babai's algorithm on the bad basis fails when r_int=0.5 and succeeds when r_int=0.0, it seems like the noise is working

#why did i make a random generation process for SL(2,Z)???
#i can remove that, but i might hold on to the code.

#update: the GL2N function is very touchy. values explode quickly. tuning parameters is annoying
#not super happy with the generation. there's got to be a better way to tune it
#should probably do statistical analyses on the randomness
#at least the hadamard ratio is really low

###################
#UTILITY FUNCTIONS#
###################

def hadamard_ratio(matrix):
    det = abs(np.linalg.det(matrix))
    denom = 1.0
    for row in matrix:
        denom *= np.linalg.norm(row)
    return (det/denom)**(1.0/matrix.shape[0])

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
        
#random matrix from sl2z
#(a,b)
#(c,d)
#max_int is the max int that will be used for a,b
#not sure this is the best way to pick a matrix
def rand_SL2Z(max_int = 2**10):
    #might return 0 but never returns max_int
    a = secrets.randbelow(max_int)
    #solved
    a = a+1

    b = 0
    while (b == 0) or (gcd(b,a) != 1):
        b = secrets.randbelow(max_int)
    d,c = ext_gcd(a,b)
    c = -c

    return np.array([[a,b],[c,d]])

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

def debug_rand_SL2Z(quantity=100, print_out=False):
    for x in range(0, quantity):
        cur = rand_SL2Z()
        det = np.linalg.det(cur)

        if print_out:
            print(cur)
            print("Det: {}".format(det))

        if not isclose(det, 1.0):
            raise Exception("rand_SL2Z generated matrix has determinant != 1")

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
    
#debug_rand_GLNZ(4, print_out=True)
#debug_ext_gcd(print_out=True)
#debug_babai()

test_sys = GGHSystem(chars = 5)
ciphertext = test_sys.encrypt("hello")

decoded = test_sys.decrypt(ciphertext)
print("DECODED MESSAGE:")
print(decoded)

print("DECODED MESSAGE FROM SIMULATED ATTACK:")
attack_decoded = attack(ciphertext, test_sys.bad_basis)
print(attack_decoded)
