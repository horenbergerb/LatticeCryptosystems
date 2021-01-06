import numpy as np
import secrets
from math import isclose
from Utilities import hadamard_ratio, gcd, ext_gcd, rand_GLNZ, babai
from Cryptosystem import Cryptosystem

#Notes:
#the noise, r_int, seems extremely touchy. consider a better way to tune this
#the GL2N function is also touchy. values explode quickly. tuning parameters is annoying. should probably do statistical analyses on the randomness of these matrices...
#at least the hadamard ratio is really low



##############
#CRYPTOSYSTEM#
##############
#number of characters allowed
class GGHSystem():
    def __init__(self, chars=1, r_rad = 0.5):
        self.dims = (8*chars)
        
        self.reinitialize()
        
        self.r_rad = r_rad
        self.noise_generator = secrets.SystemRandom()

    def reinitialize(self):
        vecs = []
        for cur_dim in range(0,self.dims):
            cur_vec = [0.0]*self.dims
            #sets the current dimension to 1.0 or -1.0
            cur_vec[cur_dim]=float((-1)**(secrets.randbelow(2)))
            vecs.append(cur_vec)
            
        self.good_basis = np.array(vecs)
        #matrix from good basis to bad basis
        self.in_between = rand_GLNZ(self.dims)
        self.bad_basis = np.matmul(self.in_between, self.good_basis)
        
    def get_public_key(self):
        return self.bad_basis
    def get_private_key(self):
        return self.good_basis
        
    #by default it assumes the message is a string
    def encrypt(self, message, is_binary = True):
        payload = None
        if not is_binary:
            payload = self.ascii_to_bin(message)
        else:
            payload = message

        #now a list of floats
        payload = list(map(float, payload))

        #padding and checking length conditions
        if len(payload) > self.dims:
            raise Exception("Message to be encoded by encrypt() is longer than dimensions allow")
        while len(payload) < self.dims:
            payload.append(1.0)

            
        #linear combination of bad basis vectors
        scrambled = np.matmul(payload, self.bad_basis)

        #adding noise
        for bit in range(0,self.dims):
            scrambled[bit] += self.noise_generator.uniform(-self.r_rad, self.r_rad)
            
        return scrambled

    def decrypt(self, ciphertext, to_binary = True):
        good_ciphertext = babai(ciphertext, self.good_basis)
        #we round here since we expect values of 1 or 0, then convert to an integer array
        raw_plaintext = np.rint(np.matmul(good_ciphertext, np.linalg.inv(self.bad_basis))).astype(np.int)
        raw_plaintext = raw_plaintext.tolist()
        raw_plaintext = map(str, raw_plaintext)
        raw_plaintext = ''.join(raw_plaintext)

        plaintext = None
        if not to_binary:
            plaintext = self.bin_to_ascii(raw_plaintext)
        else:
            plaintext = raw_plaintext
        return plaintext
        
    
#################
#DEBUG FUNCTIONS#
#################
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
        return "FAILED TO DECRYPT"


def general_debug(chars = 5, iterations = 5):
    for cur in range(iterations):
        print("Trial {}/{}".format(cur,iterations))
        test_sys = GGHSystem(chars = chars)
        plaintext = ""
        for cur_char in range(0, chars):
            plaintext += chr(secrets.randbelow(128))
        print("Plaintext message: {}".format(plaintext))
        ciphertext = test_sys.encrypt(plaintext)
        decoded = test_sys.decrypt(ciphertext)
        print("Decoded message: {}".format(decoded))
        print("Decoded from Babai attack: {}".format(attack(ciphertext, test_sys.bad_basis)))
