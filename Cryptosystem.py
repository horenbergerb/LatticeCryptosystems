from abc import ABC, abstractmethod
##############
#CRYPTOSYSTEM#
##############
#number of characters allowed
class Cryptosystem(ABC):

    def __init__(self, dims):
        self.dims = dims
    
    def reinitialize(self):
                raise Exception("Implementation of abstract class Cryptosystem did not implement reinitialize(self)")
    
    def get_public_key(self):
        raise Exception("Implementation of abstract class Cryptosystem did not implement get_public_key(self)")

    def get_private_key(self):
                raise Exception("Implementation of abstract class Cryptosystem did not implement get_private_key(self)")

    #if is_binary, assumes input is string of binary digits, otherwise assumes ascii string
    def encrypt(self, message, is_binary=True):
        raise Exception("Implementation of abstract class Cryptosystem did not implement encrypt(self, message)")

    #if to_binary, converts to string of binary, otherwise tries to convert to a string of ascii
    def decrypt(self, ciphertext, to_binary=True):
                raise Exception("Implementation of abstract class Cryptosystem did not implement decrypt(self, ciphertext)")
            
    #turns a string of ascii characters into a binary string
    def ascii_to_bin(self, message):
        payload = None
        
        #message is encoded as a binary vector
        #string form of binary first
        payload = ""
        for x in message:
            cur_bin = bin(ord(x)).replace('b','')
            cur_bin = "0"*(8-len(cur_bin)) + cur_bin
            if len(cur_bin) > 8:
                raise Exception("Plaintext has invalid characters (must be 0-128 ascii)")
            payload += cur_bin

    #converts binary string into ascii character string
    def bin_to_ascii(self, message):
        plaintext = ""
        for cur in range(0,(len(message)//8)):
            plaintext += chr(int(message[8*cur: 8*(cur+1)], 2))
        
        return plaintext
