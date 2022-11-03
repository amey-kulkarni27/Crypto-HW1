import matplotlib.pyplot as plt
import numpy as np
from skimage.io import imshow
from collections import defaultdict
import string
alphabet = string.ascii_lowercase
ALPHABET = 26

def decrypt(s):
    fq= defaultdict( int )
    # Alternatively, we could check if the character is in the alphabet, either lowercase or uppercase
    dont_check = ['.', ',', ' ']
    dont_check.extend([str(i) for i in range(10)])
    for w in s:
        if w not in ' ':
            fq[w] += 1
    
    most_freq = max(fq, key=fq.get)

    # 'e' is the most frequent letter in English
    shift = ord('e') - ord(most_freq)
    k = -shift % ALPHABET
    shift += ALPHABET
    shift %= ALPHABET

    t = ""
    for c in s:
        if c in dont_check:
            t += c
        else:
            add = ord(c) - ord('a')
            add += shift
            add %= ALPHABET
            add += ord('a')
            t += chr(add)

    return t, k



# Reading Q3_img from the file
lines = [s for s in open('353055/353055-params.txt') if s[0] != '#' and len(s) > 0 and s.isspace()==False]
exec(lines[-3])


Q3a_ct='ql qxhb vlr colj qeb kfdeq ql qeb phv vlr kbba ql hbbm qeb ixpq qeobb obap. qeb ixpq qtl dobbkp tfii dfsb vlr qeb qtl jlpq fjmloqxkq yirbp. cfkxiiv, tfqefk qeb yirb qeb qefoa lqebo yirb ifb ixpq, mobzbbaba yv clro dobbkp.'
Q3b_ct='diwhu91bhduvwkhwklugeoxhvnbzloouhyhdolwvvhfuhwdiwhurshqlqjwkhgrruzlwkrelzdq'

# Hint 1
m, key = decrypt(Q3a_ct)
print("Hint 1: ", m)
print(f"The key used as an integer: {key}")

print()

# Hint 2
m, key = decrypt(Q3b_ct)
print("Hint 2: ", m)
print(f"The key used as an integer: {key}")