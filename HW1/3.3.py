import matplotlib.pyplot as plt
import numpy as np
from skimage.io import imshow
from collections import defaultdict
import string
alphabet = string.ascii_lowercase

ALPHABET = 26

def decrypt(s, key):

    t = ""
    for c, i in zip(s, range(len(s))):
        diff = ord(key[i % len(key)]) - ord('a')
        v = (c - diff - ord('a')) % ALPHABET 
        t += chr(v + ord('a'))
    return t

def convert(l):
    s = []
    val = 0
    for i in range(len(l)):
        if i % 8 == 0:
            val = 0
        val *= 2
        val += (l[i] == 1)
        if i % 8 == 7:
            s.append(val)
    return s

# Reading Q3_img from the file
lines = [s for s in open('353055/353055-params.txt') if s[0] != '#' and len(s) > 0 and s.isspace()==False]
exec(lines[-3])

data = np.array(Q3_img, dtype=int) # N x M

red, green, blue = [(data >> (16 - 8 * i)) & 0xFF for i in range(3)]

# Using Hint 1 to convert "night" image to "day"
data = [(red & 7)<<5, ((blue>>1) & 15)<<4, ((green & 3) << 6) | ((blue & 1) << 5)]
img = np.dstack(data).astype(np.uint8) # N x M x 3
plt.plot(), imshow(img)
plt.show()
# Shifted red and blue right to the msb to make the image appear better

# "after 91 years"
bw3 = data[2][91, :]
bcopy = bw3.copy()

bw3 = (bw3>>5) & 1
cipher = convert(bw3)
print(decrypt(cipher, 'obiwan'))
bw3 = bcopy.copy()

# Brute force implementation, used it figure out 91 meant 91st row
'''
ls = []
for k in range(data[2].shape[0]):
    bw3 = data[2][k, :]
    for i in range(7):
        bcopy = bw3.copy()
    #     bw3[i] = (bw3[i] >> j) & 1
        bw3 = (bw3>>i) & 1
        cipher = convert(bw3)
        ls.append(max(cipher)-min(cipher))
        # print(decrypt(cipher, 'obiwan'))
        bw3 = bcopy.copy()
ls.sort()
print(ls[:5])
'''