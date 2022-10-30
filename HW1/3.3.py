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
        # v %= ALPHABET
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

data = [(data >> (16 - 8 * i)) & 0xFF for i in range(3)]
# red, blue, green = [(data >> (16 - 8 * i)) & 0xFF for i in range(3)]

# data = [(red & 7) << 5, (green >> 2) & 3, ((green & 3) << 6) | ((blue & 1) << 5)]
# data = [(red & 7) << 5, green & 0xc0, ((green & 3) << 6) | (blue >> 5) & 1 | (green & 0x4c)]

# data = [(data >> (16 - 8 * i)) & 0xFF for i in range(3)]
# green_l2 = data[1] & 3
# green_4 = (data[1] & 30)
# blue_third = data[2] & 32
# data[0] = data[0] & 7
# data[1] = data[1] & 252
# data[2] = (data[2] & 63) | (green_l2 << 6)
# data[2] = (data[2] & 254) | (blue_third == 1)
# data[2] = (data[2] & 225) | (green_4)

bw3 = data[2][91, :]
for j in range(1):
    bcopy = bw3.copy()
    # for i in range(bw3.shape[0]):
    #     bw3[i] = (bw3[i] >> j) & 1
    bw3 = (bw3>>j) & 1
    cipher = convert(bw3)
    print(max(cipher)-min(cipher))
    print(decrypt(cipher, 'obiwan'))
    bw3 = bcopy.copy()
# ctr = 0
# l = []
# for r in range(bw3.shape[0]):
#     x = convert(bw3[r])
#     l.append(max(x) - min(x))
# l.sort()
# print(l[:10])

# Brute force implementation
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