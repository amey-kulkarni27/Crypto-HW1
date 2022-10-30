import matplotlib.pyplot as plt
import numpy as np
from skimage.io import imshow
from collections import defaultdict
import string
alphabet = string.ascii_lowercase
ALPHABET = 26

def decrypt(s):
    fq= defaultdict( int )
    dont_check = ['.', ',', ' ']
    for w in s:
        if w not in ' ':
            fq[w] += 1
    
    most_freq = max(fq, key=fq.get)

    shift = ord('e') - ord(most_freq)
    print(shift)
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

    print(t)
    
    return t



# Reading Q3_img from the file
lines = [s for s in open('353055/353055-params.txt') if s[0] != '#' and len(s) > 0 and s.isspace()==False]
exec(lines[-3])

data = np.array(Q3_img, dtype=int) # N x M
red, blue, green = [(data >> (16 - 8 * i)) & 0xFF for i in range(3)]
data = [(red & 7) << 5, green & 0xc0, ((green & 3) << 6) | (blue >> 5) & 1 | (green & 0x4c)]
# data = [(red & 0x07) << 5,  ((blue >> 1) & 0x0F) << 4, ((green & 0x03)<<6) | (((blue>>5) & 0x01)<<5)]
# data = [(data >> (16 - 8 * i)) & 0xFF for i in range(3)]

# green_l2 = data[1] & 3
# green_4 = (data[1] & 30)
# blue_third = data[2] & 32
# data[0] = data[0] & 7
# data[1] = data[1]
# data[2] = (data[2] & 63) | (green_l2 << 6)
# data[2] = (data[2] & 254) | (blue_third == 32)
# data[2] = (data[2] & 225) | (green_4)


# bw = data[2][67, :]
# pt = ""
# counter_key = 0
# key = "obiwan"
# pt = ""
# for i in bw:
#     pt+= alphabet[(i + ord(key[counter_key]))%26]
#     counter_key += 1
#     if counter_key == 6:
#         counter_key = 0

bw0 = ((data[2] >> 5) & 1).astype(bool)
# print(bw0 % 26)
# fq= defaultdict( int )
# for w in bw0:
#     fq[w] += 1

# most_freq = max(fq, key=fq.get)
# print(most_freq)

# shift = ord('e') - ord(most_freq)

img = np.dstack(data).astype(np.uint8) # N x M x 3
# plt.plot(), imshow(img)
# plt.plot(), imshow(img)
plt.plot(), imshow(img)
plt.show()

Q3a_ct='ql qxhb vlr colj qeb kfdeq ql qeb phv vlr kbba ql hbbm qeb ixpq qeobb obap. qeb ixpq qtl dobbkp tfii dfsb vlr qeb qtl jlpq fjmloqxkq yirbp. cfkxiiv, tfqefk qeb yirb qeb qefoa lqebo yirb ifb ixpq, mobzbbaba yv clro dobbkp.'
Q3b_ct='diwhu91bhduvwkhwklugeoxhvnbzloouhyhdolwvvhfuhwdiwhurshqlqjwkhgrruzlwkrelzdq'

decrypt(Q3a_ct)
decrypt(Q3b_ct)

