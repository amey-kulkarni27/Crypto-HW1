import hashlib
import base64

def enc(K, pt, L):
    N = len(L)
    s = len(K[0])
    rho, rho_0, P = [None for i in range(s)], [None for i in range(s)], [None for i in range(s)] # Permutation, initial position, and turnover points
    for i in range(s):
        # print(len(K[i]))
        rho[i], rho_0[i], P[i] = K[0][i]
    pi = K[1]
    sigma = K[2]
    ct = []
    for j in range(len(pt)):
    # for j in range(1):
        mj = pt[j]
        if mj not in L:
            ct.append(mj)
            continue
        else: # This section is for the rotor movement update
            beta = set()
            rho_0_0_old = rho_0[0]
            rho_0[0] = (rho_0[0] + 1) % N
            beta.add(0)
            if rho_0_0_old in P[0]:
                rho_0[1] = (rho_0[1] + 1) % N
                beta.add(1)
            for i in range(1, s-1):
                if i not in beta and rho_0[i] in P[i]:
                    rho_0[i] = (rho_0[i] + 1) % N
                    rho_0[i+1] = (rho_0[i+1] + 1) % N
                    beta.add(i+1)
        
        x = L.index(mj)

        add_tau = lambda x, y: (x + y) % N
        sub_tau = lambda x, y: (x - y + N) % N

        x = sigma[x]
        # x = sigma[0][x]

        for i in range(s):
            x = add_tau(x, rho_0[i])
            x = rho[i][x]
            x = sub_tau(x, rho_0[i])

        x = pi[x]
        
        for i in range(s-1, -1, -1):
            x = add_tau(x, rho_0[i])
            x = rho[i].index(x)
            x = sub_tau(x, rho_0[i])
            
        x = sigma[x]
        # x = sigma[0][x]
        # print(x)

        ct.append(L[x])

    K_prime = []
    for i in range(s):
        K_prime.append([rho[i], rho_0[i], P[i]])
    K_prime.append(pi)
    K_prime.append(sigma)

    return K_prime, ct

lines = [s for s in open('353055/353055-params.txt') if s[0] != '#' and len(s) > 0 and s.isspace()==False]

for line in lines[:14]:
    exec(line)

# K = (Q1a_R, Q1a_pi, Q1a_sigma)
# # print(len(Q1a_R[4]))
# print(enc(K, Q1a_pt, Q1a_L)[1] == Q1a_ct)




def checksum(*args, sep=';'):
    data = sep.join(map(str, args)).encode()
    return hashlib.new('md5', data=data).hexdigest()


K = (Q1b_R, Q1b_pi, Q1b_sigma)
Q1b_pt = enc(K, Q1b_ct, Q1b_L)[1]

mess = ""
for i in Q1b_pt:
    mess += chr(i)
print(mess)
data = base64.b64encode(mess.encode()).decode() # type: str
print(data)
assert checksum(data) == Q1b_pt_checksum
print(checksum(data) == Q1b_pt_checksum)