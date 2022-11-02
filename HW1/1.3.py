
import itertools
from copy import deepcopy
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

        # x = sigma[x]
        x = sigma[0][x]

        for i in range(s):
            x = add_tau(x, rho_0[i])
            x = rho[i][x]
            x = sub_tau(x, rho_0[i])

        x = pi[x]
        
        for i in range(s-1, -1, -1):
            x = add_tau(x, rho_0[i])
            x = rho[i].index(x)
            x = sub_tau(x, rho_0[i])
            
        # x = sigma[x]
        x = sigma[0][x]
        # print(x)

        ct.append(L[x])

    K_prime = [[]]
    for i in range(s):
        K_prime[0].append([rho[i], rho_0[i], P[i]])
    K_prime.append(pi)
    K_prime.append(sigma)

    return K_prime, ct


lines = [s for s in open('/Users/amey2704/Documents/Fall22/Crypto/Cryptography-and-Security-COM-401/HW1/353055/353055-params.txt') if s[0] != '#' and len(s) > 0 and s.isspace()==False]

for line in lines[14:21]:
    exec(line)

print("Number of characters: ", len(Q1c_L))
print(f"Length of cipher text: {len(Q1c_ct)}. Number of plaintext characters known: {len(Q1c_pt_view) - Q1c_pt_view.count(-1)}")
print("Number of characters for which we have hints: ", (Q1c_trace.keys()))
print(f"Number of rotors: {len(Q1c_trace[1][0])}")
# print(f"{(Q1c_trace[1][2])}")
# print(Q1c_ct)
# print(Q1c_pt_view) 90 -> 73

def find_minus(l):
    return l.count(-1)

def fixed_pts(l):
    x = [l[i] - i for i in range(len(l))]
    return x.count(0)

def find_reflectors(Q1c_trace):
    reflectors = [-1 for i in range(52)]
    for i in Q1c_trace.keys():
        cur_pi = Q1c_trace[i][1]
        for key, value in cur_pi.items():
            reflectors[key] = value
            reflectors[value] = key
    # for i in Q1c_trace.keys():
    #     Q1c_trace[i][1] = reflectors
    return reflectors

def find_nfpt(Q1c_trace):
    for i in Q1c_trace.keys():
        sigmas = Q1c_trace[i][2]
        if(sigmas[1] != -1):
            return sigmas[1]

def find_fpt(Q1c_trace):
    sigma = [-1 for i in range(52)]
    for i in Q1c_trace.keys():
        sigmas = Q1c_trace[i][2]
        for key, value in sigmas[0].items():
            sigma[key] = value
            sigma[value] = key
    # for i in Q1c_trace.keys():
    #     Q1c_trace[i][2][0] = sigma
    return sigma

def fill_fixed(fpt):
    for i in range(len(fpt)):
        if fpt[i] == -1:
            fpt[i] = i
    # for i in Q1c_trace.keys():
    #     Q1c_trace[i][2][1] = fpt
    return fpt

def find_turnover_points(Q1c_trace):
    turnover_points = [set() for i in range(5)]
    for i in Q1c_trace.keys():
        rotors = Q1c_trace[i][0]
        for r_no in range(len(rotors)):
            cur_turnover = rotors[r_no][2]
            for t in cur_turnover:
                if t != -1:
                    turnover_points[r_no].add(t)
    for r_no in range(len(rotors)):
        turnover_points[r_no] = list(turnover_points[r_no])
    return turnover_points

def find_permutations(Q1c_trace):
    perms = [[-1 for i in range(52)] for i in range(5)]
    for i in Q1c_trace.keys():
        rotors = Q1c_trace[i][0]
        for r_no in range(len(rotors)):
            cur_permutation = rotors[r_no][0]
            for p, fp in cur_permutation.items():
                perms[r_no][p] = fp
    return perms

def mapToArray(mp):
    return [mp[i] if i in mp else -1 for i in range(N)]

def find_initial_positions(Q1c_trace, turnovers, Q1c_ct, Q1c_Z):
    N = 52
    m_keys = Q1c_trace.keys()
    # init_values = Q1c_trace[1][0] # init value of rotors

    counter = 0
    res = []
    # Q1c_trace[4][0][0][1]
    while counter < 52:
        rotors = [None for i in range(5)]
        init_values = Q1c_trace[1][0][:] # init value of rotors
        for i in range(5):
            rotors[i] = [[], init_values[i][1], turnovers[i]]
            
        rotors[0][1] = counter

        flag = 1
        
        for idx, character in enumerate(Q1c_ct):
            if character in Q1c_Z:
                continue
            
            if idx+1 in m_keys:
                # print(idx+1, character)
                for k in range(1, 5):
                    r_k_0 =  Q1c_trace[idx+1][0][k][1]
                    # print(r_k_0,rotors[k][1] )
                    if rotors[k][1] != r_k_0:
                        flag = 0
                        break
                        
            if flag == 0:
                break
            beta = set()
            # permutation init_pos turnover
            p_1_0_old = rotors[0][1]
            rotors[0][1] = (rotors[0][1] + 1) % N
            if p_1_0_old in rotors[0][2]:
                rotors[1][1] = (rotors[1][1] + 1) % N
                beta.add(1)
            for i in range(1, len(rotors)-1):
                if i not in beta and rotors[i][1] in rotors[i][2]:
                    rotors[i][1] = (rotors[i][1] + 1) % N
                    rotors[i+1][1] = (rotors[i+1][1] + 1) % N
                    beta.add(i+1)

        if flag:
            res.append(counter)
            
        counter += 1
        
    return res

pi = find_reflectors(Q1c_trace)
print(Q1c_trace[2][1])
# print(find_minus(pi))
nfpt = find_nfpt(Q1c_trace)
fpt = find_fpt(Q1c_trace)
# print(fixed_pts(fpt), find_minus(fpt), nfpt)
fpt = fill_fixed(fpt)
# print(fixed_pts(fpt), find_minus(fpt), nfpt)
turnover_pts = find_turnover_points(Q1c_trace)


inds = [Q1c_trace[1][0][i][1] for i in range(5)]

init_posns = find_initial_positions(Q1c_trace, turnover_pts, Q1c_ct, Q1c_Z)
all_perms = find_permutations(Q1c_trace)

# print(Q1c_ct[:20])
# print(Q1c_pt_view[:20])

# print(Q1c_trace[1][0][0][0])

# inds_uk = [i for i in range(len(pi)) if pi[i] == -1]
# print(inds_uk)
# # print(inds)
# # print(Q1c_ct[:20])
# # print(Q1c_pt_view[:20])
# best = 100
# ans = ""
# fin_st = None
# for i in range(3):
#     pfill = inds_uk[:]
#     pi_uk = deepcopy(pi)
#     pi_uk[pfill[i+1]] = pfill[0]
#     pi_uk[pfill[0]] = pfill[i+1]
#     others = []
#     for j in range(1, 4):
#         if j != i+1:
#             others.append(j)
#     pi_uk[pfill[others[0]]] = pfill[others[1]]
#     pi_uk[pfill[others[1]]] = pfill[others[0]]
#     # for pf in range(len(pfill)):
#     #     pi_uk[pfill[pf]] = pfill[pf]
#     assert(len(pi_uk) == len(set(pi_uk)))
#     assert(min(pi_uk) == 0)
#     for i in init_posns:
#         inds[0] = i
#         R = [[all_perms[k], inds[k], turnover_pts[k]] for k in range(5)]
#         K = [R, pi_uk, [fpt, nfpt]]
#         K_save = deepcopy(K)
#         st, pt = enc(K, Q1c_ct, Q1c_L)
#         flag = 1
#         for k in range(len(pt)):
#             if Q1c_pt_view[k] != -1 and pt[k] != Q1c_pt_view[k]:
#                 flag = 0
#                 break
#         if flag:
#             ans =  (''.join(chr(i) for i in pt))
#             break

inds_uk = []
for idd, t in enumerate(pi):
    if t == -1:
        inds_uk.append(idd)
        
for i in range(3):
    pfill = inds_uk[:]
    pi_uk = deepcopy(pi)
    pi_uk[pfill[i+1]] = pfill[0]
    pi_uk[pfill[0]] = pfill[i+1]
    others = []
    for j in range(1, 4):
        if j != i+1:
            others.append(j)
    pi_uk[pfill[others[0]]] = pfill[others[1]]
    pi_uk[pfill[others[1]]] = pfill[others[0]]
    last_pi = deepcopy(pi_uk)
    for s in init_posns:
        rotors = [None for i in range(5)]
        for i in range(5):
            rotors[i] = [all_perms[i], inds[i], turnover_pts[i]]
        rotors[0][1] = s
        
        K = (rotors, pi_uk, [fpt, nfpt])
        K_save = deepcopy(K)
        K_prime, res = enc(K, Q1c_ct, Q1c_L)
        flag = 1
        for idx, ch in enumerate(res):
            if Q1c_pt_view[idx] == -1:
                continue
            if Q1c_pt_view[idx] != ch:
                flag = 0
                break
                
            
        if flag:
            pt = ""
            for i in res:
                pt += chr(i)
            print(pt)
            break
            
    if flag:
        break

for i in range(5):
    K_save[0][i][2].sort()

def checksum(*args, sep=';'):
    data = sep.join(map(str, args)).encode()
    return hashlib.new('md5', data=data).hexdigest()

data = base64.b64encode(pt.encode()).decode() # type: str
print(data)
assert checksum(data) == Q1c_pt_checksum

assert checksum([K_save[0], K_save[1], K_save[2]]) == Q1c_sk_checksum
Q1c_sk = [K_save[0], K_save[1], K_save[2]]
print(Q1c_sk)