
import itertools
from copy import deepcopy

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

inds_uk = [i for i in range(len(pi)) if pi[i] == -1]
print(inds_uk)
# print(inds)
# print(Q1c_ct[:20])
# print(Q1c_pt_view[:20])
best = 100
ans = ""
fin_st = None
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
    # for pf in range(len(pfill)):
    #     pi_uk[pfill[pf]] = pfill[pf]
    assert(len(pi_uk) == len(set(pi_uk)))
    assert(min(pi_uk) == 0)
    for i in init_posns:
        cur_s = 0
        succ = True
        inds[0] = i
        R = [(all_perms[k], inds[k], turnover_pts[k]) for k in range(5)]
        K = (R, pi_uk, (fpt, nfpt))
        st, pt = enc(K, Q1c_ct, Q1c_L)
        for k in range(len(pt)):
            if Q1c_pt_view[k] != -1 and pt[k] != Q1c_pt_view[k]:
                cur_s += 1
        if cur_s < best:
            best = cur_s
            ans = pt
            fin_st = st
        # if succ:
        #     print(i)
        #     print(''.join(chr(i) for i in pt))
ans_str = (''.join(chr(i) for i in ans))
# print(fin_st)
print(ans_str)



Q1c_pt_str_c='If you keep bad-mouthing someone, you might feel like as if you are on a higher position than that person. But in actual fact, its completely wrong to do that! (Yoshioka Futaba, Ao Haru Ride)'
Q1c_sk_c=[[[[3, 0, 43, 36, 21, 48, 33, 29, 51, 8, 49, 20, 26, 18, 23, 31, 39, 4, 37, 25, 40, 6, 17, 50, 16, 11, 47, 38, 14, 10, 12, 7, 9, 46, 30, 32, 5, 2, 19, 28, 35, 22, 42, 15, 41, 24, 13, 44, 34, 45, 27, 1], 20, []], [[39, 17, 36, 50, 16, 31, 19, 1, 29, 11, 49, 37, 38, 8, 27, 12, 0, 21, 20, 32, 34, 15, 5, 40, 42, 47, 10, 46, 30, 51, 35, 2, 6, 4, 43, 45, 48, 25, 7, 18, 33, 13, 28, 41, 22, 24, 23, 14, 3, 44, 9, 26], 8, [11, 29, 50]], [[31, 1, 29, 37, 20, 21, 35, 50, 33, 0, 6, 44, 51, 11, 48, 19, 25, 41, 7, 10, 36, 14, 22, 32, 46, 39, 47, 45, 9, 18, 43, 26, 49, 17, 27, 42, 40, 23, 15, 8, 12, 2, 24, 4, 28, 3, 16, 38, 34, 30, 5, 13], 25, [1, 5]], [[25, 36, 1, 24, 19, 12, 14, 18, 8, 44, 45, 50, 40, 34, 9, 28, 0, 4, 47, 15, 10, 13, 46, 39, 5, 17, 51, 20, 41, 43, 3, 35, 38, 21, 33, 27, 26, 11, 42, 30, 37, 7, 29, 48, 32, 49, 31, 16, 23, 6, 2, 22], 30, [7, 29]], [[23, 6, 8, 34, 12, 27, 26, 11, 30, 20, 7, 4, 16, 43, 15, 51, 13, 33, 18, 32, 25, 41, 0, 24, 44, 29, 47, 17, 1, 36, 38, 2, 48, 39, 14, 45, 35, 46, 3, 42, 19, 22, 21, 37, 40, 49, 9, 10, 50, 5, 31, 28], 3, []]], [49, 45, 47, 41, 16, 38, 42, 18, 31, 28, 34, 12, 11, 46, 33, 36, 4, 43, 7, 51, 30, 27, 26, 50, 37, 29, 22, 21, 9, 25, 20, 8, 48, 14, 10, 40, 15, 24, 5, 44, 35, 3, 6, 17, 39, 1, 13, 2, 32, 0, 23, 19], [[0, 1, 2, 3, 18, 45, 44, 30, 35, 9, 11, 10, 16, 13, 14, 15, 12, 17, 4, 19, 41, 21, 40, 23, 24, 25, 26, 27, 28, 29, 7, 31, 34, 36, 32, 8, 33, 37, 38, 39, 22, 20, 42, 43, 6, 5, 46, 47, 48, 49, 50, 51], 30]]
Q1c_pt_c='SWYgeW91IGtlZXAgYmFkLW1vdXRoaW5nIHNvbWVvbmUsIHlvdSBtaWdodCBmZWVsIGxpa2UgYXMgaWYgeW91IGFyZSBvbiBhIGhpZ2hlciBwb3NpdGlvbiB0aGFuIHRoYXQgcGVyc29uLiBCdXQgaW4gYWN0dWFsIGZhY3QsIGl0cyBjb21wbGV0ZWx5IHdyb25nIHRvIGRvIHRoYXQhIChZb3NoaW9rYSBGdXRhYmEsIEFvIEhhcnUgUmlkZSk='
# print(len(fin_st[0][0]), len(Q1c_sk_c[0][0]))
# print(len(fin_st))
# print(len(Q1c_sk_c[0][1]))
# for i in range(52):
#     if fin_st[0][0][i] != Q1c_sk_c[0][0][0][i]:
#         print(i, fin_st[0][0][i], Q1c_sk_c[0][0][i])
ctr = 0
# for i in range(len(ans_str)):
#     if(ans_str[i] != Q1c_pt_str_c[i]):
#         ctr += 1
#         print(i, ans_str[i], Q1c_pt_str_c[i], chr(Q1c_ct[i]))


for i in range(52):
    if fin_st[1][i] != Q1c_sk_c[1][i]:
        print(i, fin_st[1][i], Q1c_sk_c[1][i])
# print(fin_st[2][1] == Q1c_sk_c[2][1])
