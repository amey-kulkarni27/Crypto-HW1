L = []
N = len(L)

def enc(K, pt):
    s = len(K) - 2
    rho, rho_0, P = [None for i in range(s)], [None for i in range(s)], [None for i in range(s)] # Permutation, initial position, and turnover points
    for i in range(s):
        rho[i], rho_0[i], P[i] = K[i]
    pi = K[-2]
    sigma = K[-1]
    ct = ""
    for j in range(len(pt)):
        mj = pt[j]
        if mj not in L:
            ct += mj
        else: # This section is for the rotor movement update
            beta = set()
            rho_0_old = rho_0
            rho_0 = (rho_0 + 1) % N
            beta.add(0)
            if rho_0_old in P[0]:
                rho_1 = (rho_1 + 1) % N
                beta.add(1)
            for i in range(1, s-1):
                if i not in beta and rho_0[i] in P[i]:
                    rho_0[i] = (rho_0[i] + 1) % N
                    rho_0[i+1] = (rho_0[i+1] + 1) % N
                    beta.add(i)
        
        x = L.index(mj)

        add_tau = lambda x, y: (x + y) % N
        sub_tau = lambda x, y: (x - y + N) % N

        x = sigma[0][x]

        for i in range(s):
            x = add_tau(x, rho_0[i])
            x = rho[i][x]
            x = sub_tau(x, rho_0[i])

        x = pi[x]
        
        for i in range(s-1, -1, -1):
            x = sub_tau(x, rho_0[i])
            x = rho[i].index(x)
            x = add_tau(x, rho_0[i])
            
        x = sigma[0][x]

        ct += L[x]

        K_prime = []
        for i in range(s):
            K_prime.append([rho[i], rho_0[i], P[i]])
        K_prime.append(pi)
        K_prime.append(sigma)

        return K_prime, ct

lines = [s for s in open('353055/353055-params.txt') if s[0] != '#' and len(s) > 0 and s.isspace()==False]

for line in lines[:7]:
    exec(line)
