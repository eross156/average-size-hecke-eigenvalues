import time, sys 

ARG = sys.argv[1]
if ARG == 'testTr':
    M_VALUES = list(range(1,17))
    N_MAX = 20
elif ARG == 'testAv2':
    M_VALUES = list(range(1,17))
    N_MAX = 20
elif ARG == 'compute':
    M_VALUES = [1,4]
    N_MAX = 150000
elif ARG == 'compute_smallNk':
    M_VALUES = [1,4]
    N_MAX = 1000
else:
    assert False, 'invalid argument: ' + ARG



################################################################
# The following implements the Eichler Selberg Trace Formula
# See "Signs of the Second Coefficients of Hecke Polynomials"
# for more details on the Eichler-Selberg trace formula.
#################################################################

def psi(N):
    ret = 1
    for pm, exp in factor(N):
        ret *= (pm+1) * pm^(exp-1)
    return ret


PSI = [0]*(N_MAX+1)
for N in range(1,N_MAX+1):
    PSI[N] = psi(N)

print('Done with psi')

def _12_A1(m, N, k):
    sqrtm = round(sqrt(m))
    if sqrtm^2 != m or gcd(sqrtm, N) != 1:
        return 0
    else:
        return (k-1) * PSI[N] * m^(k//2 - 1)
    

####################################################

def get_t_2wt(m):
    ret = []
    for t in range(4*m+1):
        if t^2 >= 4*m: 
            break
        _2wt = 1 if (t == 0) else 2
        ret.append((t, _2wt))
    return ret


def get_n(m,t):
    ret = []
    for n in range(1, 4*m - t^2 + 1):
        if (t^2 - 4*m) % (n^2) == 0 and ((t^2 - 4*m)//(n^2))%4 in [0,1]:
            ret.append(n)
    return ret




class U_seq:
    def __init__(self, t, m):
        self.t = t
        self.m = m
        self.d = [0,1]
        self.cur = 1

    def get(self,k):
        if k != self.cur: 
            if k < self.cur:
                self.d = [0,1]
                self.cur = 1
            # build up to k
            for i in range(self.cur+1, k+1):
                self.d[i%2] = self.t * self.d[(i-1)%2]  -  self.m * self.d[(i-2)%2]
                self.cur += 1
        return self.d[self.cur%2]


U = {}
for m in M_VALUES:
    for t,_ in get_t_2wt(m):
        U[(t,m)] = U_seq(t,m)



_6_h_w = {-3: 2, -4: 3, -7: 6, -8: 6, -11: 6, -12: 6, -15: 12, -16: 6, -19: 6, -20: 12, 
        -23: 18, -24: 12, -27: 6, -28: 6, -31: 18, -32: 12, -35: 12, -36: 12, -39: 24,
        -40: 12, -43: 6, -44: 18, -47: 30, -48: 12, -51: 12, -52: 12, -55: 24, -56: 24,
        -59: 18, -60: 12, -63: 24, -64: 12, -67: 6, -68: 24, -71: 42, -72: 12, -75: 12, 
        -76: 18, -79: 30, -80: 24, -83: 18, -84: 24, -87: 36, -88: 12, -91: 12, -92: 18, 
        -95: 48, -96: 24, -99: 12, -100: 12, -103: 30, -104: 36, -107: 18, -108: 18, -111: 48,
        -112: 12, -115: 12, -116: 36, -119: 60, -120: 24, -123: 12, -124: 18, -127: 30, 
        -128: 24, -131: 30, -132: 24, -135: 36, -136: 24, -139: 18, -140: 36, -143: 60, -144: 24}


def mu_sum(N,t,n,m):
    R.<x>=PolynomialRing(Integers(N))
    Nn = gcd(N,n)
    ret = 0
    poly = x^2-t*x+m 
    roots_modN = poly.roots(multiplicities=False)
    for soln_ in roots_modN:
        soln = int(soln_)
        if gcd(soln, N) != 1:
            continue
        # check if it solves the eqn for some lifting of soln
        soln_lifts = False
        for i in range(Nn):
            soln_lftd = soln + i*N
            value = soln_lftd^2 - t*soln_lftd + m
            if value % (N*Nn) == 0:
                soln_lifts = True
                break
        if soln_lifts:
            ret += 1
    return ret




MU_SUM = {}
for m in M_VALUES:
    for t,_ in get_t_2wt(m):
        for n in get_n(m, t):
            # mu_sum is a multiplicative function of N; see Cohen+Stromberg Remark 12.4.12
            mu_sm = [0]*(N_MAX+1)
            for N in range(1, N_MAX+1):
                N_fact = factor(N)
                if len(N_fact) == 1:
                    mu_sm[N] = mu_sum(N,t,n,m)
                else:
                    mu_sm[N] = product(mu_sm[pm^exp] for (pm,exp) in N_fact)
            MU_SUM[(t,n,m)] = mu_sm

print('Done with mu')


def mu(N,t,n,m):
    Nn = gcd(N,n)
    return (PSI[N] // PSI[N//Nn]) * MU_SUM[(t,n,m)][N]


def _12_A2(m,N,k):
    ret = 0 
    for t, _2wt in get_t_2wt(m):
        for n in get_n(m,t):
            ret += _2wt * U[(t,m)].get(k-1) * _6_h_w[(t^2 - 4*m)//(n^2)] * mu(N,t,n,m)
    return ret



######################################################################


# d <= sqrt(m) with weight 2 if d != sqrt(m)
def get_d_2wt(m):
    ret = []
    for d in divisors(m):
        if d^2 < m:
            ret.append((d,2))
        elif d^2 == m:
            ret.append((d,1))
    return ret



def _12_A3(m,N,k):
    ret = 0
    for d,_2_d_wt in get_d_2wt(m):
        for tau in divisors(N):
            g1 = gcd(tau, N//tau)
            g2 = gcd(N, m//d - d)
            if g2 % g1 != 0: 
                continue
            y = CRT([d,m//d], [tau,N//tau])
            if gcd(y, N) > 1:
                continue
            ret += _2_d_wt * d^(k-1) * euler_phi(gcd(tau, N//tau)) 
    return 6*ret
    

def _12_A4(m,N,k):
    if k > 2:
        return 0
    ret = 0
    for c in divisors(m):
        if gcd(N, m//c) == 1:
            ret += c
    return 12 * ret





## ErrorTerm_A ################################################


def omega(N):
    return len(factor(N))



def ErrorTerm(N):
    ret = 0.0
    ret += 14 * 2^omega(N) / PSI[N]
    ret += (1/2) * sqrt(N) * 2^omega(N) / PSI[N]
    return ret.n()  



# return kUB where    ErrorTerm(N) < (k-1)/24  for all k > kUB
# i.e. we only need to check k <= kUB
def get_kUB(N):
    ET = ErrorTerm(N)
    for k in range(2,1000,2):
        if ET < (k-1)/24:
            return k-2 
    assert False



######################################################################

def get_trace(m,N,k):
    A1 = _12_A1(m,N,k) 
    A2 = _12_A2(m,N,k) 
    A3 = _12_A3(m,N,k) 
    A4 = _12_A4(m,N,k)
    ret = A1 - A2 - A3 + A4
    # print(N,m,k, '|', A1, A2, A3, A4, '|',  ret//12)
    assert ret % 12 == 0
    return ret // 12


def get_dim(N,k):
    return get_trace(1,N,k)




def get_Av2(m,N,k):
    assert gcd(m,N)==1
    ret = 0
    for d in divisors(m):
        ret += get_trace((m//d)^2,N,k) / (m//d)^(k-1) # normalized
    return ret / get_dim(N,k)
    


def test_all_traces():
    for N in range(1, N_MAX+1):
        for m in M_VALUES:
            for k in range(2, 17, 2):
                tr = get_trace(m,N,k)
                print(N,m,k, 'Tr', tr)
                Sk = ModularForms(Gamma0(N),k).cuspidal_subspace()
                assert tr == Sk.hecke_operator(m).trace()


def test_Av2():
    for N in range(1, N_MAX+1):
        for m in [1,2,3,4]:
            if gcd(m,N) != 1:
                continue
            for k in range(2, 17, 2):
                dim = get_dim(N,k)
                if dim == 0:
                    continue
                Av2 = get_Av2(m,N,k)
                print(N,m,k, 'Av2', Av2)
                Sk = ModularForms(Gamma0(N),k).cuspidal_subspace()
                hecke_poly = [0,0,0] + Sk.hecke_operator(m).charpoly().list()
                sum_a2 = hecke_poly[-2]^2 - 2*hecke_poly[-3]
                Av2_ = (sum_a2 / m^(k-1)) / dim
                print(Av2_)
                assert Av2 == Av2_



def find_all_Av2_le_1():
    for N in range(1, N_MAX+1, 2):
        if N % 1000 == 1:
            print('checking ', N)
        kUB = get_kUB(N)
        for k in range(2, kUB+1, 2):
            dim = get_dim(N,k)
            if dim == 0:
                continue
            Av2 = get_Av2(2,N,k)
            if Av2 <= 1:
                print('Av2 <= 1:  (N,k):',N,k, '\tdim', dim, '\tAv^2', Av2)


def find_Av2_le_1_smallNk():
    for N in range(1, N_MAX+1, 2):
        for k in range(2, 100, 2):
            dim = get_dim(N,k)
            if dim == 0:
                continue
            Av2 = get_Av2(2,N,k)
            if Av2 <= 1:
                print('Av2 <= 1:  (N,k):',N,k, '\tdim', dim, '\tAv^2', Av2)


##############################################################################

if ARG == 'testTr':
    test_all_traces()
elif ARG == 'testAv2':
    test_Av2()
elif ARG == 'compute':
    find_all_Av2_le_1()
elif ARG == 'compute_smallNk':
    find_Av2_le_1_smallNk()




