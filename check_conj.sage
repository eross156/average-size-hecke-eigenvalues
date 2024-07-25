# CM: 9,4   9,10
# Supersingular: 95,4
# weird: 5,4


M = 400
N_ub = 100
K_ub = 20


for N in range(1,N_ub+1):
    for k in range(4,K_ub+1,2):
        print(N,k,M)
        NF = ModularForms(N,k).cuspidal_subspace().newforms(names='a')
        for ff in NF:
            print('computing', ff)
            print('#'*100)
            q_exp = ff.qexp(M)
            # compute running minimum
            K = ff.hecke_eigenvalue_field()
            phi_s= K.embeddings(CC)
            for phi in phi_s:
                print('computing', phi)
                cur_min = 10.0
                for m in range(1,M):
                    val = abs(phi(  q_exp[m]  )) 
                    if val == 0.0:
                        m_fact = factor(m)
                        if len(m_fact) == 1: # prime power
                            print('vanishing', m_fact) 
                        continue
                    val_nm = ( val / m^((k-1)/2) ).n()
                    if val_nm < cur_min:
                        cur_min = val_nm
                        print(m, '\t', val_nm, '\t', val_nm*m)
    




