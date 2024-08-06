# Average size of Hecke Eigenvalues

In the paper [On the Average Size of the Eigenvalues of the Hecke Operators](https://arxiv.org/abs/2407.19076), we showed that `Av_m(N,k)` tends to `sqrt(sigma1(m)/m)` as `N+k --> infinity`
This repository contains code to compute the finitely many pairs `(N,k)` such that `Av_m(N,k) <= 1`

- `compute_Avle1.sage` is the code to compute the finitely many pairs `(N,k)` such that `Av_m(N,k) <= 1`
- `output.txt` is the output of `compute_Avle1.sage`
- `check_conj.sage` is code to compute the sizes of Fourier coefficients of newforms. This is to verify that the conjecture stated in the paper looks reasonable. 


