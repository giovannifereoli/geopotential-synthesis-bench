# Table 7 emulation: CPU time of the synthesis programmes (N = 360)

Model: synthetic Kaula (N=360, seed=2009). Machine: Darwin 25.6.0, rustc 1.87.0, aarch64-apple-darwin, opt-level 3 + LTO, scalar kernels (no fast-math), single thread. Timing: CPU time of the synthesis loop inside the binary, every point, best of 2.

| Programme | GRID [s] | paper | ORBIT1 [s] | paper | ORBIT2 [h] | paper [h] | GRID / LEG_a | paper | ORBIT1 / LEG_a | paper |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| LEGENDRE_a | 10.3 | 2390 | 10.3 | 2079 | 0.09 | 17.9 | 1.000 | 1.000 | 1.00 | 1.00 |
| LEGENDRE_b | 32.0 | 9793 | 32.0 | 9975 | 0.28 | 85.9 | 3.109 | 4.097 | 3.11 | 4.80 |
| PINES | 26.4 | 5275 | 26.5 | 6306 | 0.23 | 54.3 | 2.565 | 2.207 | 2.58 | 3.03 |
| METRIS | 20.1 | 3271 | 20.3 | 3418 | 0.18 | 29.4 | 1.958 | 1.369 | 1.98 | 1.64 |
| CLENSHAW | 0.2 | 8 | 38.0 | 992 | 0.33 | 8.5 | 0.017 | 0.003 | 3.70 | 0.48 |

Paper: Pentium IV 1.3 GHz, Fortran90/77 (LEGENDRE_b: Matlab + DLL), binary I/O included. ORBIT2 = ORBIT1 x 2e6 / 64,440 (paper: x 31.04).
