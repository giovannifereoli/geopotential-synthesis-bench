# Table 7 emulation: CPU time of the synthesis programmes (N = 360)

Model: synthetic Kaula (N=360, seed=2009). Machine: x86_64 Darwin 25.6.0, Python 3.13.2, NumPy 2.5.1, single thread. Timing: process CPU time on 6000 evenly spaced points per scenario (600 for LEGENDRE_b), extrapolated; CLENSHAW/GRID per parallel on all points.

| Programme | GRID [s] | paper | ORBIT1 [s] | paper | ORBIT2 [h] | paper [h] | GRID / LEG_a | paper | ORBIT1 / LEG_a | paper |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| LEGENDRE_a | 71.8 | 2390 | 71.5 | 2079 | 0.62 | 17.9 | 1.000 | 1.000 | 1.00 | 1.00 |
| LEGENDRE_b | 182.0 | 9793 | 182.0 | 9975 | 1.57 | 85.9 | 2.535 | 4.097 | 2.54 | 4.80 |
| PINES | 99.0 | 5275 | 98.0 | 6306 | 0.84 | 54.3 | 1.379 | 2.207 | 1.37 | 3.03 |
| METRIS | 37.9 | 3271 | 36.6 | 3418 | 0.32 | 29.4 | 0.528 | 1.369 | 0.51 | 1.64 |
| CLENSHAW | 0.3 | 8 | 99.3 | 992 | 0.86 | 8.5 | 0.004 | 0.003 | 1.39 | 0.48 |

Paper: Pentium IV 1.3 GHz, Fortran90/77 (LEGENDRE_b: Matlab + DLL), binary I/O included. ORBIT2 = ORBIT1 x 2e6 / 64,440 (paper: x 31.04).
