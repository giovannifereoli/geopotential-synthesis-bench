# Table 7 emulation: CPU time of the synthesis programmes (N = 360)

Model: synthetic Kaula (N=360, seed=2009). Machine: x86_64 Darwin 25.6.0, Python 3.13.2, NumPy 2.5.1, single thread. Timing: process CPU time on 6000 evenly spaced points per scenario (600 for LEGENDRE_b), extrapolated; CLENSHAW/GRID per parallel on all points.

| Programme | GRID [s] | paper | ORBIT1 [s] | paper | ORBIT2 [h] | paper [h] | GRID / LEG_a | paper | ORBIT1 / LEG_a | paper |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| LEGENDRE_a | 100.1 | 2390 | 99.7 | 2079 | 0.86 | 17.9 | 1.000 | 1.000 | 1.00 | 1.00 |
| LEGENDRE_b | 247.2 | 9793 | 240.2 | 9975 | 2.07 | 85.9 | 2.470 | 4.097 | 2.41 | 4.80 |
| PINES | 127.9 | 5275 | 129.5 | 6306 | 1.12 | 54.3 | 1.278 | 2.207 | 1.30 | 3.03 |
| GOTTLIEB | 114.2 | - | 118.7 | - | 1.02 | - | 1.141 | - | 1.19 | - |
| METRIS | 55.0 | 3271 | 48.5 | 3418 | 0.42 | 29.4 | 0.549 | 1.369 | 0.49 | 1.64 |
| CLENSHAW | 0.4 | 8 | 131.2 | 992 | 1.13 | 8.5 | 0.004 | 0.003 | 1.32 | 0.48 |

Paper: Pentium IV 1.3 GHz, Fortran90/77 (LEGENDRE_b: Matlab + DLL), binary I/O included. ORBIT2 = ORBIT1 x 2e6 / 64,440 (paper: x 31.04). GOTTLIEB is not in the paper (added for comparison).
