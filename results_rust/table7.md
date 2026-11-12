# Table 7 emulation: CPU time of the synthesis programmes (N = 360)

Model: synthetic Kaula (N=360, seed=2009). Machine: Darwin 25.6.0, rustc 1.87.0, aarch64-apple-darwin, opt-level 3 + LTO, scalar kernels (no fast-math), single thread. Timing: CPU time of the synthesis loop inside the binary, every point, best of 2.

| Programme | GRID [s] | paper | ORBIT1 [s] | paper | ORBIT2 [h] | paper [h] | GRID / LEG_a | paper | ORBIT1 / LEG_a | paper |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| LEGENDRE_a | 11.1 | 2390 | 11.2 | 2079 | 0.10 | 17.9 | 1.000 | 1.000 | 1.00 | 1.00 |
| LEGENDRE_b | 32.2 | 9793 | 32.1 | 9975 | 0.28 | 85.9 | 2.908 | 4.097 | 2.86 | 4.80 |
| PINES | 26.4 | 5275 | 27.0 | 6306 | 0.23 | 54.3 | 2.383 | 2.207 | 2.41 | 3.03 |
| GOTTLIEB | 25.8 | - | 25.8 | - | 0.22 | - | 2.328 | - | 2.30 | - |
| METRIS | 21.3 | 3271 | 21.5 | 3418 | 0.19 | 29.4 | 1.927 | 1.369 | 1.92 | 1.64 |
| CLENSHAW | 0.2 | 8 | 38.0 | 992 | 0.33 | 8.5 | 0.016 | 0.003 | 3.39 | 0.48 |

Paper: Pentium IV 1.3 GHz, Fortran90/77 (LEGENDRE_b: Matlab + DLL), binary I/O included. ORBIT2 = ORBIT1 x 2e6 / 64,440 (paper: x 31.04). GOTTLIEB is not in the paper (added for comparison).
