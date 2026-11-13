# Spherical-harmonic gravity synthesis: accuracy and CPU time of nine methods

This repo compares nine ways to evaluate a spherical-harmonic gravity model at scattered body-fixed points, such as the points of an orbit. Each method computes two outputs:

- **Acceleration:** the gravitational acceleration (3 components).
- **Acceleration + gravity gradient:** the acceleration plus the gravity-gradient tensor (6 more components). Orbit determination needs this second output for the variational equations.

All nine methods are Rust kernels that follow one contract. Each is checked against a double-double reference and timed with the same protocol. The main timings use a build with the compiler's auto-vectorization turned off, so that every kernel runs as plain scalar code; a second build with it turned on shows how much the compiler alone changes the ranking (see [Auto-vectorization](#auto-vectorization)). Five of the methods come from Fantino & Casotto (2009), [Fantino_Casotto_2009_JGeod_Methods_of_harmonic_synthesis.pdf](Fantino_Casotto_2009_JGeod_Methods_of_harmonic_synthesis.pdf). The other four are Gottlieb (1993), Lear (1986), Holmes & Featherstone (2002), the method used by Orekit, and Cunningham (1970) in the form given by Montenbruck & Gill (2000), the method used by Scarabaeus. Lear's method gives only the acceleration; the gravity gradient for LEAR is our own extension (see [LEAR gravity gradient](#lear-gravity-gradient-our-extension)). For HOLMES we use the recursion of Holmes & Featherstone with our own pole-safe derivatives (see [HOLMES](#holmes-holmesfeatherstone)). `CUNNINGHAM` uses the same solid harmonics as `METRIS`, but applies the derivative operators at every point instead of to the coefficients (see [CUNNINGHAM](#cunningham-montenbruck--gill)).

## Main findings

These findings are for the build with auto-vectorization turned off (see the next point). The time ratios below are medians over the three point sets (as in the figures); the times in µs are at the ORBIT1 points (as in the tables).

- **All methods are accurate.** Where they do not overflow, every method has a max. relative error between 2.2 × 10⁻¹³ and 3.1 × 10⁻¹², for both outputs, at every N from 10 to 2000 and on all three point sets, including points 10⁻⁸ degrees from the poles. The error grows slowly with N. For the acceleration, `CUNNINGHAM` has the smallest error from N = 100 on, 1.3–1.6 times smaller than that of the next best method. Otherwise no method is consistently more accurate than the others.
- **Three methods overflow at high degree in plain double precision.** `PINES` and `GOTTLIEB` overflow at N = 1500 and 2000, and `CLENSHAW` at N = 2000, at high-latitude points. See [why](#why-pines-gottlieb-and-clenshaw-overflow).
- **The speed ranking depends on the compiler.** The timings below are from a build with LLVM auto-vectorization turned off, so that every method runs one operation at a time. With the compiler default (auto-vectorization on), the methods gain very unequally: at N = 360 from nothing (`LEAR` with the gradient, `CLENSHAW` without) to 1.4–1.7 times (`LEGENDRE_a`, `LEGENDRE_b`, `METRIS`) and 2.3 times (`CLENSHAW` with the gradient). This reorders the methods: `LEAR` with the gradient goes from 12% faster than `LEGENDRE_a` to 39% slower, and `HOLMES` for the acceleration from 14% faster to 6% slower. Both sets of results are given; see [Auto-vectorization](#auto-vectorization). The CPU-time ratios were confirmed with the processor's cycle counter, to within 3%.
- **Acceleration only: `LEAR` and `HOLMES` are the fastest.** From N = 100 on, both need 13–18% less CPU time than `LEGENDRE_a` (104 vs 122 µs per point at N = 360), and they are tied with each other from N = 50 on. At N = 10, and on two of the three point sets at N = 20, `LEAR` is fastest alone (15% less than `LEGENDRE_a`).
- **Acceleration + gravity gradient: `HOLMES` is the fastest, by a wide margin.** From N = 100 on it needs about a third less CPU time than `LEGENDRE_a` (0.64–0.67 times; 170 vs 257 µs per point at N = 360). It wins outright at every N from 20 on, on all three point sets; at N = 10 it is tied with `LEAR`. It recurs on the same functions as `PINES` and `GOTTLIEB` but does not overflow up to N = 2000, because it scales them by a fixed factor.
- **`METRIS` and `LEAR` come next with the gradient.** From N = 100 to 1500 both take 0.86–0.90 times as long as `LEGENDRE_a` (222 and 227 µs per point at N = 360); at N = 2000, 0.93 and 0.98 times. Below N = 100 `METRIS` loses its lead (1.00–1.16 times `LEGENDRE_a`), while `LEAR` stays at 0.90–0.94.
- **`METRIS` for the acceleration alone** is about as fast as `LEGENDRE_a` from N = 100 to 1500 (2–5% slower), and 8–15% slower at N ≤ 50 and at N = 2000.
- **`CUNNINGHAM` and `CLENSHAW` are the slowest.** `CUNNINGHAM` takes 1.4–1.6 times as long as `LEGENDRE_a` for the acceleration and 1.6–1.7 times with the gradient, at every N. It uses the same solid harmonics as `METRIS` but applies the derivative operators at every point; `METRIS` does this once, at setup. At N = 360 `CUNNINGHAM` takes 1.4 times as long as `METRIS` for the acceleration and 1.9 times as long with the gradient. `CLENSHAW` takes 1.5–1.7 times as long as `LEGENDRE_a` for both outputs. Its speed advantage needs a grid, which is not tested here.
- **`PINES` and `GOTTLIEB`** are close to `LEGENDRE_a` from N = 100 on: 5–10% slower for the acceleration, and with the gradient 2–5% (`PINES`) and 1% (`GOTTLIEB`) slower. At N = 10 they take 10–24% longer.
- **`LEGENDRE_b`** (the whole Legendre triangle is stored first, then summed) takes 2–16% longer than `LEGENDRE_a`, most at N = 10.
- The ranking is the same on all three point sets.
- **The ratio curves are not monotone in N, but the times are.** Every method's CPU time per point grows steadily with N. The ratios to `LEGENDRE_a` move up and down by a few per cent from one N to the next. Part of this is real (fixed costs at small N, the cache at large N), part is timing noise of about ±5% per cell. See [Reading the CPU-time curves](#reading-the-cpu-time-curves).

## Results

The results below come from one run on 12 November 2026, on one performance core of an Apple M3 Pro (macOS 26.7, rustc 1.87.0, opt-level 3, LTO, auto-vectorization off). The timing phase took 21.7 minutes. Other jobs were running on the machine (load average 4–5). The A/A probe differed from 1 by 1.0% (median) and 4.5% (maximum). Differences below about 5% are therefore not resolved; the winner rule takes this into account.

The same run with the auto-vectorized build (same day, 20.0 minutes, load average 4.5–7.5, A/A probe 0.3% and 4.4%) is in [results_simd/](results_simd/), with the same tables and figures. The two builds are compared directly in [Auto-vectorization](#auto-vectorization).

### Figures

**[results/fig_accuracy.pdf](results/fig_accuracy.pdf)** has one panel per method (nine panels). Each panel shows the max. relative error over all test points against the maximum degree N:

- blue solid line: the acceleration;
- orange dashed line: the gravity gradient;
- dotted line: the pass threshold 10⁻¹¹;
- grey ×: a degree at which the method overflowed.

All panels share the same axes, so the heights of the curves can be compared from one panel to the next. The legend is below the panels.

**[results/fig_speed.pdf](results/fig_speed.pdf)** shows the CPU time. Each curve is the median over the three point sets, and the shaded band is their range (lowest to highest).

- **Panel (a)** shows the CPU time per point of `LEGENDRE_a`, in µs: blue for the acceleration, orange for the acceleration + gravity gradient. The time grows roughly as N².
- **Panels (b)**, one per other method (eight panels), show its CPU time divided by that of `LEGENDRE_a`, for the same output and the same N. Below 1 means faster than `LEGENDRE_a`. The light grey band is ±5%, the size of difference the timing cannot resolve. A grey × marks an overflow; those cases are not timed.

To get a method's time in µs, multiply its ratio in (b) by the value in (a), or read the tables below.

Two more figures put all the methods on one axis, so they can be compared directly. Each has one panel per output, and each method has its own colour and marker.

- **[results/fig_speed_all.pdf](results/fig_speed_all.pdf)** shows every method's CPU time divided by that of `LEGENDRE_a` (the black line at 1): the median over the three point sets, with a vertical bar for their range. Each curve is labelled at its right end, and the grey band is again ±5%. Panel (a) is the acceleration; panel (b) is the acceleration + gravity gradient.
- **[results/fig_accuracy_all.pdf](results/fig_accuracy_all.pdf)** shows every method's max. relative error. Panel (a) is the acceleration; panel (b) is the gravity gradient. `LEGENDRE_a` and `LEGENDRE_b` give identical results, so they share one curve. The curves cross each other often, because no method is consistently more accurate than the others.

In both figures, a curve stops where its method overflows.

**[results/fig_simd.pdf](results/fig_simd.pdf)** has one panel per method. Each panel shows the speed-up from auto-vectorization (the CPU time of the scalar build divided by that of the auto-vectorized build) against N: blue for the acceleration, orange dashed for the acceleration + gravity gradient, with the 95% interval as a band. Above 1 means the auto-vectorized build is faster. A grey × marks an overflow.

[results_simd/](results_simd/) has the first four figures again, for the auto-vectorized build.

Run `python plots.py` to redraw all the figures from the saved tables, without rerunning anything.

### Accuracy

Max. relative error = largest value, over all test points, of |computed - reference| / |reference - point-mass term| (vector norm for the acceleration, Frobenius norm for the gravity gradient). Dividing by the field minus its point-mass part GM/r^2 measures the error in the part the harmonics add. Test points: 3 sets x 48 points (ORBIT1, UNIFORM, POLAR). Reference = double-double synthesis (sph_rs ref). Pass: error <= 1e-11 and <= 10 x the best method; overflow (Inf/NaN) fails the method at that N.

**Acceleration**

| Method | N=10 | N=20 | N=50 | N=100 | N=200 | N=360 | N=720 | N=1000 | N=1500 | N=2000 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| LEGENDRE_a | 4.6e-13 | 6.7e-13 | 9.2e-13 | 1.1e-12 | 1.7e-12 | 2.0e-12 | 2.4e-12 | 2.3e-12 | 2.3e-12 | 2.3e-12 |
| LEGENDRE_b | 4.6e-13 | 6.7e-13 | 9.2e-13 | 1.1e-12 | 1.7e-12 | 2.0e-12 | 2.4e-12 | 2.3e-12 | 2.3e-12 | 2.3e-12 |
| PINES | 4.1e-13 | 5.5e-13 | 7.3e-13 | 1.5e-12 | 1.5e-12 | 2.0e-12 | 2.0e-12 | 2.1e-12 | overflow | overflow |
| GOTTLIEB | 6.6e-13 | 5.6e-13 | 7.3e-13 | 9.5e-13 | 1.5e-12 | 1.8e-12 | 2.6e-12 | 2.6e-12 | overflow | overflow |
| METRIS | 3.6e-13 | 4.7e-13 | 6.7e-13 | 9.2e-13 | 1.6e-12 | 1.7e-12 | 3.1e-12 | 3.1e-12 | 3.1e-12 | 3.1e-12 |
| CLENSHAW | 5.0e-13 | 6.8e-13 | 1.1e-12 | 1.4e-12 | 1.7e-12 | 2.2e-12 | 2.7e-12 | 2.7e-12 | 2.7e-12 | overflow |
| LEAR | 5.0e-13 | 7.4e-13 | 7.5e-13 | 9.9e-13 | 1.3e-12 | 1.7e-12 | 2.7e-12 | 2.7e-12 | 2.7e-12 | 2.7e-12 |
| HOLMES | 4.6e-13 | 6.7e-13 | 8.1e-13 | 1.0e-12 | 1.2e-12 | 1.7e-12 | 2.6e-12 | 2.6e-12 | 2.6e-12 | 2.6e-12 |
| CUNNINGHAM | 4.1e-13 | 5.2e-13 | 6.8e-13 | 5.8e-13 | 8.0e-13 | 1.1e-12 | 1.5e-12 | 1.5e-12 | 1.5e-12 | 1.5e-12 |

**Gravity gradient**

| Method | N=10 | N=20 | N=50 | N=100 | N=200 | N=360 | N=720 | N=1000 | N=1500 | N=2000 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| LEGENDRE_a | 2.3e-13 | 2.9e-13 | 3.2e-13 | 4.2e-13 | 5.5e-13 | 6.2e-13 | 1.1e-12 | 1.0e-12 | 1.0e-12 | 1.0e-12 |
| LEGENDRE_b | 2.3e-13 | 2.9e-13 | 3.2e-13 | 4.2e-13 | 5.5e-13 | 6.2e-13 | 1.1e-12 | 1.0e-12 | 1.0e-12 | 1.0e-12 |
| PINES | 2.5e-13 | 2.2e-13 | 4.4e-13 | 5.1e-13 | 6.1e-13 | 9.4e-13 | 1.2e-12 | 1.1e-12 | overflow | overflow |
| GOTTLIEB | 2.5e-13 | 3.1e-13 | 2.9e-13 | 5.1e-13 | 6.2e-13 | 9.1e-13 | 1.1e-12 | 1.2e-12 | overflow | overflow |
| METRIS | 3.0e-13 | 3.4e-13 | 3.4e-13 | 3.1e-13 | 5.6e-13 | 6.0e-13 | 8.2e-13 | 9.6e-13 | 9.6e-13 | 9.6e-13 |
| CLENSHAW | 2.5e-13 | 3.1e-13 | 4.7e-13 | 4.5e-13 | 5.7e-13 | 9.1e-13 | 9.1e-13 | 1.1e-12 | 1.1e-12 | overflow |
| LEAR | 2.3e-13 | 3.2e-13 | 4.0e-13 | 4.1e-13 | 4.3e-13 | 6.5e-13 | 1.2e-12 | 1.2e-12 | 1.2e-12 | 1.2e-12 |
| HOLMES | 2.3e-13 | 2.5e-13 | 3.6e-13 | 4.7e-13 | 6.1e-13 | 7.0e-13 | 8.0e-13 | 8.3e-13 | 8.3e-13 | 8.3e-13 |
| CUNNINGHAM | 3.6e-13 | 3.9e-13 | 4.9e-13 | 4.7e-13 | 5.4e-13 | 5.9e-13 | 1.3e-12 | 1.2e-12 | 1.2e-12 | 1.2e-12 |

Reference check (two independent double-double formulations): max difference 4e-28.

### CPU time

CPU time per point in microseconds (ORBIT1 points; median of the blocks, normalised by the calibration loop). **bold** = fastest; tied methods (difference not resolvable, see winners.csv) marked with *. x = failed the accuracy check, - = overflow (not timed).

**Acceleration**

| Method | N=10 | N=20 | N=50 | N=100 | N=200 | N=360 | N=720 | N=1000 | N=1500 | N=2000 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| LEGENDRE_a | 0.160 | 0.482 | 2.48* | 10.2 | 40.3 | 122 | 485 | 937 | 2267 | 3993 |
| LEGENDRE_b | 0.187 | 0.546 | 2.82 | 11.2 | 41.8 | 130 | 515 | 1001 | 2367 | 4226 |
| PINES | 0.197 | 0.555 | 2.77 | 11.1 | 42.1 | 131 | 515 | 1022 | - | - |
| GOTTLIEB | 0.179 | 0.544 | 2.82 | 10.9 | 42.2 | 132 | 520 | 1009 | - | - |
| METRIS | 0.187 | 0.555 | 2.69 | 10.6 | 40.2 | 127 | 499 | 982 | 2322 | 4366 |
| CLENSHAW | 0.240 | 0.771 | 4.35 | 16.0 | 64.1 | 209 | 800 | 1551 | 3530 | - |
| LEAR | **0.136** | **0.411** | **2.17** | **8.53** | **32.7** | 104* | 419* | **806** | **1864** | 3397* |
| HOLMES | 0.156 | 0.453 | 2.25* | 8.60* | 32.8* | **104** | **418** | 810* | 1892* | **3352** |
| CUNNINGHAM | 0.230 | 0.686 | 3.75 | 15.0 | 57.6 | 181 | 723 | 1444 | 3312 | 6272 |

**Acceleration + gravity gradient**

| Method | N=10 | N=20 | N=50 | N=100 | N=200 | N=360 | N=720 | N=1000 | N=1500 | N=2000 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| LEGENDRE_a | 0.295 | 0.964 | 5.30 | 21.2 | 81.9 | 257 | 1021 | 2054 | 4689 | 8214 |
| LEGENDRE_b | 0.324 | 1.03 | 5.66 | 22.8 | 86.4 | 266 | 1069 | 2090 | 4857 | 8394 |
| PINES | 0.366 | 1.07 | 5.69 | 21.6 | 84.8 | 263 | 1054 | 2233 | - | - |
| GOTTLIEB | 0.334 | 1.06 | 5.42 | 21.1 | 83.8 | 259 | 1031 | 2072 | - | - |
| METRIS | 0.341 | 0.967 | 5.51 | 18.4 | 72.7 | 222 | 904 | 1787 | 4113 | 7618 |
| CLENSHAW | 0.501 | 1.57 | 8.43 | 32.8 | 136 | 424 | 1691 | 3222 | 7088 | - |
| LEAR | 0.277* | 0.876 | 4.77 | 18.2 | 70.4 | 227 | 920 | 1852 | 4089 | 8010 |
| HOLMES | **0.269** | **0.765** | **3.74** | **14.0** | **53.7** | **170** | **679** | **1335** | **2981** | **5395** |
| CUNNINGHAM | 0.472 | 1.54 | 8.56 | 34.4 | 133 | 424 | 1701 | 3299 | 7577 | 13858 |

### Fastest method per N (all point sets)

| Output | Points | N=10 | N=20 | N=50 | N=100 | N=200 | N=360 | N=720 | N=1000 | N=1500 | N=2000 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Acc. | ORBIT1 | LEAR | LEAR | tie: LEAR+HOLMES+LEGENDRE_a | tie: LEAR+HOLMES | tie: LEAR+HOLMES | tie: HOLMES+LEAR | tie: HOLMES+LEAR | tie: LEAR+HOLMES | tie: LEAR+HOLMES | tie: HOLMES+LEAR |
| Acc. | UNIFORM | LEAR | tie: LEAR+HOLMES | tie: LEAR+HOLMES | tie: LEAR+HOLMES | tie: LEAR+HOLMES | tie: LEAR+HOLMES | tie: HOLMES+LEAR | tie: LEAR+HOLMES | tie: HOLMES+LEAR | tie: HOLMES+LEAR |
| Acc. | POLAR | LEAR | LEAR | tie: LEAR+HOLMES | tie: LEAR+HOLMES | tie: LEAR+HOLMES | tie: HOLMES+LEAR | tie: LEAR+HOLMES | tie: LEAR+HOLMES | tie: HOLMES+LEAR | tie: LEAR+HOLMES |
| Acc. + grad. | ORBIT1 | tie: HOLMES+LEAR | HOLMES | HOLMES | HOLMES | HOLMES | HOLMES | HOLMES | HOLMES | HOLMES | HOLMES |
| Acc. + grad. | UNIFORM | tie: HOLMES+LEAR | HOLMES | HOLMES | HOLMES | HOLMES | HOLMES | HOLMES | HOLMES | HOLMES | HOLMES |
| Acc. + grad. | POLAR | tie: HOLMES+LEAR+LEGENDRE_a | HOLMES | HOLMES | HOLMES | HOLMES | HOLMES | HOLMES | HOLMES | HOLMES | HOLMES |

A/A probe (LEGENDRE_a timed against an identical second copy of itself): median |ratio - 1| = 1.0%, max 4.5%; a win needs the 95% interval of the time ratio above 1 + max(5%, 2 |A/A - 1|).

Point-set check (the cost of a point does not depend on where it is, so the 3 point sets are repeat measurements): spread (max - min) / median of the time ratio to LEGENDRE_a over the sets = 2.1% (median), 14.0% (max, LEAR level 2 N=2000).

### Reading the CPU-time curves

The CPU time per point grows steadily with N for every method (panel (a) of `fig_speed`, and the tables above). The ratios to `LEGENDRE_a` do not: they move up and down from one N to the next. Four things cause this.

- **Fixed costs at small N.** At N ≤ 50 a good part of the time goes to work done once per point (trigonometric functions, starting values of the recursions), and this work differs between methods. As N grows, the work per coefficient takes over. This is why, for the acceleration, `PINES`, `LEGENDRE_b` and `METRIS` start 15–23% above `LEGENDRE_a` at N = 10 and fall to 2–6% above it by N = 200, and why `HOLMES` with the gradient falls from 0.91 to 0.66 times `LEGENDRE_a`.
- **The cache at large N.** The coefficients C_nm and S_nm that a method reads at every point take about 8 N² bytes: 8 MB at N = 1000, 18 MB at N = 1500. The L2 cache of the performance cores holds 16 MB, so from N = 1500 on they no longer fit and are read from main memory at every point. `METRIS` reads more: its derivative operators are applied to the coefficients at setup, so it reads 3 coefficient sets for the acceleration and 9 with the gradient (about 24 N² and 72 N² bytes). These leave the cache from N ≈ 800 and N ≈ 470. Its ratio to `LEGENDRE_a` rises at large N (acceleration from 1.02–1.05 at N = 200–1500 to 1.11 at N = 2000; with the gradient from 0.86 at N = 360 to 0.93 at N = 2000), most likely because of this extra memory traffic.
- **`LEAR` with the gradient at N = 2000.** At the POLAR points it takes about 16% longer than at the other two point sets. This is the largest point-set spread in the run (14.0% in the ratio to `LEGENDRE_a`), and every earlier run, with either build, showed it too. The kernel has no data-dependent branches and the slowdown is not caused by subnormal numbers; we have not found the cause.
- **Timing noise.** Other jobs on the machine share the L2 cache and change the clock speed. The calibration loop removes most of the clock-speed changes: the three point sets, which do identical work, disagree by 2.1% (median) after calibration, against 12% without it. It cannot see cache sharing, so each cell keeps a noise of about ±5% (A/A probe up to 4.5%).

To read the curves:

- Treat a ratio within ±5% of 1 as a tie with `LEGENDRE_a`.
- Treat a step of less than about 5% between neighbouring N as noise, unless it is wider than the range band and continues over several N.
- Use the median over the three point sets (the curves in the figures) rather than one point set. The tables give ORBIT1 only, so they show more of the noise.

### Auto-vectorization

The M3 processor can add or multiply two double-precision numbers with one instruction (SIMD). The Rust compiler (LLVM) decides by itself where to use such instructions; this is called auto-vectorization. It applies them to the inner loops of some methods and not of others, so with the compiler default a part of the difference between two methods comes from the compiler and not from the algorithm.

For this reason the code is built twice:

- **Scalar build** (main results, [results/](results/)): auto-vectorization off (`RUSTFLAGS="-C no-vectorize-loops -C no-vectorize-slp"`). The binary contains no two-wide floating-point instruction.
- **Auto-vectorized build** ([results_simd/](results_simd/)): the compiler default.

The two builds return bit-identical outputs: the accuracy files of the two runs are byte-identical.

The tables below compare the two builds at N = 360. "Timed" is the CPU-time ratio to `LEGENDRE_a` from the full timing run of each build (median over the three point sets). "Cycles" is the same ratio measured independently with the processor's cycle counter (`/usr/bin/time -l`, ORBIT1 points, [hw_counters.py](hw_counters.py), [results/hw_counters.csv](results/hw_counters.csv)). The last column is the speed-up from auto-vectorization: the cycles per point of the scalar build divided by those of the auto-vectorized build.

**Acceleration, N = 360**

| Method | Scalar: timed | Scalar: cycles | Auto-vectorized: timed | Auto-vectorized: cycles | Speed-up |
|---|---:|---:|---:|---:|---:|
| `LEGENDRE_a` | 1.00 | 1.00 | 1.00 | 1.00 | 1.46 |
| `LEGENDRE_b` | 1.07 | 1.06 | 1.05 | 1.06 | 1.47 |
| `PINES` | 1.07 | 1.06 | 1.15 | 1.14 | 1.37 |
| `GOTTLIEB` | 1.08 | 1.07 | 1.11 | 1.11 | 1.41 |
| `METRIS` | 1.04 | 1.03 | 1.05 | 1.05 | 1.44 |
| `CLENSHAW` | 1.68 | 1.64 | 2.34 | 2.34 | 1.02 |
| `LEAR` | 0.85 | 0.85 | 0.87 | 0.88 | 1.41 |
| `HOLMES` | 0.86 | 0.85 | 1.06 | 1.07 | 1.17 |
| `CUNNINGHAM` | 1.50 | 1.49 | 1.85 | 1.86 | 1.17 |

**Acceleration + gravity gradient, N = 360**

| Method | Scalar: timed | Scalar: cycles | Auto-vectorized: timed | Auto-vectorized: cycles | Speed-up |
|---|---:|---:|---:|---:|---:|
| `LEGENDRE_a` | 1.00 | 1.00 | 1.00 | 1.00 | 1.53 |
| `LEGENDRE_b` | 1.04 | 1.04 | 1.03 | 1.03 | 1.54 |
| `PINES` | 1.02 | 1.02 | 1.24 | 1.24 | 1.25 |
| `GOTTLIEB` | 1.01 | 1.00 | 1.19 | 1.19 | 1.29 |
| `METRIS` | 0.86 | 0.87 | 0.78 | 0.78 | 1.69 |
| `CLENSHAW` | 1.65 | 1.64 | 1.10 | 1.10 | 2.29 |
| `LEAR` | 0.88 | 0.88 | 1.39 | 1.39 | 0.97 |
| `HOLMES` | 0.66 | 0.65 | 0.92 | 0.92 | 1.09 |
| `CUNNINGHAM` | 1.65 | 1.65 | 2.35 | 2.35 | 1.07 |

What the tables show:

- **The timings are confirmed.** The timed ratios and the cycle ratios agree to within 3% (scalar build) and 1% (auto-vectorized build) for every method and both outputs.
- **The gain from auto-vectorization is very unequal.** `LEGENDRE_a`, `LEGENDRE_b` and `METRIS` gain 1.4–1.7 times; `CUNNINGHAM` and `HOLMES` gain 1.1–1.2 times; `LEAR` gains 1.4 times for the acceleration and nothing with the gradient; `CLENSHAW` gains nothing for the acceleration and 2.3 times with the gradient.
- **The ranking changes with the build.** With the gradient, `HOLMES` is fastest in the scalar build (0.66) and `METRIS` in the auto-vectorized build (0.78). `LEAR` with the gradient moves from 0.88 to 1.39, and `CUNNINGHAM` from 1.65 to 2.35. For the acceleration, `LEAR` is fastest in both builds, but `HOLMES` moves from 0.86 to 1.06.
- **`CUNNINGHAM` is slow in both builds.** In the scalar build its cost matches its arithmetic: it carries the real and the imaginary part of the solid harmonics through the recursion and pairs each coefficient with three neighbouring harmonics (nine with the gradient), about 1.6 times (1.7 with the gradient) the floating-point operations of `LEGENDRE_a` per (n, m) term.

The scalar build is the main result because it compares the algorithms as written. The auto-vectorized build shows what this compiler makes of this code; another compiler, compiler version or processor can vectorize different loops and give a different ranking. A method that gains little here (`HOLMES`, `CUNNINGHAM`) could gain more if its loops were rewritten for SIMD by hand; that was not done.

[results/fig_simd.pdf](results/fig_simd.pdf) and [results/simd.csv](results/simd.csv) give the speed-up at every N, from a separate run in which the two builds were timed side by side (ORBIT1 points, 3 rounds, raw CPU time, 12.4 minutes at load average 5). Raw time is used there because the calibration loop is itself auto-vectorized (8.7 ns per iteration in the scalar build, 7.8 ns in the other), so it cannot serve as a common unit. These speed-ups are noisy, ±5–14% per point; their levels agree with the cycle counts above, and the steps between neighbouring N are noise. The speed-up is smaller at N = 10 and 20 (0.9–1.5 for all methods but one), where the loops are short.

### Why `PINES`, `GOTTLIEB` and `CLENSHAW` overflow

`PINES` and `GOTTLIEB` recur on the Helmholtz polynomials Ā_nm(u) = P̄_nm(sin φ) / cos^m φ, with u = sin φ. Near the poles these polynomials grow very fast with the degree. Their largest values, found numerically, are:

| N | Largest log10 \|Ā_nm\| | Where it passes the double-precision limit (about 10³⁰⁸) |
|---|---|---|
| 1000 | 209, at the pole | Nowhere |
| 1500 | 314 at the pole; 308 at u = 0.99 | Latitudes above about 82° |
| 2000 | 418 at the pole; 343 at u = 0.9; 274 at u = 0.8 | Latitudes above about 58° |

All three point sets have points above 82° latitude (ORBIT1 is a near-polar orbit), so both methods fail in every set at N = 1500 and 2000.

`CLENSHAW` forms, for each order m, the sum over n of C_nm (a/r)^(n−m) P̄_nm / P̄_mm, and multiplies it by the sectoral term (a/r)^m P̄_mm only at the end. The ratio P̄_nm / P̄_mm grows like Ā_nm. Inside the sum it is multiplied by the small coefficients and by (a/r)^(n−m), which delays the overflow: in these tests the sum stays finite up to N = 1500 and overflows at N = 2000.

The other six methods do not overflow up to N = 2000:

- `LEGENDRE_a`, `LEGENDRE_b` and `LEAR` keep the factor cos^m φ inside their Legendre values (they divide by at most cos² φ).
- `METRIS` and `CUNNINGHAM` work with solid harmonics, which include the factor (a/r)^(n+1).
- `HOLMES` recurs on the same function as `PINES`, P̄_nm / cos^m φ, but every value is multiplied by the fixed factor 2⁻⁹³⁰ ≈ 10⁻²⁸⁰. The largest value at N = 2000, 10⁴¹⁸ near the pole, is therefore stored as about 10¹³⁸. The factor is removed once, at the end.

A fixed scale factor, as in `HOLMES`, is the simplest cure. Scaled numbers with their own exponent, such as the X-numbers of Fukushima (2012), work at any degree but cost more time. No kernel here uses X-numbers.

### Differences from Fantino & Casotto (2009)

The paper timed complete programmes on 2009 hardware, including file input and output, and its `LEGENDRE_b` ran through Matlab. Here only the synthesis kernels are timed, and all of them follow the same contract. The main differences:

- **`CLENSHAW`:** the paper found it about 2 times faster than `LEGENDRE_a` at scattered points. Here it is slower: it takes 1.5–1.7 times as long, for both outputs.
- **The other methods:** in the paper's Table 7, `LEGENDRE_b`, `PINES` and `METRIS` take 4.8, 3.0 and 1.6 times as long as `LEGENDRE_a`. Here, for the acceleration + gravity gradient at N = 360, they take 1.04, 1.02 and 0.86 times as long.
- **New methods:** `GOTTLIEB`, `LEAR`, `HOLMES` and `CUNNINGHAM` are not in the paper. (The paper's `METRIS` uses the same solid harmonics as `CUNNINGHAM`.)

These differences most likely come from the 2009 implementations, compilers and hardware rather than from the algorithms. The compiler alone can move a ratio a long way: see [Auto-vectorization](#auto-vectorization).

## Methods

| Method | Functions | Coordinates | How derivatives are obtained | Source |
|---|---|---|---|---|
| `LEGENDRE_a` | Fully normalised Legendre functions, column recursion over the degree | Spherical (λ, φ, r) | Analytic derivatives of each term in the local east/north/radial frame, then rotated to Cartesian. Two rolling rows, terms summed as they are computed. | Fantino & Casotto (2009) |
| `LEGENDRE_b` | Same as `LEGENDRE_a` | Spherical | Same mathematics, but the whole triangle of Legendre values is stored first and summed afterwards (the paper's "Matlab + DLL" programme) | Fantino & Casotto (2009) |
| `PINES` | Normalised Helmholtz polynomials A_nm(u) | Cartesian direction cosines (s, t, u) | Derivatives from rows n+1 and n+2 of A_nm | Pines (1973), normalised as in Eckman et al. (2016) |
| `GOTTLIEB` | Same A_nm(u) as `PINES` | Cartesian (s, t, u) | ∂A_nm/∂u = k_nm A_n,m+1, so only row n is needed | Gottlieb (1993) |
| `METRIS` | Fully normalised solid harmonics (Cunningham) | Cartesian (x, y, z) | The derivative operators are applied to the coefficients once, at setup: 3 coefficient sets for the acceleration and 9 for acceleration + gradient. Each point then needs only plain solid-harmonic sums. | Metris et al. (1999) |
| `CLENSHAW` | Clenshaw summation over the degree, for each order m | Spherical | Clenshaw sums of the derivative series, then the same local-frame terms as `LEGENDRE_a`. This is the scattered-point version; the per-parallel reuse that makes Clenshaw fast on grids does not apply here. | Fantino & Casotto (2009) |
| `LEAR` | P̄_n0 and dP̄_n0/d(sin φ) for m = 0; P̄_nm / cos φ and dP̄_nm/dφ for m ≥ 1 | Spherical | Lear's own recursions for the acceleration. The gravity gradient is our extension (below). | Lear (1986), normalised as in Eckman et al. (2016, app. C.7) |
| `HOLMES` | P̄_nm / cos^m φ, scaled by 2⁻⁹³⁰, recurred across the orders m at fixed degree n (m = N down to 0) | Spherical | The same local-frame terms as `LEGENDRE_a`. The orders are summed by Horner's rule in cos φ, which puts the factor cos^m φ back (below). | Holmes & Featherstone (2002), organised as in Orekit |
| `CUNNINGHAM` | Same solid harmonics as `METRIS` | Cartesian (x, y, z) | The derivative operators are applied at every point, to the original coefficients C̄_nm, S̄_nm: each term is formed from the harmonics of degree n+1 (acceleration) and n+2 (gradient) next to it. No transformed coefficient sets. | Cunningham (1970), as in Montenbruck & Gill (2000, Sect. 3.2.4–3.2.5) |

**Pole safety.** None of the methods divides by cos φ:

- The spherical methods `LEGENDRE_a`, `LEGENDRE_b` and `CLENSHAW` recur on Z_nm = P̄_nm / cos(φ)^min(m,2). With this choice the 1/cos φ and 1/cos² φ factors of the local-frame formulas cancel analytically. The cost is about two extra multiplications per term.
- `LEAR` works with P̄_nm / cos φ for the acceleration and adds P̄_nm / cos² φ for the gravity gradient.
- `HOLMES` works with P̄_nm / cos^m φ and multiplies the factor cos^m φ back in by Horner's rule. The derivative terms that would need cos^(m−1) φ or cos^(m−2) φ have their own Horner chains.
- `PINES`, `GOTTLIEB`, `METRIS` and `CUNNINGHAM` are Cartesian and have no singularity.

### LEAR gravity gradient (our extension)

Lear (1986) and Eckman et al. (2016) give only the acceleration. To compare LEAR on both outputs, we derived the gravity gradient in Lear's style. The extension keeps three properties of Lear's method:

- It uses only column recursions, as Lear does.
- It contains no division by cos φ.
- It does not use the Legendre differential equation.

In the local frame, each (n, m) term contributes six tensor components:

| Component | Term (× μ (a/r)^n / r³, with P = P̄_nm, ′ = d/dφ, c = cos φ, s = sin φ) |
|---|---|
| ee | −m² P/c² − (s/c) P′ − (n+1) P |
| en | m e₅ (P′/c + s P/c²) |
| er | −(n+2) m e₅ P/c |
| nn | P″ − (n+1) P |
| nr | −(n+2) P′ |
| rr | (n+1)(n+2) P |

The components er, nr and rr, and the (n+1)P parts, reuse the sums already formed for the acceleration. The 1/c and 1/c² factors of ee, en and nn are removed separately for each kind of column:

- **m = 0.** Take one more derivative with respect to s, by the same recursion Lear uses for dP̄_n0/ds.
- **m = 1.** P̄_n1 / c is a polynomial Q(s). Differentiating its column recursion gives Q′ and Q″.
- **m ≥ 2.** Carry q = P̄_nm / c² with the same column recursion as Lear's P̄_nm / c. Dividing Lear's derivative identity by c² gives G = P′/c. Differentiating that identity once more gives P″. Each needs only the neighbouring degree n−1.

The kernel is `Lear` in [sph_rs/src/bench.rs](sph_rs/src/bench.rs); the comment block above it has the full formulas. It passes the accuracy check at every N and on every point set, including points 10⁻⁸ degrees from the poles. Its errors are of the same size as those of `LEGENDRE_a`.

### HOLMES (Holmes–Featherstone)

Orekit evaluates the geopotential with the method of Holmes & Featherstone (2002). `HOLMES` follows Orekit's organisation of the recursion, with s = sin φ and c = cos φ:

- **Functions.** It recurs on P̃_nm = P̄_nm / c^m, the function `PINES` calls Ā_nm. Without the factor c^m these values grow very fast near the poles, so every value is multiplied by 2⁻⁹³⁰ ≈ 10⁻²⁸⁰ (Orekit's factor; the paper uses 10⁻²⁸⁰). With this factor nothing overflows up to N = 2000.
- **Recursion.** For each order m, from N down to 0, the values for n = m+1 … N follow from the two previous orders: P̃_nm = g_nm s P̃_n,m+1 − h_nm c² P̃_n,m+2, starting from the sectoral value P̃_mm. Holmes & Featherstone call this the modified forward row method. The constants g_nm and h_nm are stored in the order the loop reads them, as in Orekit.
- **Sums.** For each order, the sums over n are formed first. The orders are then combined by Horner's rule in c (Σ_m c^m T_m), which puts the factor c^m back. The scale factor is removed once, at the end.
- **Derivatives.** They use the same local-frame terms as `LEGENDRE_a`. Some terms carry c^(m−1) or c^(m−2) instead of c^m; they always come with a factor m or m(m−1). Those terms go into their own Horner chains in c, which start at m = 1 and m = 2. Nothing is divided by c, so the kernel works at the poles.

**Differences from Orekit.**

- Orekit divides by c: tan φ appears in its derivative of P̄_nm, and 1/c and 1/c² in its conversion from spherical to Cartesian derivatives. Its gradient is therefore singular at the poles. Our Horner chains avoid this.
- Orekit skips the degrees n = 0 and 1. `HOLMES` includes them, as every other method here does.

The kernel is `Hf` in [sph_rs/src/bench.rs](sph_rs/src/bench.rs). It passes the accuracy check at every N and on every point set, including points 10⁻⁸ degrees from the poles. Its errors are of the same size as those of `LEGENDRE_a`.

### CUNNINGHAM (Montenbruck & Gill)

Montenbruck & Gill (2000, Sect. 3.2.4–3.2.5) give Cunningham's (1970) method in the form most orbit programmes use; Scarabaeus is one of them. `CUNNINGHAM` follows their formulas, fully normalised:

- **Functions.** The solid harmonics V_nm + i W_nm = (a/r)^(n+1) P̄_nm(sin φ) e^(imλ), from the same Cartesian recursion as `METRIS`. No trigonometric functions are needed, and nothing is divided by cos φ.
- **Acceleration.** Each term of degree n and order m is formed directly from the original coefficients C̄_nm, S̄_nm and the harmonics V, W of degree n+1 and orders m−1, m and m+1 (M&G Eq. 3.33, normalised). The normalisation factors are square roots of small integers; they come from two short tables, √k and √(k(k+1)).
- **Gravity gradient.** Montenbruck & Gill write out the acceleration. Every derivative of a solid harmonic is again a combination of solid harmonics of the next degree, so we apply the same operators a second time. Each term then also needs the harmonics of degree n+2 and orders m−2 … m+2.
- **One sum for xx, yy and zz.** In this formulation the factor of ∂²/∂x² + ∂²/∂y² of each term is exactly that of −∂²/∂z², so the six tensor components come from five sums. This is not the Laplace shortcut forbidden by the contract (yy is not computed as −xx − zz): a separate sum for ∂²/∂x² + ∂²/∂y² would repeat the same products and give the same result.
- **Difference from `METRIS`.** `METRIS` applies the derivative operators to the coefficients once, at setup, so at each point a term is one stored coefficient times one harmonic per output component. `CUNNINGHAM` does this work at every point, which costs more multiplications per term, but it stores nothing besides C̄ and S̄: it reads 8 N² bytes of coefficients per point, like `LEGENDRE_a`, against 24 N² (acceleration) and 72 N² (with the gradient) for `METRIS`.

The kernel is `Cunn` in [sph_rs/src/bench.rs](sph_rs/src/bench.rs). It passes the accuracy check at every N and on every point set, including points 10⁻⁸ degrees from the poles. For the acceleration it has the smallest error of all methods from N = 100 on (1.1 × 10⁻¹² at N = 360, against 1.7–2.2 × 10⁻¹² for the others); for the gravity gradient its errors are of the same size as those of `LEGENDRE_a`. It is the slowest method that does not overflow: at N = 360 it takes 181 µs per point for the acceleration and 424 µs with the gradient, against 127 and 222 µs for `METRIS`.

## The contract

Every method obeys the same rules, so that the timings compare algorithms and not implementations:

- **Input:** body-fixed x, y, z.
- **Output:** the Cartesian body-fixed acceleration (ax, ay, az), or the acceleration plus the gravity-gradient tensor (xx, xy, xz, yy, yz, zz). The potential is not output.
- **Output choice:** made at compile time (a const generic), so the acceleration-only code does none of the gradient work. The acceleration-only output must be bit-identical to the first three components of the acceleration + gradient output.
- **No shortcuts:** all six tensor components are synthesised, with no Laplace shortcut and no use of the Legendre differential equation.
- **Shared building blocks:**
  - cos mλ and sin mλ come from the same complex recursion in every spherical method.
  - (a/r)^n comes from a running product.
  - Terms are first summed per degree, then the per-degree sums are added up. A single running sum over all terms loses about one digit at N = 360.
- **Plain code:** no allocation in the timed loop, no hand-written fused multiply-add, no fast-math and no flush-to-zero. Each kernel runs on one thread.
- **Setup excluded:** recursion tables and the `METRIS` coefficient operators are built once per model. This setup is timed separately (`setup_s` in `timing.csv`) and is not part of the time per point.

## Accuracy check

**Reference.** `sph_rs ref` computes the acceleration and the gravity gradient in double-double arithmetic, which carries about 32 significant digits. The Legendre values are stored as X-numbers (Fukushima 2012), so they cannot underflow or overflow at any degree or latitude. The reference uses two independent formulations:

- **A:** a spherical formulation.
- **B:** a Cartesian formulation, after Pines.

Their relative difference is below 4 × 10⁻²⁸ at every test point.

**Test points.** There are three sets of 48 points each:

- **ORBIT1:** samples of a 23-day, 250 km near-polar orbit.
- **UNIFORM:** uniform on the sphere, at altitudes of 250–1000 km.
- **POLAR:** the same altitudes, at colatitudes from 10⁻⁸ to 1 degree.

**Error measure.** For each method and test point:

relative error = |computed − reference| / |reference − point-mass term|

- |·| is the vector norm for the acceleration and the Frobenius norm for the gravity gradient.
- The point-mass term is −GM x / r³ (and its gradient). It is about 1000 times larger than the rest of the field. Dividing by the full field would therefore make every method look about 1000 times more accurate than it is in the part the harmonics add.
- The tables and figures report the **maximum relative error over all test points**.

**Pass rule.** A method passes at a given N if all four conditions hold:

- its maximum relative error is ≤ 10⁻¹¹;
- that error is ≤ 10 times the error of the best method at that N;
- every output is finite (an overflow to Inf or NaN fails the method at that N);
- its acceleration-only output equals the first three components of its acceleration + gradient output, bit for bit.

Only methods that pass can win the CPU-time comparison.

## Timing protocol

- **Clock:** thread CPU time (`clock_gettime(CLOCK_THREAD_CPUTIME_ID)`). On macOS the thread is given QoS class USER_INTERACTIVE, which keeps it on a performance core.
- **Build:** the timed kernels are compiled with opt-level 3 and LTO, and with LLVM's auto-vectorization turned off (`-C no-vectorize-loops -C no-vectorize-slp`). The compiler vectorizes some kernels and not others, and the gain differs a lot between methods, so with it turned on the ranking partly measures the compiler rather than the algorithm (see [Auto-vectorization](#auto-vectorization)). The two builds give bit-identical results.
- **Cells:** a cell is one combination of degree N, output and point set. Each method is timed in its own fresh process. Cells run in random order, and so do the methods within a cell, so that slow drift in the machine cannot line up with N or with a method.
- **Rounds:** the whole list of cells is run 3 times, each time in a new random order. Each (cell, method) is therefore measured in 3 separate processes spread over the whole run, and its blocks are pooled. This spreads slow changes in the machine's load over all cells and methods, so they cannot line up with one cell or one method. It does not remove them. Other jobs that share the 16 MB L2 cache of the performance cores slow the large-N kernels for seconds to minutes at a time, and the calibration loop (which runs in registers) cannot see this. What remains is about ±5% per cell; see [Reading the CPU-time curves](#reading-the-cpu-time-curves).
- **Inside a process:**
  - One warm-up block, which is discarded.
  - 4 timed blocks of at least 100 ms each (12 per cell and method over the 3 rounds). A block is a number of passes over K points, with K = 4·10⁶ / (number of coefficients), clamped to 32–4096 points.
  - A fixed calibration loop runs before and after every block: 8 independent multiply-add chains with no FMA.
  - Each block time is multiplied by c_ref / c_block. Here c_block is the mean of the two calibration runs around that block, and c_ref is the 5th percentile of all calibration runs in the whole timing phase (the machine at full speed). This removes changes in clock speed between blocks and between processes. It matters: without it, the three point sets disagree about 5 times more, 12% instead of 2.1% (median; see [Reading the CPU-time curves](#reading-the-cpu-time-curves)).
  - All blocks are kept. An earlier version dropped blocks whose two calibration runs differed by more than 3%; this made the worst cases slightly worse, not better.
- **Reported time:** the median over blocks, in µs per point, with a 95% bootstrap interval.
- **Point-set check:** the kernels have no data-dependent branches, so the cost of a point does not depend on where it is (one exception, below). The three point sets are therefore three repeat measurements of each (N, output, method). The figures show their median, with their range as a band; `summary.md` reports how far they disagree.
- **A/A probe:** a second compiled copy of `LEGENDRE_a` is timed as if it were a tenth method. Any difference it shows from `LEGENDRE_a` is measurement noise.
- **Winner rule:** the fastest passing method wins outright only if every other method is slower by more than δ, with δ = max(5%, 2 × |A/A ratio − 1|). "Slower by more than δ" means that the lower end of the 95% bootstrap interval of the time ratio is above 1 + δ. Otherwise the methods within δ are reported as tied.

## What operational software uses

This is from public documentation and source code, as of November 2026.

| Tool | Formulation | Closest method here |
|---|---|---|
| Scarabaeus (CCAR-ORCCA) | Normalised Cunningham solid-harmonic recursion, B_nm = V_nm + iW_nm (Montenbruck & Gill 2000, Sect. 3.2), with an analytic Jacobian. Its docstrings say "Pines/Eckman", but the code is the Cunningham recursion. | `CUNNINGHAM` |
| Basilisk | Normalised Pines | `PINES` |
| GMAT | Normalised Pines, with an optional gradient for the state transition matrix | `PINES` |
| Nyx | Pines, ported from GMAT | `PINES` |
| Tudat | Fully normalised Legendre functions in spherical coordinates, then a rotation | `LEGENDRE_a` |
| Orekit | Holmes–Featherstone (modified forward row method) | `HOLMES` (with pole-safe derivatives) |

None of these tools uses Gottlieb's or Lear's formulation. Both come from NASA Johnson Space Center; Eckman et al. (2016) normalise and compare them with Pines.

## Running it

**Requirements**

- Python 3 with NumPy and Matplotlib.
- A Rust toolchain (`cargo`). The crate has no dependencies.

`benchmark.py` builds the Rust binary on first use, and again whenever a source file is newer than the binary. There are two builds of the same source:

- **`scalar`** (the default): LLVM's auto-vectorization is turned off with `RUSTFLAGS="-C no-vectorize-loops -C no-vectorize-slp"`. It goes to `~/.cache/sph_rs/target`.
- **`simd`** (`--simd`): the compiler defaults, so LLVM may use 2-wide NEON instructions where it finds them. It goes to `~/.cache/sph_rs/target_simd`.

Set `SPH_RS_CACHE` to use another folder than `~/.cache/sph_rs`. `meta.json` records which build a run used.

**Commands**

```sh
python benchmark.py                   # accuracy check + timing + analysis -> results/
python benchmark.py --quick           # smoke test: N = 10, 50, 360, short blocks -> results_quick/
python benchmark.py --phase gate      # accuracy check only
python benchmark.py --phase time      # timing only (reuses results/accuracy.csv)
python benchmark.py --phase analyze   # tables and figures from the saved raw data
python benchmark.py --simd            # the same run with the auto-vectorized build -> results_simd/
python benchmark.py --phase simd      # both builds timed side by side (about 12 minutes) -> results/simd.csv, fig_simd.pdf
python hw_counters.py /tmp/hw 360 500 4000   # cycle and instruction counts per point (macOS), about 3 minutes
python plots.py                       # figures only, from the saved tables in results/ (no Rust needed)
```

| Option | Default | Meaning |
|---|---|---|
| `--N` | 10 20 50 100 200 360 720 1000 1500 2000 | Maximum degrees to test |
| `--blocks` | 4 | Timed blocks per process |
| `--rounds` | 3 | Passes over all cells; each (cell, method) gets one process per round |
| `--target-ms` | 100 | Minimum length of one block, in ms |
| `--k-gate` | 48 | Test points per point set in the accuracy check |
| `--threads` | 4 | Threads used by the reference (the timed kernels always use one) |
| `--seed` | 2009 | Seed for the run order and the timing points |
| `--simd` | off | Time the auto-vectorized build instead of the scalar one (default output folder `results_simd/`) |
| `--out` | `results/` | Output folder |

Every run saves its raw data in the output folder. `--phase analyze` rebuilds the tables and figures from `raw_timing.jsonl.gz` and `accuracy.csv`, and `plots.py` redraws the figures from the tables, so neither needs a new run or a Rust build. Use `--out` to keep several runs side by side.

The accuracy check takes about 1 minute. The timing phase starts about 1710 processes (60 cells × 3 rounds × up to 10 methods, including the A/A copy) and takes about 22 minutes. For clean timings, run on mains power, keep the machine awake (on macOS: `caffeinate -i -s python benchmark.py`) and close other heavy jobs. The A/A probe line at the end of `summary.md` shows how noisy the run was.

The binary (`~/.cache/sph_rs/target/release/sph_rs`, or `target_simd/release/sph_rs` for the other build) can also be called directly. To build the scalar one by hand, run `RUSTFLAGS="-C no-vectorize-loops -C no-vectorize-slp" cargo build --release` in `sph_rs/`.

```sh
sph_rs eval  coef.bin xyz.bin out.bin --method LEAR --level 2 --nmax 360   # one method
sph_rs ref   coef.bin xyz.bin out.bin --nmax 360                           # double-double reference
sph_rs bench coef.bin xyz.bin --level 1 --nmax 360 --methods LEAR,METRIS   # timing, JSON lines
```

The binary file formats are given at the top of [sph_rs/src/main.rs](sph_rs/src/main.rs) and [sph_rs/src/bench.rs](sph_rs/src/bench.rs).

## Files

| Path | Contents |
|---|---|
| [benchmark.py](benchmark.py) | Driver: model and test points, accuracy check, timing runs, analysis, tables |
| [plots.py](plots.py) | The figures, drawn from `accuracy.csv`, `timing.csv` and `simd.csv` |
| [sph_rs/src/main.rs](sph_rs/src/main.rs) | Command line, model loading, recursion tables, `METRIS` coefficient operators |
| [sph_rs/src/bench.rs](sph_rs/src/bench.rs) | The nine kernels, the `eval` and `bench` commands, the timing loop |
| [sph_rs/src/reference.rs](sph_rs/src/reference.rs) | Double-double reference with X-numbers (two formulations) |
| [results/accuracy.csv](results/accuracy.csv) | Max. relative error per N, point set and method; pass or fail |
| [results/timing.csv](results/timing.csv) | CPU time per point (µs) with 95% interval, per N, output, point set and method |
| [results/winners.csv](results/winners.csv) | Fastest method per cell, ties, A/A ratio |
| [results/summary.md](results/summary.md) | The tables shown above |
| [results/fig_accuracy.pdf](results/fig_accuracy.pdf), [results/fig_speed.pdf](results/fig_speed.pdf) | Figures, one panel per method |
| [results/fig_accuracy_all.pdf](results/fig_accuracy_all.pdf), [results/fig_speed_all.pdf](results/fig_speed_all.pdf) | Figures, all methods on one axis |
| [results/fig_simd.pdf](results/fig_simd.pdf), [results/simd.csv](results/simd.csv) | Speed-up from auto-vectorization, per method, N and output |
| [results/raw_timing.jsonl.gz](results/raw_timing.jsonl.gz), [results/raw_simd.jsonl.gz](results/raw_simd.jsonl.gz) | Every timed block, for re-analysis |
| [results/meta.json](results/meta.json) | Machine, compiler, build flags, run settings |
| [results_simd/](results_simd/) | The same files for the full run with the auto-vectorized build |
| [hw_counters.py](hw_counters.py), [results/hw_counters.csv](results/hw_counters.csv) | Cycles and instructions per point at N = 360, both builds (independent check of the timings) |
| [Fantino_Casotto_2009_JGeod_Methods_of_harmonic_synthesis.pdf](Fantino_Casotto_2009_JGeod_Methods_of_harmonic_synthesis.pdf) | The paper the comparison started from |

## Limits

- **One compiler.** The ranking by CPU time changes between the scalar and the auto-vectorized build (see [Auto-vectorization](#auto-vectorization)). Another compiler or compiler version can change it again. The differences of 10–15% between methods should be read with this in mind; the larger ones (`CUNNINGHAM` and `CLENSHAW` slow, `HOLMES` fast with the gradient in the scalar build) are outside it.
- **Busy machine.** The runs were made with other jobs on the machine (load average 4–5). The protocol corrects for this (calibration loop, three rounds, A/A probe), and the cycle counts agree with the timings, but the runs were not repeated on an idle machine.
- **One machine.** All times are from one Apple M3 Pro performance core. Ratios between methods may differ on other processors, especially between the Cartesian and the spherical methods.
- **Scalar code at scattered points.** Each point is evaluated on its own, on one thread. The main timings use no SIMD instructions at all. The auto-vectorized build uses 2-wide NEON inside some kernels, but no kernel is vectorized across points, and no kernel has hand-written SIMD. A method that vectorizes well by hand, or across several points at once, may rank differently. Grids, where Clenshaw summation and per-parallel reuse pay off, are not tested.
- **Synthetic model.** The coefficients are random numbers whose size follows Kaula's rule (10⁻⁵/n²), with the low-degree terms of EGM96. The CPU times do not depend on the coefficient values. The errors of a real model (EGM2008, GRACE) should be similar, but this is not tested.
- **Plain double precision.** No kernel uses X-numbers, and only `HOLMES` uses a scale factor. That is why `PINES`, `GOTTLIEB` and `CLENSHAW` overflow at high degree. Degrees above 2000 are not tested.
- **The gradient for `LEAR` is our own.** Lear (1986) and Eckman et al. (2016) give only the acceleration.

## References

- Cunningham, L. E. (1970). On the computation of the spherical harmonic terms needed during the numerical integration of the orbital motion of an artificial satellite. *Celestial Mechanics* 2, 207–216.
- Eckman, R. A., Brown, A. J., Adamo, D. R. (2016). *Normalization and implementation of three gravitational acceleration models.* NASA/TP-2016-218604.
- Fantino, E., Casotto, S. (2009). Methods of harmonic synthesis for global geopotential models and their first-, second- and third-order gradients. *Journal of Geodesy* 83, 595–619.
- Fukushima, T. (2012). Numerical computation of spherical harmonics of arbitrary degree and order by extending exponent of floating point numbers. *Journal of Geodesy* 86, 271–285.
- Gottlieb, R. G. (1993). *Fast gravity, gravity partials, normalized gravity, gravity gradient torque and magnetic field: derivation, code and data.* NASA Contractor Report 188243.
- Holmes, S. A., Featherstone, W. E. (2002). A unified approach to the Clenshaw summation and the recursive computation of very high degree and order normalised associated Legendre functions. *Journal of Geodesy* 76, 279–299.
- Lear, W. M. (1986). *The gravitational acceleration equations.* JSC Internal Note 86-FM-15 (JSC-22080), NASA Johnson Space Center.
- Metris, G., Xu, J., Wytrzyszczak, I. (1999). Derivatives of the gravity potential with respect to rectangular coordinates. *Celestial Mechanics and Dynamical Astronomy* 71, 137–151.
- Montenbruck, O., Gill, E. (2000). *Satellite Orbits: Models, Methods and Applications.* Springer.
- Pines, S. (1973). Uniform representation of the gravitational potential and its derivatives. *AIAA Journal* 11(11), 1508–1511.
