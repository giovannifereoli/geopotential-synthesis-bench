# Spherical-harmonic gravity synthesis: accuracy and CPU time of eight methods

This repo compares eight ways to evaluate a spherical-harmonic gravity model at scattered body-fixed points, such as the points of an orbit. Each method computes two outputs:

- **Acceleration:** the gravitational acceleration (3 components).
- **Acceleration + gravity gradient:** the acceleration plus the gravity-gradient tensor (6 more components). Orbit determination needs this second output for the variational equations.

All eight methods are scalar Rust kernels that follow one contract. Each is checked against a double-double reference and timed with the same protocol. Five of the methods come from Fantino & Casotto (2009), [Fantino_Casotto_2009_JGeod_Methods_of_harmonic_synthesis.pdf](Fantino_Casotto_2009_JGeod_Methods_of_harmonic_synthesis.pdf). The other three are Gottlieb (1993), Lear (1986) and Holmes & Featherstone (2002), the method used by Orekit. Lear's method gives only the acceleration; the gravity gradient for LEAR is our own extension (see [LEAR gravity gradient](#lear-gravity-gradient-our-extension)). For HOLMES we use the recursion of Holmes & Featherstone with our own pole-safe derivatives (see [HOLMES](#holmes-holmesfeatherstone)).

## Main findings

The time ratios below are medians over the three point sets (as in the figures); the times in µs are at the ORBIT1 points (as in the tables).

- **All methods are accurate.** Where they do not overflow, every method has a max. relative error between 2.2 × 10⁻¹³ and 3.1 × 10⁻¹², for both outputs, at every N from 10 to 2000 and on all three point sets, including points 10⁻⁸ degrees from the poles. The error grows slowly with N. No method is consistently more accurate than the others.
- **Three methods overflow at high degree in plain double precision.** `PINES` and `GOTTLIEB` overflow at N = 1500 and 2000, and `CLENSHAW` at N = 2000, at high-latitude points. See [why](#why-pines-gottlieb-and-clenshaw-overflow).
- **Acceleration only: `LEAR` is the fastest.** From N = 100 on, it needs 7–14% less CPU time than `LEGENDRE_a`, the second fastest (72.5 vs 83.9 µs per point at N = 360). At N ≤ 50 the two are tied within timing noise.
- **Acceleration + gravity gradient: `METRIS` is the fastest.** From N = 20 to 1000 it needs 12–22% less CPU time than `LEGENDRE_a` (131 vs 175 µs per point at N = 360). From N = 1500 on its lead shrinks (0.85 and 0.92 times `LEGENDRE_a`), and from N = 720 on it is often tied with `HOLMES`.
- **`HOLMES` is second with the gradient.** With the gradient it needs 2–10% less CPU time than `LEGENDRE_a` from N = 20 on (154 µs per point at N = 360). For the acceleration alone it takes 4–13% longer than `LEGENDRE_a` from N = 50 on. It recurs on the same functions as `PINES` and `GOTTLIEB` but does not overflow up to N = 2000, because it scales them by a fixed factor.
- **`LEAR` is the slowest once the gradient is needed.** With the gradient it takes 1.3–1.5 times as long as `LEGENDRE_a` from N = 50 on. Our gradient extension carries extra recursions (P̄_nm / cos² φ and second derivatives) so that it never divides by cos φ. `LEAR` is therefore fastest only when the gradient is not needed.
- **`CLENSHAW` is slow at scattered points.** For the acceleration it takes 2.0–2.5 times as long as `LEGENDRE_a` from N = 100 on; with the gradient, 1.04–1.14 times. Its speed advantage needs a grid, which is not tested here.
- **`PINES` and `GOTTLIEB`** take 11–25% longer than `LEGENDRE_a` for the acceleration, and 13–30% longer with the gradient. With the gradient, `GOTTLIEB` is 2–11% faster than `PINES`; for the acceleration alone the two are about equal.
- **`LEGENDRE_b`** (the whole Legendre triangle is stored first, then summed) takes 1–20% longer than `LEGENDRE_a`, most at N = 10. **`METRIS` for the acceleration alone** takes 2–22% longer than `LEGENDRE_a`, most at N = 10 and N = 2000.
- The ranking is the same on all three point sets.
- **The ratio curves are not monotone in N, but the times are.** Every method's CPU time per point grows steadily with N. The ratios to `LEGENDRE_a` move up and down by a few per cent from one N to the next. Part of this is real (fixed costs at small N, the cache at N ≥ 1500), part is timing noise of about ±5% per cell. See [Reading the CPU-time curves](#reading-the-cpu-time-curves).

## Results

The results below come from one run on 12 November 2026, on one performance core of an Apple M3 Pro (macOS 26.7, rustc 1.87.0, opt-level 3, LTO). The timing phase took 17.3 minutes. Other jobs were running on the machine (load average 3–6). The A/A probe differed from 1 by 1.2% (median) and 8.3% (maximum). Differences below about 5% are therefore not resolved, and a few cells are noisier; the winner rule takes this into account.

### Figures

**[results/fig_accuracy.pdf](results/fig_accuracy.pdf)** has one panel per method (eight panels). Each panel shows the max. relative error over all test points against the maximum degree N:

- blue solid line: the acceleration;
- orange dashed line: the gravity gradient;
- dotted line: the pass threshold 10⁻¹¹;
- grey ×: a degree at which the method overflowed.

All panels share the same axes, so the heights of the curves can be compared from one panel to the next. The legend is below the panels.

**[results/fig_speed.pdf](results/fig_speed.pdf)** shows the CPU time. Each curve is the median over the three point sets, and the shaded band is their range (lowest to highest).

- **Panel (a)** shows the CPU time per point of `LEGENDRE_a`, in µs: blue for the acceleration, orange for the acceleration + gravity gradient. The time grows roughly as N².
- **Panels (b)**, one per other method (seven panels), show its CPU time divided by that of `LEGENDRE_a`, for the same output and the same N. Below 1 means faster than `LEGENDRE_a`. The light grey band is ±5%, the size of difference the timing cannot resolve. A grey × marks an overflow; those cases are not timed.

To get a method's time in µs, multiply its ratio in (b) by the value in (a), or read the tables below.

Two more figures put all the methods on one axis, so they can be compared directly. Each has one panel per output, and each method has its own colour and marker.

- **[results/fig_speed_all.pdf](results/fig_speed_all.pdf)** shows every method's CPU time divided by that of `LEGENDRE_a` (the black line at 1): the median over the three point sets, with a vertical bar for their range. Each curve is labelled at its right end, and the grey band is again ±5%. Panel (a) is the acceleration; panel (b) is the acceleration + gravity gradient.
- **[results/fig_accuracy_all.pdf](results/fig_accuracy_all.pdf)** shows every method's max. relative error. Panel (a) is the acceleration; panel (b) is the gravity gradient. `LEGENDRE_a` and `LEGENDRE_b` give identical results, so they share one curve. The curves cross each other often, because no method is consistently more accurate than the others.

In both figures, a curve stops where its method overflows.

Run `python plots.py` to redraw all four figures from the saved tables, without rerunning anything.

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

Reference check (two independent double-double formulations): max difference 4e-28.

### CPU time

CPU time per point in microseconds (ORBIT1 points; median of the blocks, normalised by the calibration loop). **bold** = fastest; tied methods (difference not resolvable, see winners.csv) marked with *. x = failed the accuracy check, - = overflow (not timed).

**Acceleration**

| Method | N=10 | N=20 | N=50 | N=100 | N=200 | N=360 | N=720 | N=1000 | N=1500 | N=2000 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| LEGENDRE_a | **0.148** | **0.413** | 1.88* | 6.85* | 26.2 | 83.9 | 334 | 651 | 1561 | 2786 |
| LEGENDRE_b | 0.174 | 0.459 | 2.02* | 7.28* | 28.3 | 87.5 | 354 | 685 | 1723 | 3160 |
| PINES | 0.184 | 0.509 | 2.08 | 8.05 | 30.1 | 95.6 | 383 | 749 | - | - |
| GOTTLIEB | 0.175 | 0.491 | 2.13 | 8.10 | 29.3 | 93.9 | 371 | 722 | - | - |
| METRIS | 0.180 | 0.490 | 2.07 | 7.40* | 27.5 | 86.8 | 357 | 696 | 1745 | 3457 |
| CLENSHAW | 0.163* | 0.525 | 3.35 | 13.9 | 59.4 | 197 | 816 | 1557 | 3536 | - |
| LEAR | 0.151* | 0.434* | **1.84** | **6.56** | **23.4** | **72.5** | **289** | **559** | **1374** | **2455** |
| HOLMES | 0.180 | 0.473 | 2.07 | 7.77 | 28.9 | 90.4 | 363 | 686 | 1617 | 2890 |

**Acceleration + gravity gradient**

| Method | N=10 | N=20 | N=50 | N=100 | N=200 | N=360 | N=720 | N=1000 | N=1500 | N=2000 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| LEGENDRE_a | 0.261* | 0.781 | 3.77 | 13.5 | 52.3 | 175 | 730 | 1346 | 3026 | 5471 |
| LEGENDRE_b | 0.292 | 0.817 | 3.81 | 14.3 | 54.6 | 184 | 687 | 1357 | 3315 | 5631 |
| PINES | 0.338 | 0.905 | 4.42 | 16.9 | 63.7 | 207 | 849 | 1665 | - | - |
| GOTTLIEB | 0.294 | 0.891 | 4.34 | 16.2 | 62.3 | 200 | 809 | 1531 | - | - |
| METRIS | **0.249** | **0.634** | **3.25** | **11.0** | **41.0** | **131** | **537** | **1068** | **2568** | **4947** |
| CLENSHAW | 0.275* | 0.807 | 4.12 | 15.5 | 58.7 | 192 | 747 | 1442 | 3303 | - |
| LEAR | 0.274* | 0.898 | 4.95 | 19.4 | 73.6 | 232 | 922 | 1803 | 4131 | 7956 |
| HOLMES | 0.278* | 0.763 | 3.55 | 13.2 | 48.7 | 154 | 627* | 1226* | 2738* | 4958* |

### Fastest method per N (all point sets)

| Output | Points | N=10 | N=20 | N=50 | N=100 | N=200 | N=360 | N=720 | N=1000 | N=1500 | N=2000 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Acc. | ORBIT1 | tie: LEGENDRE_a+LEAR+CLENSHAW | tie: LEGENDRE_a+LEAR | tie: LEAR+LEGENDRE_a+LEGENDRE_b | tie: LEAR+LEGENDRE_a+LEGENDRE_b+METRIS | LEAR | LEAR | LEAR | LEAR | LEAR | LEAR |
| Acc. | UNIFORM | tie: LEGENDRE_a+LEAR+CLENSHAW | tie: LEGENDRE_a+LEAR | tie: LEAR+LEGENDRE_a | LEAR | LEAR | LEAR | LEAR | LEAR | tie: LEAR+LEGENDRE_a | LEAR |
| Acc. | POLAR | tie: LEAR+LEGENDRE_a+CLENSHAW | tie: LEGENDRE_a+LEAR | tie: LEAR+LEGENDRE_a+LEGENDRE_b | LEAR | LEAR | LEAR | LEAR | LEAR | LEAR | LEAR |
| Acc. + grad. | ORBIT1 | tie: METRIS+LEGENDRE_a+LEAR+CLENSHAW+HOLMES | METRIS | METRIS | METRIS | METRIS | METRIS | tie: METRIS+HOLMES | tie: METRIS+HOLMES | tie: METRIS+HOLMES | tie: METRIS+HOLMES |
| Acc. + grad. | UNIFORM | tie: METRIS+LEGENDRE_a+CLENSHAW+HOLMES | METRIS | tie: METRIS+HOLMES | METRIS | METRIS | METRIS | METRIS | METRIS | tie: METRIS+HOLMES | tie: HOLMES+METRIS |
| Acc. + grad. | POLAR | tie: METRIS+LEGENDRE_a+HOLMES+CLENSHAW+LEAR | METRIS | tie: METRIS+HOLMES | METRIS | METRIS | METRIS | METRIS | tie: METRIS+HOLMES | tie: METRIS+HOLMES | tie: HOLMES+METRIS |

A/A probe (LEGENDRE_a timed against an identical second copy of itself): median |ratio - 1| = 1.2%, max 8.3%; a win needs the 95% interval of the time ratio above 1 + max(5%, 2 |A/A - 1|).

Point-set check (the cost of a point does not depend on where it is, so the 3 point sets are repeat measurements): spread (max - min) / median of the time ratio to LEGENDRE_a over the sets = 3.7% (median), 14.4% (max, LEAR level 2 N=2000).

### Reading the CPU-time curves

The CPU time per point grows steadily with N for every method (panel (a) of `fig_speed`, and the tables above). The ratios to `LEGENDRE_a` do not: they move up and down from one N to the next. Four things cause this.

- **Fixed costs at small N.** At N ≤ 50 a good part of the time goes to work done once per point (trigonometric functions, starting values of the recursions), and this work differs between methods. As N grows, the work per coefficient takes over. This is why, for the acceleration, `CLENSHAW` starts near 1 and climbs to 2.4, while `HOLMES`, `LEGENDRE_b` and `METRIS` start at about 1.2 and fall.
- **The cache at large N.** The coefficients C_nm and S_nm that a method reads at every point take about 8 N² bytes: 8 MB at N = 1000, 18 MB at N = 1500. The L2 cache of the performance cores holds 16 MB, so from N = 1500 on they no longer fit and are read from main memory at every point. `METRIS` reads more: its derivative operators are applied to the coefficients at setup, so it reads 3 coefficient sets for the acceleration and 9 with the gradient (about 24 N² and 72 N² bytes). These leave the cache from N ≈ 800 and N ≈ 470. Its ratio to `LEGENDRE_a` climbs steadily at large N (acceleration from 1.04 at N = 360 to 1.22 at N = 2000; with the gradient from 0.78 to 0.92), most likely because of this extra memory traffic. `LEGENDRE_b`, which also stores the whole Legendre triangle, rises at N = 1500 (acceleration 1.08 → 1.15).
- **`LEAR` with the gradient at N = 2000.** At the POLAR points it takes 10–13% longer than at the other two point sets. This is the largest point-set spread in the run (14.4% in the ratio to `LEGENDRE_a`), and an earlier run showed it too. The kernel has no data-dependent branches and the slowdown is not caused by subnormal numbers; we have not found the cause.
- **Timing noise.** Other jobs on the machine share the L2 cache and change the clock speed. The calibration loop removes most of the clock-speed changes: the three point sets, which do identical work, disagree by 3.7% (median) after calibration, against 10% without it. It cannot see cache sharing, so each cell keeps a noise of about ±5% (A/A probe up to 8.3%).

To read the curves:

- Treat a ratio within ±5% of 1 as a tie with `LEGENDRE_a`.
- Treat a step of less than about 5% between neighbouring N as noise, unless it is wider than the range band and continues over several N.
- Use the median over the three point sets (the curves in the figures) rather than one point set. The tables give ORBIT1 only, so they show more of the noise.

### Why `PINES`, `GOTTLIEB` and `CLENSHAW` overflow

`PINES` and `GOTTLIEB` recur on the Helmholtz polynomials Ā_nm(u) = P̄_nm(sin φ) / cos^m φ, with u = sin φ. Near the poles these polynomials grow very fast with the degree. Their largest values, found numerically, are:

| N | Largest log10 \|Ā_nm\| | Where it passes the double-precision limit (about 10³⁰⁸) |
|---|---|---|
| 1000 | 209, at the pole | Nowhere |
| 1500 | 314 at the pole; 308 at u = 0.99 | Latitudes above about 82° |
| 2000 | 418 at the pole; 343 at u = 0.9; 274 at u = 0.8 | Latitudes above about 58° |

All three point sets have points above 82° latitude (ORBIT1 is a near-polar orbit), so both methods fail in every set at N = 1500 and 2000.

`CLENSHAW` forms, for each order m, the sum over n of C_nm (a/r)^(n−m) P̄_nm / P̄_mm, and multiplies it by the sectoral term (a/r)^m P̄_mm only at the end. The ratio P̄_nm / P̄_mm grows like Ā_nm. Inside the sum it is multiplied by the small coefficients and by (a/r)^(n−m), which delays the overflow: in these tests the sum stays finite up to N = 1500 and overflows at N = 2000.

The other five methods do not overflow up to N = 2000:

- `LEGENDRE_a`, `LEGENDRE_b` and `LEAR` keep the factor cos^m φ inside their Legendre values (they divide by at most cos² φ).
- `METRIS` works with solid harmonics, which include the factor (a/r)^(n+1).
- `HOLMES` recurs on the same function as `PINES`, P̄_nm / cos^m φ, but every value is multiplied by the fixed factor 2⁻⁹³⁰ ≈ 10⁻²⁸⁰. The largest value at N = 2000, 10⁴¹⁸ near the pole, is therefore stored as about 10¹³⁸. The factor is removed once, at the end.

A fixed scale factor, as in `HOLMES`, is the simplest cure. Scaled numbers with their own exponent, such as the X-numbers of Fukushima (2012), work at any degree but cost more time. No kernel here uses X-numbers.

### Differences from Fantino & Casotto (2009)

The paper timed complete programmes on 2009 hardware, including file input and output, and its `LEGENDRE_b` ran through Matlab. Here only the synthesis kernels are timed, and all of them follow the same contract. The main differences:

- **`CLENSHAW`:** the paper found it about 2 times faster than `LEGENDRE_a` at scattered points. Here it is slower: 2.0–2.5 times for the acceleration, and 1.04–1.14 times with the gradient.
- **The other methods:** in the paper's Table 7, `LEGENDRE_b`, `PINES` and `METRIS` take 4.8, 3.0 and 1.6 times as long as `LEGENDRE_a`. Here, for the acceleration + gravity gradient at N = 360, they take 1.04, 1.22 and 0.78 times as long.
- **New methods:** `GOTTLIEB`, `LEAR` and `HOLMES` are not in the paper.

These differences most likely come from the 2009 implementations and hardware rather than from the algorithms.

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

**Pole safety.** None of the methods divides by cos φ:

- The spherical methods `LEGENDRE_a`, `LEGENDRE_b` and `CLENSHAW` recur on Z_nm = P̄_nm / cos(φ)^min(m,2). With this choice the 1/cos φ and 1/cos² φ factors of the local-frame formulas cancel analytically. The cost is about two extra multiplications per term.
- `LEAR` works with P̄_nm / cos φ for the acceleration and adds P̄_nm / cos² φ for the gravity gradient.
- `HOLMES` works with P̄_nm / cos^m φ and multiplies the factor cos^m φ back in by Horner's rule. The derivative terms that would need cos^(m−1) φ or cos^(m−2) φ have their own Horner chains.
- `PINES`, `GOTTLIEB` and `METRIS` are Cartesian and have no singularity.

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
- **Cells:** a cell is one combination of degree N, output and point set. Each method is timed in its own fresh process. Cells run in random order, and so do the methods within a cell, so that slow drift in the machine cannot line up with N or with a method.
- **Rounds:** the whole list of cells is run 3 times, each time in a new random order. Each (cell, method) is therefore measured in 3 separate processes spread over the whole run, and its blocks are pooled. This spreads slow changes in the machine's load over all cells and methods, so they cannot line up with one cell or one method. It does not remove them. Other jobs that share the 16 MB L2 cache of the performance cores slow the large-N kernels for seconds to minutes at a time, and the calibration loop (which runs in registers) cannot see this. What remains is about ±5% per cell; see [Reading the CPU-time curves](#reading-the-cpu-time-curves).
- **Inside a process:**
  - One warm-up block, which is discarded.
  - 4 timed blocks of at least 100 ms each (12 per cell and method over the 3 rounds). A block is a number of passes over K points, with K = 4·10⁶ / (number of coefficients), clamped to 32–4096 points.
  - A fixed calibration loop runs before and after every block: 8 independent multiply-add chains with no FMA.
  - Each block time is multiplied by c_ref / c_block. Here c_block is the mean of the two calibration runs around that block, and c_ref is the 5th percentile of all calibration runs in the whole timing phase (the machine at full speed). This removes changes in clock speed between blocks and between processes. It matters: without it, the three point sets disagree about 3 times more, 10% instead of 3.7% (median; see [Reading the CPU-time curves](#reading-the-cpu-time-curves)).
  - All blocks are kept. An earlier version dropped blocks whose two calibration runs differed by more than 3%; this made the worst cases slightly worse, not better.
- **Reported time:** the median over blocks, in µs per point, with a 95% bootstrap interval.
- **Point-set check:** the kernels have no data-dependent branches, so the cost of a point does not depend on where it is (one exception, below). The three point sets are therefore three repeat measurements of each (N, output, method). The figures show their median, with their range as a band; `summary.md` reports how far they disagree.
- **A/A probe:** a second compiled copy of `LEGENDRE_a` is timed as if it were a ninth method. Any difference it shows from `LEGENDRE_a` is measurement noise.
- **Winner rule:** the fastest passing method wins outright only if every other method is slower by more than δ, with δ = max(5%, 2 × |A/A ratio − 1|). "Slower by more than δ" means that the lower end of the 95% bootstrap interval of the time ratio is above 1 + δ. Otherwise the methods within δ are reported as tied.

## What operational software uses

This is from public documentation and source code, as of November 2026.

| Tool | Formulation | Closest method here |
|---|---|---|
| Scarabaeus (CCAR-ORCCA) | Normalised Cunningham solid-harmonic recursion, B_nm = V_nm + iW_nm (Montenbruck & Gill 2000, Sect. 3.2), with an analytic Jacobian. Its docstrings say "Pines/Eckman", but the code is the Cunningham recursion. | `METRIS` |
| Basilisk | Normalised Pines | `PINES` |
| GMAT | Normalised Pines, with an optional gradient for the state transition matrix | `PINES` |
| Nyx | Pines, ported from GMAT | `PINES` |
| Tudat | Fully normalised Legendre functions in spherical coordinates, then a rotation | `LEGENDRE_a` |
| Orekit | Holmes–Featherstone (modified forward row method) | `HOLMES` (with pole-safe derivatives) |
| MONTE (JPL) | Not public | – |

None of these tools uses Gottlieb's or Lear's formulation. Both come from NASA Johnson Space Center; Eckman et al. (2016) normalise and compare them with Pines.

## Running it

**Requirements**

- Python 3 with NumPy and Matplotlib.
- A Rust toolchain (`cargo`). The crate has no dependencies.

`benchmark.py` builds the Rust binary on first use, and again whenever a source file is newer than the binary. The build goes to `~/.cache/sph_rs/target`; set `SPH_RS_TARGET` to use another folder.

**Commands**

```sh
python benchmark.py                   # accuracy check + timing + analysis -> results/
python benchmark.py --quick           # smoke test: N = 10, 50, 360, short blocks -> results_quick/
python benchmark.py --phase gate      # accuracy check only
python benchmark.py --phase time      # timing only (reuses results/accuracy.csv)
python benchmark.py --phase analyze   # tables and figures from the saved raw data
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
| `--out` | `results/` | Output folder |

Every run saves its raw data in the output folder. `--phase analyze` rebuilds the tables and figures from `raw_timing.jsonl.gz` and `accuracy.csv`, and `plots.py` redraws the figures from the tables, so neither needs a new run or a Rust build. Use `--out` to keep several runs side by side.

The accuracy check takes about 1 minute. The timing phase starts about 1530 processes (60 cells × 3 rounds × up to 9 methods, including the A/A copy) and takes about 17 minutes. For clean timings, run on mains power, keep the machine awake (on macOS: `caffeinate -i -s python benchmark.py`) and close other heavy jobs. The A/A probe line at the end of `summary.md` shows how noisy the run was.

The binary (`~/.cache/sph_rs/target/release/sph_rs`) can also be called directly:

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
| [plots.py](plots.py) | The four figures, drawn from `accuracy.csv` and `timing.csv` |
| [sph_rs/src/main.rs](sph_rs/src/main.rs) | Command line, model loading, recursion tables, `METRIS` coefficient operators |
| [sph_rs/src/bench.rs](sph_rs/src/bench.rs) | The eight kernels, the `eval` and `bench` commands, the timing loop |
| [sph_rs/src/reference.rs](sph_rs/src/reference.rs) | Double-double reference with X-numbers (two formulations) |
| [results/accuracy.csv](results/accuracy.csv) | Max. relative error per N, point set and method; pass or fail |
| [results/timing.csv](results/timing.csv) | CPU time per point (µs) with 95% interval, per N, output, point set and method |
| [results/winners.csv](results/winners.csv) | Fastest method per cell, ties, A/A ratio |
| [results/summary.md](results/summary.md) | The tables shown above |
| [results/fig_accuracy.pdf](results/fig_accuracy.pdf), [results/fig_speed.pdf](results/fig_speed.pdf) | Figures, one panel per method |
| [results/fig_accuracy_all.pdf](results/fig_accuracy_all.pdf), [results/fig_speed_all.pdf](results/fig_speed_all.pdf) | Figures, all methods on one axis |
| [results/raw_timing.jsonl.gz](results/raw_timing.jsonl.gz) | Every timed block, for re-analysis |
| [results/meta.json](results/meta.json) | Machine, compiler, run settings |
| [Fantino_Casotto_2009_JGeod_Methods_of_harmonic_synthesis.pdf](Fantino_Casotto_2009_JGeod_Methods_of_harmonic_synthesis.pdf) | The paper the comparison started from |

## Limits

- **One machine.** All times are from one Apple M3 Pro performance core. Ratios between methods may differ on other processors, especially between the Cartesian and the spherical methods.
- **Scalar code at scattered points.** Each point is evaluated on its own, on one thread, without SIMD across points. Grids, where Clenshaw summation and per-parallel reuse pay off, are not tested.
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
