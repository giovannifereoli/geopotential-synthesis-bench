# Spherical-harmonic gravity synthesis: accuracy and CPU time of seven methods

This repo compares seven ways to evaluate a spherical-harmonic gravity model at scattered body-fixed points, such as the points of an orbit. Each method computes two outputs:

- **Acceleration:** the gravitational acceleration (3 components).
- **Acceleration + gravity gradient:** the acceleration plus the gravity-gradient tensor (6 more components). Orbit determination needs this second output for the variational equations.

All seven methods are scalar Rust kernels that follow one contract. Each is checked against a double-double reference and timed with the same protocol. Five of the methods come from Fantino & Casotto (2009), [Fantino_Casotto_2009_JGeod_Methods_of_harmonic_synthesis.pdf](Fantino_Casotto_2009_JGeod_Methods_of_harmonic_synthesis.pdf). The other two are Gottlieb (1993) and Lear (1986). Lear's method gives only the acceleration; the gravity gradient for LEAR is our own extension (see [LEAR gravity gradient](#lear-gravity-gradient-our-extension)).

## Main findings

- **All methods are accurate.** Where they do not overflow, every method has a max. relative error between 2.2 × 10⁻¹³ and 3.1 × 10⁻¹², for both outputs, at every N from 10 to 2000 and on all three point sets, including points 10⁻⁸ degrees from the poles. The error grows slowly with N. No method is consistently more accurate than the others.
- **Three methods overflow at high degree in plain double precision.** `PINES` and `GOTTLIEB` overflow at N = 1500 and 2000, and `CLENSHAW` at N = 2000, at high-latitude points. See [why](#why-pines-gottlieb-and-clenshaw-overflow).
- **Acceleration only: `LEAR` is the fastest.** From N = 100 on, it needs 7–15% less CPU time than `LEGENDRE_a`, the second fastest (74.9 vs 85.8 µs per point at N = 360). At N ≤ 50 the two are tied within timing noise.
- **Acceleration + gravity gradient: `METRIS` is the fastest** at every N. From N = 100 to 1000 it needs 13–22% less CPU time than `LEGENDRE_a`, the second fastest (135 vs 173 µs per point at N = 360). At N ≤ 50 and N ≥ 1500 the margin is smaller (2–17%), and a few of those cases are ties.
- **`LEAR` is the slowest once the gradient is needed.** With the gradient it takes 1.3–1.6 times as long as `LEGENDRE_a` from N = 50 on. Our gradient extension carries extra recursions (P̄_nm / cos² φ and second derivatives) so that it never divides by cos φ. `LEAR` is therefore fastest only when the gradient is not needed.
- **`CLENSHAW` is slow at scattered points.** For the acceleration it takes 2.0–2.5 times as long as `LEGENDRE_a` from N = 100 on; with the gradient, 1.04–1.16 times. Its speed advantage needs a grid, which is not tested here.
- **`PINES` and `GOTTLIEB`** take 9–20% longer than `LEGENDRE_a` for the acceleration, and 10–29% longer with the gradient. With the gradient, `GOTTLIEB` is 2–8% faster than `PINES`; for the acceleration alone the two are about equal.
- **`LEGENDRE_b`** (the whole Legendre triangle is stored first, then summed) takes up to 20% longer than `LEGENDRE_a`, and about 5% longer at N = 360. **`METRIS` for the acceleration alone** takes 3–23% longer than `LEGENDRE_a`.
- The ranking is the same on all three point sets.

## Results

The results below come from one run on 12 November 2026, on one performance core of an Apple M3 Pro (macOS 26.7, rustc 1.87.0, opt-level 3, LTO). The timing phase took 12.5 minutes. The A/A probe differed from 1 by 0.6% (median) and 3.0% (maximum), so time differences above about 5% are real.

### Figures

**[results/fig_accuracy.pdf](results/fig_accuracy.pdf)** has one panel per method. Each panel shows the max. relative error against the maximum degree N:

- blue solid line: the acceleration;
- orange dashed line: the gravity gradient;
- dotted line: the pass threshold 10⁻¹¹;
- grey ×: a degree at which the method overflowed.

All panels share the same axes, so the heights of the curves can be compared from one panel to the next.

**[results/fig_speed.pdf](results/fig_speed.pdf)** shows the CPU time at the ORBIT1 points:

- **Panel (a)** shows the CPU time per point of `LEGENDRE_a`, in µs: blue for the acceleration, orange for the acceleration + gravity gradient. The time grows roughly as N².
- **Panels (b)** show, for each other method, its CPU time divided by that of `LEGENDRE_a`, for the same output and the same N. Below 1 means faster than `LEGENDRE_a`. The grey band is ±5%, the size of difference the timing cannot resolve. A grey × marks an overflow; those cases are not timed.

To get a method's time in µs, multiply its ratio in (b) by the value in (a), or read the tables below.

Two more figures put all the methods on one axis, so they can be compared directly. Each has one panel per output, and each method has its own colour and marker.

- **[results/fig_speed_all.pdf](results/fig_speed_all.pdf)** shows every method's CPU time divided by that of `LEGENDRE_a` (the black line at 1), at the ORBIT1 points. Each curve is labelled at its right end, and the grey band is again ±5%. Panel (a) is the acceleration; panel (b) is the acceleration + gravity gradient.
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

Reference check (two independent double-double formulations): max difference 4e-28.

### CPU time

CPU time per point in microseconds (ORBIT1 points; median of the blocks, normalised by the calibration loop). **bold** = fastest; tied methods (difference not resolvable, see winners.csv) marked with *. x = failed the accuracy check, - = not run (overflow).

**Acceleration**

| Method | N=10 | N=20 | N=50 | N=100 | N=200 | N=360 | N=720 | N=1000 | N=1500 | N=2000 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| LEGENDRE_a | **0.153** | **0.430** | 1.95* | 7.06 | 27.1 | 85.8 | 348 | 667 | 1663 | 2919 |
| LEGENDRE_b | 0.181 | 0.479 | 2.04 | 7.75 | 29.1 | 91.0 | 367 | 732 | 1903 | 3217 |
| PINES | 0.182 | 0.487 | 2.17 | 8.18 | 30.9 | 98.2 | 402 | 794 | - | - |
| GOTTLIEB | 0.179 | 0.477 | 2.12 | 8.40 | 30.5 | 96.1 | 387 | 774 | - | - |
| METRIS | 0.188 | 0.487 | 2.20 | 7.63 | 28.0 | 90.2 | 370 | 733 | 1826 | 3504 |
| CLENSHAW | 0.160* | 0.545 | 3.39 | 14.5 | 60.5 | 199 | 825 | 1611 | 3752 | - |
| LEAR | 0.154* | 0.437* | **1.88** | **6.57** | **24.0** | **74.9** | **299** | **571** | **1452** | **2574** |

**Acceleration + gravity gradient**

| Method | N=10 | N=20 | N=50 | N=100 | N=200 | N=360 | N=720 | N=1000 | N=1500 | N=2000 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| LEGENDRE_a | 0.274 | 0.780 | 3.82 | 14.0 | 54.2 | 173 | 694 | 1362 | 3234* | 5664 |
| LEGENDRE_b | 0.296 | 0.832 | 3.92 | 14.6 | 56.2 | 181 | 720 | 1389 | 3359* | 5993 |
| PINES | 0.326 | 0.928 | 4.57 | 17.4 | 66.1 | 213 | 879 | 1727 | - | - |
| GOTTLIEB | 0.304 | 0.893 | 4.47 | 16.6 | 64.4 | 206 | 823 | 1596 | - | - |
| METRIS | **0.257** | **0.651** | **3.35** | **11.3** | **42.5** | **135** | **565** | **1132** | **3182** | **5235** |
| CLENSHAW | 0.284 | 0.833 | 4.20 | 15.8 | 60.6 | 192 | 770 | 1501 | 3737 | - |
| LEAR | 0.288 | 0.911 | 4.96 | 19.2 | 74.6 | 239 | 965 | 1876 | 4389 | 8207 |

### Fastest method per N (all point sets)

| Output | Points | N=10 | N=20 | N=50 | N=100 | N=200 | N=360 | N=720 | N=1000 | N=1500 | N=2000 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Acc. | ORBIT1 | tie: LEGENDRE_a+LEAR+CLENSHAW | tie: LEGENDRE_a+LEAR | tie: LEAR+LEGENDRE_a | LEAR | LEAR | LEAR | LEAR | LEAR | LEAR | LEAR |
| Acc. | UNIFORM | tie: LEGENDRE_a+LEAR+CLENSHAW | tie: LEGENDRE_a+LEAR | tie: LEAR+LEGENDRE_a | LEAR | LEAR | LEAR | LEAR | LEAR | LEAR | LEAR |
| Acc. | POLAR | tie: LEGENDRE_a+LEAR+CLENSHAW | tie: LEGENDRE_a+LEAR | tie: LEAR+LEGENDRE_a | LEAR | LEAR | LEAR | LEAR | LEAR | LEAR | LEAR |
| Acc. + grad. | ORBIT1 | METRIS | METRIS | METRIS | METRIS | METRIS | METRIS | METRIS | METRIS | tie: METRIS+LEGENDRE_a+LEGENDRE_b | METRIS |
| Acc. + grad. | UNIFORM | tie: METRIS+LEGENDRE_a | METRIS | METRIS | METRIS | METRIS | METRIS | METRIS | METRIS | METRIS | METRIS |
| Acc. + grad. | POLAR | METRIS | METRIS | tie: METRIS+LEGENDRE_a | METRIS | METRIS | METRIS | METRIS | METRIS | METRIS | tie: METRIS+LEGENDRE_a |

A/A probe (LEGENDRE_a timed against an identical second copy of itself): median |ratio - 1| = 0.6%, max 3.0%; a win needs the 95% interval of the time ratio above 1 + max(5%, 2 |A/A - 1|).

### Why `PINES`, `GOTTLIEB` and `CLENSHAW` overflow

`PINES` and `GOTTLIEB` recur on the Helmholtz polynomials Ā_nm(u) = P̄_nm(sin φ) / cos^m φ, with u = sin φ. Near the poles these polynomials grow very fast with the degree. Their largest values, found numerically, are:

| N | Largest log10 \|Ā_nm\| | Where it passes the double-precision limit (about 10³⁰⁸) |
|---|---|---|
| 1000 | 209, at the pole | Nowhere |
| 1500 | 314 at the pole; 308 at u = 0.99 | Latitudes above about 82° |
| 2000 | 418 at the pole; 343 at u = 0.9; 274 at u = 0.8 | Latitudes above about 58° |

All three point sets have points above 82° latitude (ORBIT1 is a near-polar orbit), so both methods fail in every set at N = 1500 and 2000.

`CLENSHAW` forms, for each order m, the sum over n of C_nm (a/r)^(n−m) P̄_nm / P̄_mm, and multiplies it by the sectoral term (a/r)^m P̄_mm only at the end. The ratio P̄_nm / P̄_mm grows like Ā_nm. Inside the sum it is multiplied by the small coefficients and by (a/r)^(n−m), which delays the overflow: in these tests the sum stays finite up to N = 1500 and overflows at N = 2000.

The other four methods do not overflow up to N = 2000:

- `LEGENDRE_a`, `LEGENDRE_b` and `LEAR` keep the factor cos^m φ inside their Legendre values (they divide by at most cos² φ).
- `METRIS` works with solid harmonics, which include the factor (a/r)^(n+1).

Scaled numbers, such as the X-numbers of Fukushima (2012), would remove the overflow at some cost in speed. No kernel here uses them.

### Differences from Fantino & Casotto (2009)

The paper timed complete programmes on 2009 hardware, including file input and output, and its `LEGENDRE_b` ran through Matlab. Here only the synthesis kernels are timed, and all of them follow the same contract. The main differences:

- **`CLENSHAW`:** the paper found it about 2 times faster than `LEGENDRE_a` at scattered points. Here it is slower: 2.0–2.5 times for the acceleration, and 1.04–1.16 times with the gradient.
- **The other methods:** in the paper's Table 7, `LEGENDRE_b`, `PINES` and `METRIS` take 4.8, 3.0 and 1.6 times as long as `LEGENDRE_a`. Here, for the acceleration + gravity gradient at N = 360, they take 1.05, 1.24 and 0.78 times as long.
- **New methods:** `GOTTLIEB` and `LEAR` are not in the paper.

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

**Pole safety.** None of the methods divides by cos φ:

- The spherical methods `LEGENDRE_a`, `LEGENDRE_b` and `CLENSHAW` recur on Z_nm = P̄_nm / cos(φ)^min(m,2). With this choice the 1/cos φ and 1/cos² φ factors of the local-frame formulas cancel analytically. The cost is about two extra multiplications per term.
- `LEAR` works with P̄_nm / cos φ for the acceleration and adds P̄_nm / cos² φ for the gravity gradient.
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
- **Rounds:** the whole list of cells is run 3 times, each time in a new random order. Each (cell, method) is therefore measured in 3 separate processes spread over the whole run, and its blocks are pooled. A single round is not enough on a shared machine: other jobs that share the 16 MB L2 cache of the performance cores slow the large-N kernels for seconds to minutes at a time, and the calibration loop (which runs in registers) cannot see this. With one process per method, such an episode biased single cells by up to 13% in an earlier run and made the ratio curves zig-zag.
- **Inside a process:**
  - One warm-up block, which is discarded.
  - 4 timed blocks of at least 100 ms each (12 per cell and method over the 3 rounds). A block is a number of passes over K points, with K = 4·10⁶ / (number of coefficients), clamped to 32–4096 points.
  - A fixed calibration loop runs before and after every block: 8 independent multiply-add chains with no FMA. A block is dropped if the two calibration runs differ by more than 3%. If fewer than 5 blocks are left, all blocks are kept.
  - Each block time is multiplied by c_ref / c_block. Here c_block is the mean of the two calibration runs around that block, and c_ref is the 5th percentile of all calibration runs in the whole timing phase (the machine at full speed). This removes changes in clock speed between blocks and between processes.
- **Reported time:** the median over blocks, in µs per point, with a 95% bootstrap interval.
- **Point-set check:** the kernels have no data-dependent branches, so the cost of a point does not depend on where it is (one exception, below). The three point sets are therefore three repeat measurements of each (N, output, method). The figures show their median, with their range as a band; `summary.md` reports how far they disagree.
- **A/A probe:** a second compiled copy of `LEGENDRE_a` is timed as if it were an eighth method. Any difference it shows from `LEGENDRE_a` is measurement noise.
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
| Orekit | Holmes–Featherstone (modified forward row method) | – |
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

The accuracy check takes about 1 minute. The timing phase starts about 450 processes and takes about 13 minutes. For clean timings, run on mains power, keep the machine awake (on macOS: `caffeinate -i -s python benchmark.py`) and close other heavy jobs. The A/A probe line at the end of `summary.md` shows how noisy the run was.

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
| [sph_rs/src/bench.rs](sph_rs/src/bench.rs) | The seven kernels, the `eval` and `bench` commands, the timing loop |
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
- **Plain double precision.** No kernel uses X-numbers or scaling. That is why `PINES`, `GOTTLIEB` and `CLENSHAW` overflow at high degree. Degrees above 2000 are not tested.
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
