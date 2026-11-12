# geopotential-synthesis-bench

This repo benchmarks spherical-harmonic synthesis algorithms for the geopotential V, its gradient ∇V and its gradient tensor ∇∇V.

It re-implements the five programmes compared by Fantino & Casotto (2009) and adds Gottlieb's (1993) formulation. It then reproduces the paper's three studies at degree and order N = 360:

- **efficiency:** Table 7 and Fig. 4;
- **precision:** Figs. 5–10;
- **accuracy:** Eqs. 31–33 and Fig. 3.

It also adds a conditioning analysis of the Cartesian → local-frame rotation.

Every algorithm runs on two backends:

- **NumPy:** vectorised over the points. It works in `float64`, and an 80-bit `longdouble` version serves as the extended-precision reference.
- **Rust:** compiled scalar kernels that evaluate one point at a time, as the paper's Fortran did. Use this backend to compare speed.

The paper is in [Fantino_Casotto_2009_JGeod_Methods_of_harmonic_synthesis.pdf](Fantino_Casotto_2009_JGeod_Methods_of_harmonic_synthesis.pdf).

## Methods

| Programme | Functions | Coordinates | Origin |
|---|---|---|---|
| `LEGENDRE_a` | Fully normalised ALFs, forward column recursion (FCR) | spherical (λ, φ, r) | paper (Fortran90) |
| `LEGENDRE_b` | Same algorithm, point by point, with whole-triangle array operations | spherical | paper (Matlab + DLL) |
| `PINES` | Normalised Helmholtz polynomials Ā_nm(u), FCR | Cartesian direction cosines (s, t, u) | Pines (1973) |
| `GOTTLIEB` | Same Ā_nm(u), Gottlieb's derivative relations | Cartesian (s, t, u) | Gottlieb (1993). **Not in the paper.** |
| `METRIS` | Fully normalised solid harmonics (Cunningham); derivatives applied to the coefficients | Cartesian (x, y, z) | Metris et al. (1999) |
| `CLENSHAW` | Clenshaw summation over the degree on fully normalised ALFs | spherical; scattered points or per parallel on a grid | paper |

Every programme returns the paper's ten local-frame quantities, in the north-oriented east/north/radial frame:

V, V_λ, V_φ, V_r, V_λλ, V_λφ, V_φφ, V_λr, V_φr, V_rr

`PINES`, `GOTTLIEB` and `METRIS` also return Cartesian ∇V and ∇∇V natively.

**How GOTTLIEB differs from PINES.** Both use the same Helmholtz polynomials and the same normalised recursion. They differ in how they take derivatives:

- PINES needs the rows n+1 and n+2 of Ā.
- GOTTLIEB uses ∂_u Ā_nm = k_nm Ā_n,m+1 and ∂²_u Ā_nm = j_nm Ā_n,m+2, so it needs only row n.
- GOTTLIEB accumulates 14 sums per point and assembles ∇V and ∇∇V from them in closed form; see `_gottlieb_block` in `sph_core.py`.
- The relations follow NASA CR-188243 (Ada routine `Gotpot`). The implementation is checked numerically against PINES and CLENSHAW*16.

## Requirements

- **Python 3 with NumPy and Matplotlib.** Tested with Python 3.13 and NumPy 2.5. LaTeX is optional: without it, the figures use mathtext.
- **An 80-bit `np.longdouble` for the `*16` reference.** That means x86-64 Linux or macOS. On Apple silicon, use an x86-64 Python under Rosetta. Where `longdouble` is plain double, `sph_precision.py` prints a warning, and the reference is no better than the methods.
- **A Rust toolchain (`cargo`) for `--backend rust`.** Tested with rustc 1.87. The crate needs no external dependencies.
  - It is built on first use, into `~/.cache/sph_rs/target`. Set `SPH_RS_TARGET` to change that location.
  - The kernels use IEEE double, with no fast-math, no FMA contraction and no `target-cpu=native`.

## Running

```bash
python sph_core.py                  # self-test: six NumPy programmes vs CLENSHAW*16 at N = 30, Laplace, Cartesian cross-checks
python sph_rust.py                  # same for the Rust kernels (builds the crate)

python run_all.py --quick           # smoke test of every experiment (~2 min) -> ./quick_results
python run_all.py                   # NumPy backend                           -> ./results
python run_all.py --backend rust    # Rust backend                            -> ./results_rust
python run_all.py --gfc EGM96.gfc   # an ICGEM model instead of the synthetic one
python run_all.py --N 180 --skip efficiency --out my_results
```

Each experiment also runs on its own: `sph_efficiency.py`, `sph_precision.py` and `sph_accuracy.py`. Each accepts `--quick` and `--backend rust`. `sph_efficiency.py` and `sph_precision.py` also accept `--plot-only` to redraw from cached results. The module docstrings list all options.

On one core, a full run takes about 14 min with the Rust backend and about 8 min with NumPy, which times a subsample of the points.

## Layout

| File | Content |
|---|---|
| `sph_core.py` | The six NumPy programmes, the recursion tables, the gravity model (synthetic Kaula or ICGEM `.gfc`), local ↔ Cartesian rotations and the digit metric |
| `sph_rs/` | Rust crate: one scalar kernel per programme, plus CLENSHAW per parallel |
| `sph_rust.py` | Python wrapper of the Rust binary. It is a drop-in replacement for `sph_core.synthesize`. |
| `sph_efficiency.py` | Table 7, Fig. 4, cost per point vs N |
| `sph_precision.py` | Figs. 5–10, precision summary, conditioning analysis |
| `sph_accuracy.py` | Eqs. 31–33 (Laplace and ALF identities), Fig. 3 (instability of Pines' original recursion) |
| `sph_style.py` | Figure style: column widths, Okabe–Ito colours, vector PDF output |
| `run_all.py` | Runs the three experiments in order |
| `results/`, `results_rust/` | Outputs of the NumPy and Rust backends; `run_all.log` and `run_all_rust.log` are their logs |

Outputs in each results folder:

- `table7.md`, `table7.csv`, `efficiency.json`, `fig04_ground_track.pdf`, `fig_efficiency.pdf`
- `fig05`–`fig09_precision_<METHOD>.pdf`, where `fig07b` is GOTTLIEB
- `fig10_precision_PINES16_METRIS8.pdf`, `fig_precision_summary.pdf`, `fig_precision_conditioning.pdf`, `precision_summary.csv`, `precision_cache.npz`
- `fig03_msr_instability.pdf`, `fig_accuracy_laplace.pdf`, `fig_accuracy_identities.pdf`, `accuracy_summary.csv`

## Results

All results below use N = 360 and the synthetic Kaula-rule model, which carries the EGM96 values of C20, C22, S22 and C30. They were run on one core of an Apple M3 Pro.

### Speed

| Method | Rust ORBIT1 [ms/pt] | Rust [cycles/coef] | Rust / LEG_a | paper / LEG_a | paper [cycles/coef] | NumPy ORBIT1 [ms/pt] | Rust GRID [s] | paper GRID [s] |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| LEGENDRE_a | 0.174 | 10.8 | 1.00 | 1.00 | 642 | 1.55 | 11.1 | 2390 |
| LEGENDRE_b | 0.497 | 30.8 | 2.86 | 4.80 | 3080 | 3.73 | 32.2 | 9793 |
| PINES | 0.419 | 26.0 | 2.41 | 3.03 | 1947 | 2.01 | 26.4 | 5275 |
| GOTTLIEB | 0.400 | 24.8 | 2.30 | – | – | 1.84 | 25.8 | – |
| METRIS | 0.333 | 20.7 | 1.92 | 1.64 | 1055 | 0.75 | 21.3 | 3271 |
| CLENSHAW | 0.590 | 36.5 | 3.39 | 0.48 | 306 | 2.04 | 0.173 | 8 |

How to read the table:

- **ORBIT1:** 64,440 points along a near-polar LEO.
- **GRID:** 64,440 points on a 1° × 1° grid at 250 km.
- **cycles/coef:** CPU time per point × 4.05 GHz ÷ 65,341 coefficients.
- **Ratio and cycle columns:** ORBIT1.
- **NumPy ORBIT1:** extrapolated from 6,000 evenly spaced points (600 for LEGENDRE_b); the Rust kernels time every point.
- **Paper columns:** a Pentium 4 at 1.3 GHz. The paper's timings include binary I/O.
- **CLENSHAW on the GRID:** computed per parallel, which reuses the recursions along each latitude.

### Precision

Coinciding decimal digits against CLENSHAW*16 (80-bit), on 4,296 points of the paper's GRID. Each entry is median / minimum over the points, for the worst component of the group.

NumPy backend:

| Method | V | ∇V local (worst) | ∇∇V local (worst) | ∇V Cartesian (worst) | ∇∇V Cartesian (worst) |
|---|---:|---:|---:|---:|---:|
| LEGENDRE_a | 15.3 / 14.3 | 14.9 / 11.6 | 14.5 / 10.4 | 15.3 / 14.4 | 15.1 / 12.3 |
| LEGENDRE_b | 15.9 / 15.1 | 15.0 / 11.6 | 14.5 / 10.3 | – | – |
| PINES | 15.3 / 14.3 | 11.9 / 8.8 | 12.2 / 8.3 | 15.3 / 14.4 | 14.9 / 12.4 |
| GOTTLIEB | 15.3 / 14.3 | 11.8 / 8.4 | 12.2 / 8.4 | 15.3 / 14.4 | 15.1 / 12.4 |
| METRIS | 15.3 / 14.6 | 10.9 / 7.5 | 11.4 / 7.1 | 15.2 / 14.3 | 15.1 / 11.4 |
| CLENSHAW | 15.5 / 14.5 | 14.9 / 11.2 | 14.5 / 10.1 | – | – |

<details><summary>Rust backend</summary>

| Method | V | ∇V local (worst) | ∇∇V local (worst) | ∇V Cartesian (worst) | ∇∇V Cartesian (worst) |
|---|---:|---:|---:|---:|---:|
| LEGENDRE_a | 15.3 / 14.3 | 14.9 / 11.3 | 14.5 / 10.3 | 15.3 / 14.4 | 15.1 / 12.3 |
| LEGENDRE_b | 15.3 / 14.3 | 14.9 / 11.1 | 14.5 / 10.4 | – | – |
| PINES | 15.3 / 14.3 | 11.9 / 8.8 | 12.3 / 8.3 | 15.3 / 14.4 | 14.9 / 12.4 |
| GOTTLIEB | 15.3 / 14.3 | 11.8 / 8.4 | 12.2 / 8.1 | 15.2 / 14.4 | 15.1 / 12.4 |
| METRIS | 15.3 / 14.6 | 10.9 / 7.5 | 11.4 / 7.1 | 15.2 / 14.3 | 15.2 / 11.4 |
| CLENSHAW | 15.4 / 14.4 | 14.9 / 11.2 | 14.5 / 10.1 | – | – |

</details>

### Accuracy

These are exact identities, so they test the arithmetic.

- **Laplace (Eq. 31), all programmes.** The residual ∇²V (trace of the tensor) has a median of about 4e-22–3e-21 s⁻² and a maximum of about 2e-20 s⁻². Relative to |V_rr| that is at most 8e-15.
- **ALF identities (Eqs. 32–33, FCR in double).** The relative residual is about 1e-14 median and below 8e-13 maximum. The paper reports 1e-14 and 1e-13.
- **Pines' original modified sub-sectorial recursion (Fig. 3).** In double precision it diverges, with relative errors up to 10^120 at 85° latitude. The normalised column recursion used by PINES and GOTTLIEB does not.

## Findings

**1. The speed ranking of the algorithms is reproduced, but not the paper's absolute or CLENSHAW-on-orbit numbers.**

- The Rust kernels reproduce the ordering LEGENDRE_a < METRIS < PINES ≈ GOTTLIEB < LEGENDRE_b. LEGENDRE_b, PINES and METRIS cost 2.9, 2.4 and 1.9 times LEGENDRE_a; the paper has 4.8, 3.0 and 1.6.
- Absolute times are 160–310× shorter. That comes from the hardware, from the I/O included in the paper's timings, and from the Matlab layer in LEGENDRE_b.

**2. CLENSHAW at scattered orbit points is slower than LEGENDRE_a here, while the paper found it about 2× faster.**

- The paper's codes spent 300–3000 cycles per coefficient. The kernels here spend 10–40, so the 2009 ranking reflects those implementations (x87, small caches, stored derivative triangles) rather than the algorithms.
- Turning off auto-vectorisation only moves the CLENSHAW / LEGENDRE_a ratio from about 3.4–3.7 to 2.6.
- On a grid, CLENSHAW per parallel is by far the fastest method: about 64× faster than LEGENDRE_a here (11.1 s vs 0.17 s), and 300× in the paper.

**3. NumPy timings rank NumPy implementations.**

- METRIS maps onto BLAS matrix products and looks fastest.
- For the other methods, memory traffic dominates.
- The NumPy runs here were also x86-64 under Rosetta. Compare algorithms with the Rust backend.

**4. Cartesian methods are not less accurate; the rotation to the local frame is ill-conditioned.**

- PINES, GOTTLIEB and METRIS give about 15 digits (median) in their native Cartesian components. Their worst points (11.4–12.7 digits) are where a tensor component such as V_xy passes near zero, and LEGENDRE_a rotated to Cartesian has the same worst points.
- They lose digits only in the small local components (V_λ, V_λφ, V_λr) after rotation. The loss follows −log₁₀(u·κ), where κ is the condition number of the component with respect to the Cartesian vector or tensor (`fig_precision_conditioning.pdf`).
- Rounding the *exact* Cartesian result to double and rotating it reproduces the same loss.
- This explains the precision gap that the paper leaves open in Sect. 7.2.

**5. GOTTLIEB is as precise and as fast as PINES.**

- Its precision matches PINES: against the 80-bit reference, the medians of every component agree within 0.2 digits, and the minima within 0.4.
- Compiled, it costs the same per coefficient (~25 cycles).
- In NumPy it is about 10 % faster than PINES. METRIS is faster than both there, because its sums map onto BLAS.
- Its practical advantage is one ALF row instead of three and simpler indexing, not arithmetic.

**Which method to use:**

- **Global grids:** CLENSHAW per parallel.
- **Scattered points with spherical output, away from the poles:** LEGENDRE_a (fastest compiled). It needs pole handling and a rotation to produce Cartesian accelerations.
- **Orbit propagation and orbit determination in body-fixed Cartesian coordinates:** METRIS / Cunningham, the fastest singularity-free method (~20 cycles per coefficient), then PINES ≈ GOTTLIEB (~25). All three deliver the Cartesian acceleration to 14+ digits and the gradient tensor to about 15 digits (median).
- **Very high degree (e.g. EGM2008 at N = 2190):** double-precision recursions underflow near the poles. Use X-numbers or Holmes–Featherstone scaling. This range is **not tested here.**

## What operational software uses

This is from public documentation and source code, as of 2026-11.

| Tool | Formulation | Family here |
|---|---|---|
| Scarabaeus (CCAR-ORCCA) | Normalised Cunningham solid-harmonic recursion B_nm = V_nm + iW_nm (Montenbruck & Gill Sect. 3.2), analytic Jacobian and C/S partials. Its docstrings say "Pines/Eckman". | METRIS |
| Basilisk | Normalised Pines | PINES |
| GMAT | Normalised Pines, with optional gradient for the state transition matrix | PINES |
| Nyx | Pines, ported from GMAT | PINES |
| Tudat | Fully normalised ALFs in spherical coordinates, then rotation | LEGENDRE_a |
| Orekit | Holmes–Featherstone (modified forward row method) | – |
| MONTE (JPL) | Not public | – |

None of these tools uses Gottlieb's formulation. Gottlieb's formulation comes from NASA Johnson Space Center; Eckman et al. (2016) compare it with Pines and Lear.

## Caveats

- The `*16` reference is 80-bit extended precision (about 19 digits), not the paper's 128-bit REAL*16. That is enough to resolve double-precision errors.
- The default gravity model is synthetic. Use `--gfc` with an ICGEM file (EGM96, EGM2008) for a real field.
- Cycles per coefficient assume the 4.05 GHz maximum clock of the M3 Pro performance cores.

## References

- Fantino, E., Casotto, S. (2009). Methods of harmonic synthesis for global geopotential models and their first-, second- and third-order gradients. *J. Geod.* 83, 595–619.
- Pines, S. (1973). Uniform representation of the gravitational potential and its derivatives. *AIAA J.* 11(11), 1508–1511.
- Gottlieb, R. G. (1993). *Fast gravity, gravity partials, normalized gravity, gravity gradient torque and magnetic field: derivation, code and data.* NASA CR-188243.
- Metris, G., Xu, J., Wytrzyszczak, I. (1999). Derivatives of the gravity potential with respect to rectangular coordinates. *Celest. Mech. Dyn. Astron.* 71, 137–151.
- Cunningham, L. E. (1970). On the computation of the spherical harmonic terms needed during the numerical integration of the orbital motion of an artificial satellite. *Celest. Mech.* 2, 207–216.
- Eckman, R. A., Brown, A. J., Adamo, D. R. (2016). *Normalization and implementation of three gravitational acceleration models.* NASA/TP-2016-218604.
- Holmes, S. A., Featherstone, W. E. (2002). A unified approach to the Clenshaw summation and the recursive computation of very high degree and order normalised associated Legendre functions. *J. Geod.* 76, 279–299.
- Montenbruck, O., Gill, E. (2000). *Satellite Orbits.* Springer, Sect. 3.2.
