"""
sph_precision.py -- relative numerical precision of the five synthesis programmes
(Figs. 5-10 of Fantino & Casotto 2009).

Precision = number of coinciding decimal digits, -log10(|x*8 - x_ref| / |x_ref|), of each of
the ten quantities (V, its gradient and its gradient tensor in the local frame) computed in
double precision (*8), against a reference computed in extended precision (*16):

    Figs. 5-9   CLENSHAW*16 vs LEGENDRE_a*8, LEGENDRE_b*8, PINES*8, METRIS*8, CLENSHAW*8
    Fig. 10     PINES*16 vs METRIS*8 (Cartesian components)

"*16" is numpy.longdouble: 80-bit x87 extended precision on x86_64 (64-bit mantissa, 19-20
significant digits), not the 128-bit REAL*16 of the paper; it is ~3 digits above double, so the
double-precision digits (up to ~16) are resolved. On platforms where longdouble == double
(Windows, Apple-silicon native builds) the reference is no better than the methods: a warning
is printed.

Test points: the GRID of the paper (h = 250 km, phi = -89..89 deg) on every --lon-step degrees
of longitude (default 15 deg: 179 x 24 = 4,296 points; --lon-step 1 is the full 64,440-point
GRID). Fig. 10 uses a reduced set (--cart-lat-step, --cart-lon-step) because PINES*16 is slow.

Beyond the paper, fig_precision_summary also compares the Cartesian components with
CLENSHAW*16 rotated to the Cartesian frame in extended precision, which separates the error
of the synthesis from that of the local <-> Cartesian rotation (the open question of Sect. 7.2).

Outputs (in --out, default ./results):
    precision_cache.npz                      outputs and digits (re-used by sph_accuracy.py)
    precision_summary.csv                    median / 1st percentile / minimum digits
    fig05_precision_LEGENDRE_a.pdf ... fig09_precision_CLENSHAW.pdf, fig10_precision_PINES16_METRIS8.pdf
    fig_precision_summary.pdf
    fig_precision_conditioning.pdf           digits vs condition number of the local components

    python sph_precision.py                # ~2-5 min (mostly PINES*16 for Fig. 10)
    python sph_precision.py --quick        # smoke test (~1 min)
    python sph_precision.py --lon-step 1   # full GRID (~10 min)
    python sph_precision.py --plot-only    # redraw from precision_cache.npz
    python sph_precision.py --backend rust # *8 results from the compiled kernels (sph_rs/)

With --backend rust every *8 result comes from the Rust binary (IEEE double, no fast-math, no
FMA contraction); the *16 references are always NumPy longdouble.
"""

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")  # before numpy is imported

import argparse
import time
from pathlib import Path

import numpy as np

import sph_core as sc
import sph_rust as sr
import sph_style as ss
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE = Path(__file__).resolve().parent
LD = np.longdouble
EXTENDED = np.finfo(LD).eps < 1e-18

# panel order of Figs. 5-10 (row-wise, 5 x 2)
PANEL_LOCAL = ("V", "Vl", "Vp", "Vr", "Vll", "Vpp", "Vrr", "Vlp", "Vlr", "Vpr")
PANEL_CART = ("V", "Vx", "Vy", "Vz", "Vxx", "Vyy", "Vzz", "Vxz", "Vxy", "Vyz")
FIG_NUM = {"LEGENDRE_a": 5, "LEGENDRE_b": 6, "PINES": 7, "METRIS": 8, "CLENSHAW": 9}
# Cartesian comparisons of the summary figure: (label, test, reference)
CART_PAIRS = (
    ("METRIS*8 vs PINES*16 (Fig. 10)", "METRIS", "PINES16"),
    ("METRIS*8 vs CLENSHAW*16", "METRIS", "CLEN16"),
    ("PINES*8 vs CLENSHAW*16", "PINES", "CLEN16"),
    ("LEGENDRE_a*8 rotated vs CLENSHAW*16", "LEGENDRE_a", "CLEN16"),
)
PAIR_COLOR = {"METRIS": ss.METHOD_COLOR["METRIS"], "PINES": ss.METHOD_COLOR["PINES"],
              "LEGENDRE_a": ss.METHOD_COLOR["LEGENDRE_a"]}


# ============================================================
# COMPUTATION
# ============================================================
def _grid(lat_step, lon_step, h):
    lon = np.arange(-179.0, 180.0 + 1e-9, lon_step)
    lat = np.arange(-89.0, 89.0 + 1e-9, lat_step)
    lam, phi, r = sc.grid_points(h, lon, lat)
    return lon, lat, lam, phi, r


def compute_local(model, lat_step, lon_step, h=250e3, backend="numpy"):
    """Local quantities of every method (*8) and CLENSHAW*16 on the precision grid."""
    lon, lat, lam, phi, r = _grid(lat_step, lon_step, h)
    out = {"lon": lon, "lat": lat, "lam": lam, "phi": phi, "r": r}
    B = sr if backend == "rust" else sc
    t0 = time.perf_counter()
    out["ref"] = sc.clenshaw_grid(model, np.radians(lon), np.radians(lat), r[0], dtype=LD)
    print(f"[precision] CLENSHAW*16 reference, {lam.size} points: {time.perf_counter() - t0:.1f} s")
    for meth in sc.METHODS:
        t0 = time.perf_counter()
        if meth == "CLENSHAW":  # the paper's grid arrangement (per parallel)
            res = B.clenshaw_grid(model, np.radians(lon), np.radians(lat), r[0])
        else:
            res = B.synthesize(meth, model, lam, phi, r)
        out[meth] = np.asarray(res, dtype=np.float64)
        print(f"[precision] {meth:10s}*8: {time.perf_counter() - t0:6.1f} s")
    return out


def compute_cart(model, lat_step, lon_step, h=250e3, backend="numpy"):
    """Cartesian components on the reduced set: PINES*16, CLENSHAW*16 rotated in *16,
    METRIS*8, PINES*8, LEGENDRE_a*8 rotated in *8."""
    lon, lat, lam, phi, r = _grid(lat_step, lon_step, h)
    out = {"lon": lon, "lat": lat}
    t0 = time.perf_counter()
    _, out["PINES16"] = sc.pines(model, lam, phi, r, dtype=LD, cartesian=True)
    print(f"[precision] PINES*16 (Cartesian), {lam.size} points: {time.perf_counter() - t0:.1f} s")
    loc16 = sc.clenshaw_grid(model, np.radians(lon).astype(LD), np.radians(lat).astype(LD),
                             r[0], dtype=LD)
    out["CLEN16"] = sc.local_to_cart(loc16, lam.astype(LD), phi.astype(LD))
    B = sr if backend == "rust" else sc
    for meth in ("METRIS", "PINES", "LEGENDRE_a"):
        _, out[meth] = B.synthesize(meth, model, lam, phi, r, cartesian=True)
    return out


def load_or_compute(model, a, outdir):
    """Cached outputs (precision_cache.npz) unless --recompute or the setup changed."""
    path = outdir / "precision_cache.npz"
    key = np.array([model.N, a.lat_step, a.lon_step, a.cart_lat_step, a.cart_lon_step, a.seed]
                   + ([1.0] if a.backend == "rust" else []))
    if path.exists() and not a.recompute:
        z = np.load(path, allow_pickle=False)
        if "key" in z and z["key"].shape == key.shape and np.allclose(z["key"], key):
            loc = {k[4:]: z[k] for k in z.files if k.startswith("loc_")}
            car = {k[4:]: z[k] for k in z.files if k.startswith("car_")}
            print(f"[precision] loaded {path.name}")
            return loc, car
        print("[precision] cache does not match the setup: recomputing")
    loc = compute_local(model, a.lat_step, a.lon_step, backend=a.backend)
    car = compute_cart(model, a.cart_lat_step, a.cart_lon_step, backend=a.backend)
    payload = {"key": key}
    payload.update({f"loc_{k}": v for k, v in loc.items()})
    payload.update({f"car_{k}": v for k, v in car.items()})
    import io

    buf = io.BytesIO()
    np.savez_compressed(buf, **payload)
    ss.write_output(path, buf.getvalue(), tag="Data")
    return loc, car


def local_digits(loc):
    return {m: sc.digits(loc[m], loc["ref"]) for m in sc.METHODS}


def cart_digits(car):
    return {f"{t}|{r}": sc.digits(car[t], car[r]) for _, t, r in CART_PAIRS}


def per_lat(lat_pts, d, lat):
    """Median and minimum of d (K,) per parallel."""
    step = lat[1] - lat[0] if lat.size > 1 else 1.0
    idx = np.rint((lat_pts - lat[0]) / step).astype(int)
    med = np.array([np.median(d[idx == i]) for i in range(lat.size)])
    mn = np.array([np.min(d[idx == i]) for i in range(lat.size)])
    return med, mn


# ============================================================
# FIGURES
# ============================================================
def _ylim(all_digits):
    lo = min(float(np.min(d)) for d in all_digits)
    return (min(6.0, np.floor(lo) - 0.3), 16.6)


def fig_precision(lat_pts, lat, dig, quants, color, title, name, ylim):
    """5 x 2 panels of digits vs latitude (Figs. 5-10): all points, median and minimum per
    parallel."""
    fig, axes = plt.subplots(5, 2, figsize=(ss.TXT_W, 6.4), sharex=True, sharey=True,
                             layout="constrained")
    lat_deg = np.asarray(lat_pts)  # [deg]
    names = sc.QUANTS if quants is PANEL_LOCAL else sc.CART
    dark = tuple(0.6 * np.asarray(plt.matplotlib.colors.to_rgb(color)))
    for p, (ax, q) in enumerate(zip(axes.ravel(), quants)):
        d = dig[names.index(q)]
        ax.scatter(lat_deg, d, s=1.2, lw=0, color=ss._tint(color, 0.45), rasterized=True,
                   zorder=1)
        med, mn = per_lat(lat_deg, d, lat)
        ax.plot(lat, mn, color="0.25", lw=0.6, ls=(0, (2, 1)), zorder=2)
        ax.plot(lat, med, color=dark, lw=1.0, zorder=3)
        ss.panel_title(ax, "abcdefghij"[p], sc.QUANT_TEX[q])
        ax.text(0.99, 0.04, f"median {np.median(d):.1f}, min {np.min(d):.1f}",
                transform=ax.transAxes, ha="right", va="bottom", fontsize=6, color="0.2")
        ax.set_ylim(*ylim)
        ax.set_xlim(-90, 90)
        ax.set_xticks([-90, -60, -30, 0, 30, 60, 90])
        ax.set_yticks(np.arange(np.ceil(ylim[0] / 2) * 2, 16.1, 2))
        if p % 2 == 0:
            ax.set_ylabel("Digits [-]")
        if p >= 8:
            ax.set_xlabel(r"Latitude $\varphi$ [deg]")
    handles = [
        Line2D([], [], ls="none", marker="o", ms=3, mec="none", color=ss._tint(color, 0.45),
               label="grid points"),
        Line2D([], [], color=dark, lw=1.0, label="median per parallel"),
        Line2D([], [], color="0.25", lw=0.6, ls=(0, (2, 1)), label="minimum per parallel"),
    ]
    fig.legend(handles=handles, loc="outside upper center", ncol=3, title=title,
               title_fontsize=7.5)
    ss.save_figure(fig, name)


def _stats(d):
    return np.median(d, axis=-1), np.percentile(d, 1, axis=-1), np.min(d, axis=-1)


def fig_summary(dl, dc, n_loc, n_car, name="fig_precision_summary"):
    """(a) median and 1st-percentile digits per local quantity and method (Figs. 5-9 in one
    panel); (b) the same for the Cartesian comparisons (Fig. 10 + rotation test)."""
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(ss.TXT_W, 2.7), layout="constrained",
                                 gridspec_kw={"width_ratios": [1.0, 1.0]})
    xq = np.arange(len(PANEL_LOCAL))
    off = np.linspace(-0.3, 0.3, len(sc.METHODS))
    for j, meth in enumerate(sc.METHODS):
        med, p1, _ = _stats(dl[meth])
        idx = [sc.QUANTS.index(q) for q in PANEL_LOCAL]
        x = xq + off[j]
        a1.vlines(x, p1[idx], med[idx], color=ss.METHOD_COLOR[meth], lw=0.8)
        a1.plot(x, med[idx], ls="none", marker=ss.METHOD_MARKER[meth], ms=3.2,
                color=ss.METHOD_COLOR[meth], label=ss.method_label(meth) + "*8")
        a1.plot(x, p1[idx], ls="none", marker="_", ms=4, color=ss.METHOD_COLOR[meth])
    a1.set_xticks(xq)
    a1.set_xticklabels([sc.QUANT_TEX[q] for q in PANEL_LOCAL])
    a1.set_ylabel("Digits [-]")
    ss.panel_title(a1, "a", f"Local frame vs CLENSHAW*16 ({n_loc:,} points)")
    a1.grid(axis="x", visible=False)

    xq = np.arange(len(PANEL_CART))
    off = np.linspace(-0.27, 0.27, len(CART_PAIRS))
    mk = ("D", "o", "^", "s")
    for j, (lab, t, r) in enumerate(CART_PAIRS):
        med, p1, _ = _stats(dc[f"{t}|{r}"])
        idx = [sc.CART.index(q) for q in PANEL_CART]
        x = xq + off[j]
        col = PAIR_COLOR[t]
        filled = r == "PINES16"
        a2.vlines(x, p1[idx], med[idx], color=col, lw=0.8)
        a2.plot(x, med[idx], ls="none", marker=mk[j], ms=3.2, color=col,
                mfc=col if filled else "white", mew=0.8, label=ss.tex_escape(lab))
        a2.plot(x, p1[idx], ls="none", marker="_", ms=4, color=col)
    a2.set_xticks(xq)
    a2.set_xticklabels([sc.QUANT_TEX[q] for q in PANEL_CART])
    ss.panel_title(a2, "b", f"Cartesian components ({n_car:,} points)")
    a2.grid(axis="x", visible=False)
    lo = min(a1.get_ylim()[0], a2.get_ylim()[0])
    for ax in (a1, a2):
        ax.set_ylim(max(lo, 0), 16.8)
    fig.legend(*a1.get_legend_handles_labels(), loc="outside upper center", ncol=5)
    a2.legend(loc="lower left", ncol=1, fontsize=6, handletextpad=0.3)
    fig.text(0.5, -0.02, "marker: median over the points; bar: down to the 1st percentile",
             ha="center", va="top", fontsize=6.5, color="0.3")
    ss.save_figure(fig, name)


def condition_numbers(ref):
    """
    Relative condition number of each local component with respect to perturbations of the
    Cartesian vector/tensor of the same size: kappa_i = |grad V| / |V_i| for the gradient,
    ||grad grad V||_F / |V_ij| for the tensor, 1 for V. An algorithm that delivers the Cartesian
    components with a norm-wise relative error u yields V_i with relative error ~ u kappa_i.
    """
    ref = np.asarray(ref, dtype=np.float64)
    g = np.sqrt(ref[1] ** 2 + ref[2] ** 2 + ref[3] ** 2)
    T = np.sqrt(ref[4] ** 2 + ref[6] ** 2 + ref[9] ** 2
                + 2 * (ref[5] ** 2 + ref[7] ** 2 + ref[8] ** 2))
    kap = np.ones_like(ref)
    with np.errstate(divide="ignore"):
        for i in (1, 2, 3):
            kap[i] = g / np.abs(ref[i])
        for i in range(4, 10):
            kap[i] = T / np.abs(ref[i])
    return np.minimum(kap, 1e300)


def rounded_cartesian_digits(loc):
    """Digits of the local components obtained by rounding the *exact* Cartesian components
    (CLENSHAW*16 rotated in *16) to double and rotating them back in double: the loss that
    any Cartesian-frame algorithm suffers, however accurate."""
    lam16, phi16 = loc["lam"].astype(LD), loc["phi"].astype(LD)
    cart = sc.local_to_cart(loc["ref"], lam16, phi16).astype(np.float64)
    back = sc.cart_to_local(cart, loc["lam"], loc["phi"])
    return sc.digits(back, loc["ref"])


def fig_conditioning(loc, dl, name="fig_precision_conditioning"):
    """Digits of the local components vs their condition number log10(kappa): the Cartesian
    codes follow the bound -log10(u kappa); the spherical codes do not."""
    kap = np.log10(condition_numbers(loc["ref"]))
    u = np.finfo(np.float64).eps / 2
    cases = [("LEGENDRE_a*8", dl["LEGENDRE_a"], ss.METHOD_COLOR["LEGENDRE_a"]),
             ("PINES*8", dl["PINES"], ss.METHOD_COLOR["PINES"]),
             ("METRIS*8", dl["METRIS"], ss.METHOD_COLOR["METRIS"]),
             ("exact Cartesian rounded to *8, rotated", rounded_cartesian_digits(loc), "0.35")]
    fig, axes = plt.subplots(1, 4, figsize=(ss.TXT_W, 2.35), sharex=True, sharey=True,
                             layout="constrained")
    xx = np.array([0.0, 8.0])
    groups = ((slice(0, 1), "$V$"), (slice(1, 4), r"$\nabla V$"),
              (slice(4, 10), r"$\nabla\nabla V$"))
    for p, (ax, (lab, d, col)) in enumerate(zip(axes, cases)):
        for (sl, glab), qc in zip(groups, ss.QUANT_COLOR):
            ax.scatter(kap[sl].ravel(), d[sl].ravel(), s=0.8, lw=0, color=ss._tint(qc, 0.5),
                       rasterized=True, label=glab)
        ax.plot(xx, -np.log10(u) - xx, color="k", lw=0.9, label=r"$-\log_{10}(u\,\kappa)$")
        ax.fill_between(xx, -np.log10(u) - xx - 1, -np.log10(u) - xx + 1, color=ss.C_SPAN,
                        lw=0, zorder=0)
        ss.panel_title(ax, "abcd"[p], ss.tex_escape(lab) if p < 3 else "Rounding only")
        resid = np.median(d[1:] - (-np.log10(u) - kap[1:]))
        ax.text(0.04, 0.04, f"median excess over\nthe bound: {resid:+.1f} digits",
                transform=ax.transAxes, ha="left", va="bottom", fontsize=6,
                color=col if p < 3 else "0.2")
        ax.set_xlabel(r"$\log_{10}\kappa$ [-]")
    axes[0].set_ylabel("Digits [-]")
    axes[0].set_xlim(0, 7.5)
    axes[0].set_ylim(5, 17)
    h = [Line2D([], [], ls="none", marker="o", ms=3.5, color=ss._tint(qc, 0.6))
         for qc in ss.QUANT_COLOR]
    h.append(Line2D([], [], color="k", lw=0.9))
    l = [g[1] for g in groups] + [r"$-\log_{10}(u\,\kappa)$, $u = 2^{-53}$ ($\pm 1$ shaded)"]
    fig.legend(h, l, loc="outside upper center", ncol=4)
    ss.save_figure(fig, name)


def write_summary(dl, dc, outdir):
    rows = ["comparison,quantity,median,p01,min"]
    for meth in sc.METHODS:
        med, p1, mn = _stats(dl[meth])
        for i, q in enumerate(sc.QUANTS):
            rows.append(f"CLENSHAW*16 vs {meth}*8,{q},{med[i]:.2f},{p1[i]:.2f},{mn[i]:.2f}")
    for lab, t, r in CART_PAIRS:
        med, p1, mn = _stats(dc[f"{t}|{r}"])
        for i, q in enumerate(sc.CART):
            rows.append(f"{lab},{q},{med[i]:.2f},{p1[i]:.2f},{mn[i]:.2f}")
    ss.write_output(outdir / "precision_summary.csv", "\n".join(rows) + "\n", tag="Data")
    print("[precision] median digits (local frame, vs CLENSHAW*16):")
    print("            " + " ".join(f"{q:>5s}" for q in sc.QUANTS))
    for meth in sc.METHODS:
        print(f"{meth:11s} " + " ".join(f"{v:5.1f}" for v in np.median(dl[meth], axis=1)))
    print("[precision] median digits (Cartesian):")
    print(" " * 38 + " ".join(f"{q:>5s}" for q in sc.CART))
    for lab, t, r in CART_PAIRS:
        print(f"{lab:37s} " + " ".join(f"{v:5.1f}" for v in np.median(dc[f'{t}|{r}'], axis=1)))


# ============================================================
# MAIN
# ============================================================
def parse(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--N", type=int, default=360, help="maximum degree (paper: 360)")
    ap.add_argument("--gfc", type=str, default=None, help="ICGEM .gfc model (default: synthetic)")
    ap.add_argument("--seed", type=int, default=2009, help="seed of the synthetic model")
    ap.add_argument("--lat-step", type=float, default=1.0, help="latitude step [deg]")
    ap.add_argument("--lon-step", type=float, default=15.0, help="longitude step [deg]")
    ap.add_argument("--cart-lat-step", type=float, default=2.0, help="Fig. 10 latitude step")
    ap.add_argument("--cart-lon-step", type=float, default=60.0, help="Fig. 10 longitude step")
    ap.add_argument("--quick", action="store_true", help="small smoke-test run")
    ap.add_argument("--recompute", action="store_true", help="ignore precision_cache.npz")
    ap.add_argument("--plot-only", action="store_true", help="(alias) use the cache")
    ap.add_argument("--backend", choices=("numpy", "rust"), default="numpy",
                    help="origin of the *8 results: numpy or rust (compiled kernels, sph_rs/)")
    ap.add_argument("--out", type=str, default=None,
                    help="output folder (default ./results, or ./results_rust)")
    a = ap.parse_args(argv)
    a.out = a.out or str(HERE / ("results_rust" if a.backend == "rust" else "results"))
    if a.quick:
        a.lat_step, a.lon_step, a.cart_lat_step, a.cart_lon_step = 2.0, 45.0, 4.0, 120.0
    return a


def get_model(a):
    if a.gfc:
        return sc.GravityModel.from_gfc(a.gfc, a.N)
    return sc.GravityModel.synthetic(a.N, seed=a.seed)


def main(argv=None):
    a = parse(argv)
    if not EXTENDED:
        print("[precision] WARNING: numpy.longdouble is not extended precision on this platform;"
              " the *16 references are only double precision")
    outdir = ss.set_figure_dir(a.out)
    ss.apply_style()
    model = get_model(a)
    loc, car = load_or_compute(model, a, outdir)
    dl, dc = local_digits(loc), cart_digits(car)
    write_summary(dl, dc, outdir)
    ylim = _ylim(list(dl.values()) + [dc["METRIS|PINES16"]])
    for meth in sc.METHODS:
        title = ss.tex_escape(f"Relative numerical precision: CLENSHAW*16 vs {meth}*8")
        fig_precision(np.degrees(loc["phi"]), loc["lat"], dl[meth], PANEL_LOCAL, ss.METHOD_COLOR[meth],
                      title, f"fig{FIG_NUM[meth]:02d}_precision_{meth}", ylim)
    lat_c = np.arange(-89.0, 89.0 + 1e-9, a.cart_lat_step)
    LATc = np.repeat(lat_c, car["lon"].size)
    fig_precision(LATc, lat_c, dc["METRIS|PINES16"], PANEL_CART, ss.METHOD_COLOR["METRIS"],
                  "Relative numerical precision: PINES*16 vs METRIS*8 (Cartesian)",
                  "fig10_precision_PINES16_METRIS8", ylim)
    fig_summary(dl, dc, loc["lam"].size, LATc.size)
    fig_conditioning(loc, dl)
    return loc, car, dl, dc


if __name__ == "__main__":
    main()
