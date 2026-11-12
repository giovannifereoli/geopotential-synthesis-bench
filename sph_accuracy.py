"""
sph_accuracy.py -- numerical accuracy of the five synthesis programmes (Sect. 7.3 of Fantino &
Casotto 2009) and the instability of Pines' original recursion (Fig. 3).

Exact identities checked:

    Eq. 31  Laplace          V_ll + V_pp + V_rr = 0             (all five programmes)
    Eq. 32  sum_nm Pbar_nm^2           = (N + 1)^2              (fnALFs of LEGENDRE_a)
    Eq. 33  sum_nm (dPbar_nm/dphi)^2   = N (N + 1)^2 (N + 2) / 4

Eq. 31 uses the GRID outputs of sph_precision.py (precision_cache.npz, computed if missing);
the trace is summed in extended precision, so that it measures the synthesis error only.
Eqs. 32-33 are evaluated on every parallel (-89..89 deg) with
    FCR*8                   standard forward column recursion (LEGENDRE_a, CLENSHAW)
    Helmholtz FCR*8         Abar_nm(sin phi) cos^m(phi) (PINES with the FCR)
    FCR*16                  the same in extended precision (floor of the test)
and the modified sub-sectorial recursion (MSR) of Pines (1973) is reported separately: its error
is many orders of magnitude off the scale of the other variants (Fig. 3).

Paper (Sect. 7.3): Eq. 31 holds to 1e-18 - 1e-19 for all programmes, Eqs. 32-33 to 1e-14
and 1e-13 for LEGENDRE_a.

Outputs (in --out, default ./results):
    accuracy_summary.csv
    fig03_msr_instability.pdf      log10 relative error of MSR*8 vs FCR*16 (n, m) + per degree
    fig_accuracy_laplace.pdf       Eq. 31 per parallel, absolute [s^-2] and relative to |V_rr|
    fig_accuracy_identities.pdf    Eqs. 32-33 per parallel

    python sph_accuracy.py               # ~1 min once precision_cache.npz exists
    python sph_accuracy.py --quick       # smoke test (uses the --quick precision cache)
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
import sph_precision as sp
import sph_style as ss
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE = Path(__file__).resolve().parent
LD = np.longdouble
EPS = np.finfo(np.float64).eps
PAPER_LAPLACE = (1e-19, 1e-18)  # [s^-2], Sect. 7.3
PAPER_EQ32, PAPER_EQ33 = 1e-14, 1e-13
I_LL, I_PP, I_RR = (sc.QUANTS.index(q) for q in ("Vll", "Vpp", "Vrr"))


# ============================================================
# EQ. 31: LAPLACE
# ============================================================
def laplace(loc):
    """|V_ll + V_pp + V_rr| per point [s^-2] and relative to |V_rr|, for every method and the
    CLENSHAW*16 reference ("REF16"). Sums in extended precision."""
    out = {}
    vrr = np.abs(loc["ref"][I_RR].astype(LD))
    for m in sc.METHODS + ("ref",):
        x = loc[m].astype(LD)
        tr = np.abs(x[I_LL] + x[I_PP] + x[I_RR])
        key = "REF16" if m == "ref" else m
        out[key] = (tr.astype(np.float64), (tr / vrr).astype(np.float64))
    return out


def per_lat_max(lat_pts, d, lat):
    step = lat[1] - lat[0] if lat.size > 1 else 1.0
    idx = np.rint((lat_pts - lat[0]) / step).astype(int)
    return np.array([np.max(d[idx == i]) for i in range(lat.size)])


def fig_laplace(loc, lap, name="fig_accuracy_laplace"):
    lat_pts = np.degrees(loc["phi"].astype(np.float64))
    lat = loc["lat"]
    fig, axes = plt.subplots(1, 2, figsize=(ss.TXT_W, 2.5), sharex=True, layout="constrained")
    for k, ax in enumerate(axes):
        for m in sc.METHODS:
            y = per_lat_max(lat_pts, np.maximum(lap[m][k], 1e-30), lat)
            ax.plot(lat, y, color=ss.METHOD_COLOR[m], lw=0.9, label=ss.method_label(m) + "*8")
        y = per_lat_max(lat_pts, np.maximum(lap["REF16"][k], 1e-30), lat)
        ax.plot(lat, y, color=ss.C_REF, lw=0.9, ls=(0, (4, 1.5)), label="CLENSHAW*16")
        ax.set_yscale("log")
        ax.set_xlabel(r"Latitude $\varphi$ [deg]")
        ax.set_xlim(-90, 90)
        ax.set_xticks(np.arange(-90, 91, 30))
    a1, a2 = axes
    a1.axhspan(*PAPER_LAPLACE, color=ss.C_SPAN, lw=0, zorder=0)
    a1.text(-85, np.sqrt(np.prod(PAPER_LAPLACE)), "paper, Sect. 7.3", fontsize=6,
            color=ss.C_PAPER, va="center")
    a1.set_ylabel(r"$\max_\lambda |V_{\lambda\lambda} + V_{\varphi\varphi} + V_{rr}|$ [s$^{-2}$]")
    ss.panel_title(a1, "a", r"Laplace identity (Eq. 31), absolute")
    a2.axhline(EPS, color=ss.C_PAPER, lw=0.8, ls="--")
    a2.text(-85, EPS / 1.5, r"double $\epsilon = 2^{-52}$", fontsize=6, color=ss.C_PAPER,
            va="top")
    a2.set_ylabel(r"$\max_\lambda |\mathrm{trace}| \,/\, |V_{rr}|$ [-]")
    ss.panel_title(a2, "b", r"Relative to $|V_{rr}|$")
    lo = min(ax.get_ylim()[0] for ax in axes)
    a1.set_ylim(max(lo, 1e-26), 1e-17)
    a2.set_ylim(1e-20, 1e-13)
    fig.legend(*a1.get_legend_handles_labels(), loc="outside upper center", ncol=6)
    ss.save_figure(fig, name)


# ============================================================
# EQS. 32-33: SUMS OF fnALFs AND OF THEIR DERIVATIVES
# ============================================================
def _sums(P, phi, dtype):
    """Relative residuals of Eqs. 32 and 33 from a (K, N+1, N+1) triangle, accumulated in LD."""
    N = P.shape[1] - 1
    s32 = (P.astype(LD) ** 2).sum(axis=(1, 2))
    D = sc.dlegendre_dphi(P, np.asarray(phi, dtype), dtype)
    s33 = (D.astype(LD) ** 2).sum(axis=(1, 2))
    r32 = s32 / LD((N + 1) ** 2) - 1
    r33 = s33 / (LD(N) * LD(N + 1) ** 2 * LD(N + 2) / 4) - 1
    return np.abs(r32.astype(np.float64)), np.abs(r33.astype(np.float64))


def _helmholtz_P(phi, N, dtype, msr=False):
    phi = np.asarray(phi, dtype)
    u, c = np.sin(phi), np.cos(phi)
    A = (sc.helmholtz_msr if msr else sc.helmholtz_fcr)(u, N, dtype)
    with np.errstate(under="ignore", over="ignore", invalid="ignore"):
        cm = c[:, None] ** np.arange(N + 1, dtype=dtype)[None, :]
        return A * cm[:, None, :]


VARIANTS = {
    "FCR*8": dict(dtype=np.float64, helm=False, color=ss.C_FCR, ls="-"),
    "Helmholtz FCR*8": dict(dtype=np.float64, helm=True, color=ss.C_HFCR, ls=(0, (2.5, 1.2))),
    "FCR*16": dict(dtype=LD, helm=False, color=ss.C_REF, ls=(0, (4, 1.5))),
}


def identities(N, lat_deg, chunk=20):
    """{variant: (res32, res33)} on the parallels lat_deg, plus MSR*8 on the same parallels."""
    out = {}
    phi = np.radians(np.asarray(lat_deg, dtype=np.float64))
    for name, v in list(VARIANTS.items()) + [("MSR*8", dict(dtype=np.float64, helm="msr"))]:
        t0 = time.perf_counter()
        r32, r33 = [], []
        for s in range(0, phi.size, chunk):
            ph = phi[s:s + chunk].astype(v["dtype"])
            if v["helm"]:
                P = _helmholtz_P(ph, N, v["dtype"], msr=v["helm"] == "msr")
            else:
                P = sc.legendre_fcr(ph, N, v["dtype"])
            with np.errstate(all="ignore"):
                a, b = _sums(P, ph, v["dtype"])
            r32.append(a)
            r33.append(b)
        out[name] = (np.concatenate(r32), np.concatenate(r33))
        print(f"[accuracy] Eqs. 32-33 {name:16s} {time.perf_counter() - t0:6.1f} s")
    return out


def fig_identities(lat_deg, ids, name="fig_accuracy_identities"):
    fig, axes = plt.subplots(1, 2, figsize=(ss.TXT_W, 2.5), sharex=True, sharey=True,
                             layout="constrained")
    titles = (r"$\sum_{n,m} \bar P_{nm}^2 = (N+1)^2$ (Eq. 32)",
              r"$\sum_{n,m} (\partial \bar P_{nm}/\partial\varphi)^2 = N(N+1)^2(N+2)/4$ (Eq. 33)")
    paper = (PAPER_EQ32, PAPER_EQ33)
    for k, ax in enumerate(axes):
        for var, v in VARIANTS.items():
            y = np.maximum(ids[var][k], 1e-21)
            ax.plot(lat_deg, y, color=v["color"], ls=v["ls"], lw=0.9, label=ss.tex_escape(var))
        ax.axhline(paper[k], color=ss.C_PAPER, lw=1.0, ls=(0, (1, 1.2)),
                   label=ss.tex_escape("LEGENDRE_a*8, paper"))
        msr = ids["MSR*8"][k]
        fin = msr[np.isfinite(msr) & (msr > 0)]
        txt = (rf"MSR*8 (off scale): median $10^{{{np.log10(np.median(fin)):.0f}}}$, "
               rf"max $10^{{{np.log10(fin.max()):.0f}}}$" if fin.size else "MSR*8: overflow")
        if fin.size < msr.size:
            txt += f"\n(non-finite on {msr.size - fin.size} parallels)"
        ax.text(0.5, 0.96, txt, transform=ax.transAxes, ha="center", va="top", fontsize=6,
                color=ss.C_MSR)
        ax.set_yscale("log")
        ax.set_xlabel(r"Latitude $\varphi$ [deg]")
        ax.set_xlim(-90, 90)
        ax.set_xticks(np.arange(-90, 91, 30))
        ss.panel_title(ax, "ab"[k], titles[k])
    axes[0].set_ylabel("Relative residual [-]")
    axes[0].set_ylim(1e-20, 1e-9)
    fig.legend(*axes[0].get_legend_handles_labels(), loc="outside upper center", ncol=4)
    ss.save_figure(fig, name)


# ============================================================
# FIG. 3: INSTABILITY OF THE MSR
# ============================================================
def msr_errors(N, lat_deg):
    """log10 relative error of MSR*8 (and FCR*8, MSR*16) vs FCR*16, (L, N+1, N+1)."""
    u = np.sin(np.radians(np.asarray(lat_deg, dtype=np.float64)))
    ref = sc.helmholtz_fcr(u.astype(LD), N, LD)
    res = {}
    for name, A in (("MSR*8", sc.helmholtz_msr(u, N)),
                    ("MSR*16", sc.helmholtz_msr(u.astype(LD), N, LD)),
                    ("FCR*8", sc.helmholtz_fcr(u, N))):
        with np.errstate(all="ignore"):
            err = np.abs(A.astype(LD) - ref)
            rel = (err / np.abs(ref)).astype(np.float64)
            # per degree, norm-wise: max_m |dA_nm| / max_m |A_nm| (robust to zeros of A_nm)
            row = (err.max(axis=2) / np.abs(ref).max(axis=2)).astype(np.float64)
            res[name] = (np.log10(np.maximum(rel, 1e-20)), np.log10(np.maximum(row, 1e-20)))
    tri = np.tril(np.ones((N + 1, N + 1), bool))
    for name in res:
        res[name][0][:, ~tri] = np.nan
    return res


def fig_msr(N, map_lats, curve_lats, name="fig03_msr_instability"):
    lats = list(map_lats) + [x for x in curve_lats if x not in map_lats]
    t0 = time.perf_counter()
    res = msr_errors(N, lats)
    print(f"[accuracy] MSR vs FCR*16 on {len(lats)} parallels: {time.perf_counter() - t0:.1f} s")
    fig = plt.figure(figsize=(ss.TXT_W, 2.75), layout="constrained")
    gs = fig.add_gridspec(1, 4, width_ratios=[1, 1, 0.06, 1.25])
    maps = [fig.add_subplot(gs[0, 0])]
    maps.append(fig.add_subplot(gs[0, 1], sharex=maps[0], sharey=maps[0]))
    cax, ax3 = fig.add_subplot(gs[0, 2]), fig.add_subplot(gs[0, 3])
    vmax = max(np.nanmax(res["MSR*8"][0][lats.index(x)]) for x in map_lats)
    cmap = plt.get_cmap("viridis").copy()
    cmap.set_bad("white")
    for k, (ax, latv) in enumerate(zip(maps, map_lats)):
        E = res["MSR*8"][0][lats.index(latv)]
        im = ax.imshow(E, origin="upper", cmap=cmap, vmin=-17, vmax=np.ceil(vmax / 10) * 10,
                       extent=(-0.5, N + 0.5, N + 0.5, -0.5), interpolation="nearest",
                       rasterized=True)
        ax.set_xlabel("Order $m$ [-]")
        ax.grid(False)
        ss.panel_title(ax, "ab"[k], rf"MSR*8, $\varphi = {latv:g}^\circ$")
        ax.text(0.95, 0.95, f"max {np.nanmax(E):.0f}", transform=ax.transAxes, ha="right",
                va="top", fontsize=6.5)
    maps[0].set_ylabel("Degree $n$ [-]")
    plt.setp(maps[1].get_yticklabels(), visible=False)
    cb = fig.colorbar(im, cax=cax)
    cb.set_label(r"$\log_{10}$ relative error vs FCR*16 [-]")
    cb.outline.set_linewidth(0.5)

    n = np.arange(N + 1)
    shades = np.linspace(0.35, 1.0, len(curve_lats))
    for latv, a in zip(curve_lats, shades):
        ax3.plot(n, res["MSR*8"][1][lats.index(latv)], color=ss._tint(ss.C_MSR, a), lw=0.9,
                 label=rf"MSR*8, ${latv:g}^\circ$")
    l0 = map_lats[0]
    ax3.plot(n, res["MSR*16"][1][lats.index(l0)], color=ss.C_MSR, lw=0.9, ls=(0, (3, 1.5)),
             label=rf"MSR*16, ${l0:g}^\circ$")
    fcr = np.max(np.stack([res["FCR*8"][1][lats.index(x)] for x in curve_lats]), axis=0)
    ax3.plot(n, fcr, color=ss.C_FCR, lw=0.9, label="FCR*8, all")
    ax3.axhline(np.log10(EPS), color=ss.C_PAPER, lw=0.7, ls="--")
    ax3.set_xlabel("Degree $n$ [-]")
    ax3.set_ylabel(r"$\log_{10} \max_m |\delta\bar A_{nm}| / \max_m |\bar A_{nm}|$ [-]")
    ax3.set_xlim(0, N)
    ss.panel_title(ax3, "c", "Error per degree")
    ax3.set_ylim(-20, 75)
    ax3.legend(loc="upper left", fontsize=6, ncol=2, handlelength=1.3, columnspacing=0.8)
    ss.save_figure(fig, name)
    return {latv: float(np.nanmax(res["MSR*8"][0][lats.index(latv)])) for latv in lats}


# ============================================================
# MAIN
# ============================================================
def write_summary(lap, ids, msr_max, outdir):
    rows = ["test,variant,median,max"]
    print("\n[accuracy] Eq. 31 |trace| [s^-2] / relative to |V_rr| (median, max over points)")
    for m, (ab, rel) in lap.items():
        rows.append(f"Eq31_abs,{m},{np.median(ab):.3e},{ab.max():.3e}")
        rows.append(f"Eq31_rel,{m},{np.median(rel):.3e},{rel.max():.3e}")
        print(f"   {m:12s} {np.median(ab):9.2e} {ab.max():9.2e}   "
              f"{np.median(rel):9.2e} {rel.max():9.2e}")
    print("[accuracy] Eqs. 32 / 33 relative residual (median, max over parallels)")
    for name, (a, b) in ids.items():
        rows.append(f"Eq32,{name},{np.nanmedian(a):.3e},{np.nanmax(a):.3e}")
        rows.append(f"Eq33,{name},{np.nanmedian(b):.3e},{np.nanmax(b):.3e}")
        print(f"   {name:16s} {np.nanmedian(a):9.2e} {np.nanmax(a):9.2e}   "
              f"{np.nanmedian(b):9.2e} {np.nanmax(b):9.2e}")
    print("[accuracy] MSR*8 max log10 relative error vs FCR*16 per parallel")
    for latv, v in msr_max.items():
        rows.append(f"MSR_max_log10_relerr,lat={latv:g},,{v:.2f}")
        print(f"   {latv:5g} deg  {v:6.1f}")
    ss.write_output(outdir / "accuracy_summary.csv", "\n".join(rows) + "\n", tag="Data")


def parse(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--N", type=int, default=360, help="maximum degree (paper: 360)")
    ap.add_argument("--gfc", type=str, default=None, help="ICGEM .gfc model (default: synthetic)")
    ap.add_argument("--seed", type=int, default=2009, help="seed of the synthetic model")
    ap.add_argument("--quick", action="store_true", help="use the --quick precision setup")
    ap.add_argument("--id-step", type=float, default=1.0, help="Eqs. 32-33 latitude step [deg]")
    ap.add_argument("--msr-lats", type=float, nargs=2, default=(30.0, 60.0),
                    help="latitudes of the two (n, m) maps of Fig. 3 [deg]")
    ap.add_argument("--backend", choices=("numpy", "rust"), default="numpy",
                    help="origin of the *8 results of Eq. 31 (see sph_precision.py)")
    ap.add_argument("--out", type=str, default=None,
                    help="output folder (default ./results, or ./results_rust)")
    a = ap.parse_args(argv)
    a.out = a.out or str(HERE / ("results_rust" if a.backend == "rust" else "results"))
    return a


def main(argv=None):
    a = parse(argv)
    outdir = ss.set_figure_dir(a.out)
    ss.apply_style()
    # Eq. 31 on the precision GRID (same setup and cache as sph_precision.py)
    pa = sp.parse(["--N", str(a.N), "--seed", str(a.seed), "--out", a.out, "--backend", a.backend]
                  + (["--gfc", a.gfc] if a.gfc else []) + (["--quick"] if a.quick else []))
    model = sp.get_model(pa)
    loc, _ = sp.load_or_compute(model, pa, outdir)
    lap = laplace(loc)
    fig_laplace(loc, lap)
    # Eqs. 32-33
    step = 2.0 if a.quick and a.id_step == 1.0 else a.id_step
    lat_deg = np.arange(-89.0, 89.0 + 1e-9, step)
    ids = identities(a.N, lat_deg)
    fig_identities(lat_deg, ids)
    # Fig. 3
    msr_max = fig_msr(a.N, tuple(a.msr_lats), (5.0, 15.0, 30.0, 60.0, 85.0))
    write_summary(lap, ids, msr_max, outdir)
    return lap, ids, msr_max


if __name__ == "__main__":
    main()
