"""
sph_efficiency.py -- computational efficiency of the synthesis programmes
(Table 7 and Fig. 4 of Fantino & Casotto 2009).

Scenarios (paper, Sect. 7):
    GRID    64,440 points, 1 deg x 1 deg, h = 250 km
    ORBIT1  64,440 points of a near-polar LEO (a = 6628 km, e = 0.002, i = 96.5 deg), every 31 s
    ORBIT2  2,000,000 points at 1 Hz: extrapolated from ORBIT1 (as the CPU time per point is
            constant along an orbit)

Each programme is timed (process CPU time, single thread) on an evenly spaced subset of the
scenario and extrapolated to the full set (--full times every point). CLENSHAW on the GRID uses
the per-parallel arrangement (clenshaw_grid) on all 64,440 points. The cost per point is also
measured as a function of the maximum degree N.

Outputs (in --out, default ./results):
    table7.csv, table7.md, efficiency.json
    fig04_ground_track.pdf     ORBIT1 ground track (Fig. 4) + latitude distribution
    fig_efficiency.pdf         time relative to LEGENDRE_a vs the paper; cost per point vs N

    python sph_efficiency.py               # subset timing (~3-5 min)
    python sph_efficiency.py --full        # every point (~15-20 min)
    python sph_efficiency.py --quick       # smoke test (~1 min)
    python sph_efficiency.py --plot-only   # redraw from efficiency.json
    python sph_efficiency.py --backend rust   # compiled kernels (sph_rs/), every point (~10 min)

--backend rust runs the scalar, point-by-point Rust kernels of sph_rust.py (native binary): the
CPU time of the synthesis loop is measured inside the binary on every point of the scenario
(best of --reps), with no extrapolation.
"""

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")  # before numpy is imported: one thread, as in the paper

import argparse
import json
import time
from pathlib import Path

import numpy as np

import sph_core as sc
import sph_rust as sr
import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
from matplotlib.lines import Line2D
import sph_style as ss

HERE = Path(__file__).resolve().parent

# Table 7 of the paper: CPU seconds on a Pentium IV 1.3 GHz (Fortran; LEGENDRE_b Matlab + DLL)
PAPER_T7 = {
    "LEGENDRE_a": {"GRID": 2390, "ORBIT1": 2079, "ORBIT2": 64449},
    "LEGENDRE_b": {"GRID": 9793, "ORBIT1": 9975, "ORBIT2": 309225},
    "PINES": {"GRID": 5275, "ORBIT1": 6306, "ORBIT2": 195486},
    "METRIS": {"GRID": 3271, "ORBIT1": 3418, "ORBIT2": 105958},
    "CLENSHAW": {"GRID": 8, "ORBIT1": 992, "ORBIT2": 30752},
}
N_POINTS = {"GRID": 64440, "ORBIT1": 64440, "ORBIT2": 2_000_000}
SCENARIOS = ("GRID", "ORBIT1", "ORBIT2")


# ============================================================
# TIMING
# ============================================================
def _subset(K, n):
    return np.unique(np.linspace(0, K - 1, min(n, K)).round().astype(int))


def time_points(method, model, lam, phi, r, idx, reps=1, backend="numpy"):
    """Best-of-reps CPU and wall time per point of `method` on the points idx."""
    if backend == "rust":
        res = sr.run(method, model, lam[idx], phi[idx], r[idx], reps=reps)
        return res["cpu"] / idx.size, res["wall"] / idx.size
    w = idx[: min(4, idx.size)]
    sc.synthesize(method, model, lam[w], phi[w], r[w])  # warm-up: coefficient tables, caches
    best = (np.inf, np.inf)
    for _ in range(reps):
        c0, w0 = time.process_time(), time.perf_counter()
        sc.synthesize(method, model, lam[idx], phi[idx], r[idx])
        best = min(best, (time.process_time() - c0, time.perf_counter() - w0))
    return best[0] / idx.size, best[1] / idx.size


def time_clenshaw_grid(model, lon, lat, r, reps=3, backend="numpy"):
    if backend == "rust":
        res = sr.run_grid(model, np.radians(lon), np.radians(lat), r, reps=max(reps, 5))
        return res["cpu"] / res["K"], res["wall"] / res["K"]
    sc.clenshaw_grid(model, np.radians(lon[:2]), np.radians(lat[:2]), r)
    best = (np.inf, np.inf)
    for _ in range(reps):
        c0, w0 = time.process_time(), time.perf_counter()
        sc.clenshaw_grid(model, np.radians(lon), np.radians(lat), r)
        best = min(best, (time.process_time() - c0, time.perf_counter() - w0))
    K = lon.size * lat.size
    return best[0] / K, best[1] / K


def run_table7(model, n_sub, n_b, full=False, reps=1, backend="numpy"):
    """CPU seconds of each programme for GRID, ORBIT1 (measured/extrapolated) and ORBIT2."""
    lon, lat = sc.grid_axes()
    h = 250e3
    G = sc.grid_points(h, lon, lat, model.a)
    O = sc.orbit_points(N_POINTS["ORBIT1"], 31.0)[:3]
    res = {}
    for meth in sc.METHODS:
        res[meth] = {}
        for scen, (lam, phi, r) in (("GRID", G), ("ORBIT1", O)):
            K = lam.size
            if meth == "CLENSHAW" and scen == "GRID":
                cpu, wall = time_clenshaw_grid(model, lon, lat, model.a + h, backend=backend)
                n_used = K
            else:
                n = n_b if meth == "LEGENDRE_b" and not full else (K if full else n_sub)
                idx = _subset(K, n)
                cpu, wall = time_points(meth, model, lam, phi, r, idx, reps, backend)
                n_used = idx.size
            res[meth][scen] = {"cpu_per_point": cpu, "wall_per_point": wall, "n_timed": int(n_used),
                               "cpu": cpu * N_POINTS[scen], "wall": wall * N_POINTS[scen]}
            print(f"  {meth:11s} {scen:7s} {cpu * 1e3:8.3f} ms/pt (CPU, {n_used} pts) "
                  f"-> {cpu * N_POINTS[scen]:9.1f} s")
        o1 = res[meth]["ORBIT1"]
        res[meth]["ORBIT2"] = {k: (v * N_POINTS["ORBIT2"] / N_POINTS["ORBIT1"] if k in ("cpu", "wall")
                                   else v) for k, v in o1.items()}
    return res


def run_scaling(model, Ns, K, K_b, seed=11, reps=1, backend="numpy"):
    """CPU time per point (random points, h = 250 km) vs the maximum degree N."""
    rng = np.random.default_rng(seed)
    lam = rng.uniform(-np.pi, np.pi, K)
    phi = np.arcsin(rng.uniform(-np.sin(np.radians(89)), np.sin(np.radians(89)), K))
    r = np.full(K, model.a + 250e3)
    out = {m: [] for m in sc.METHODS}
    for N in Ns:
        mN = model.truncate(N)
        for meth in sc.METHODS:
            idx = np.arange(K_b if meth == "LEGENDRE_b" else K)
            cpu, _ = time_points(meth, mN, lam, phi, r, idx, reps, backend)
            out[meth].append(cpu)
        print(f"  N={N:3d}: " + ", ".join(f"{m} {out[m][-1] * 1e3:.3f}" for m in sc.METHODS) + " ms/pt")
    return {"N": list(Ns), "cpu_per_point": out, "K": K, "K_b": K_b}


# ============================================================
# TABLES
# ============================================================
def write_tables(res, outdir, meta):
    lines = ["programme," + ",".join(f"{s}_cpu_s,{s}_cpu_h,{s}_paper_s,{s}_ratio,{s}_paper_ratio"
                                      for s in SCENARIOS)]
    for meth in sc.METHODS:
        row = [meth]
        for s in SCENARIOS:
            t, p = res[meth][s]["cpu"], PAPER_T7.get(meth, {}).get(s)
            row += [f"{t:.2f}", f"{t / 3600:.3f}", "" if p is None else f"{p}",
                    f"{t / res['LEGENDRE_a'][s]['cpu']:.3f}",
                    "" if p is None else f"{p / PAPER_T7['LEGENDRE_a'][s]:.3f}"]
        lines.append(",".join(row))
    ss.write_output(outdir / "table7.csv", "\n".join(lines) + "\n", tag="Table")

    md = [f"# Table 7 emulation: CPU time of the synthesis programmes (N = {meta['N']})", "",
          f"Model: {meta['model']}. Machine: {meta['machine']}. "
          f"Timing: {meta['timing']}.", "",
          "| Programme | GRID [s] | paper | ORBIT1 [s] | paper | ORBIT2 [h] | paper [h] "
          "| GRID / LEG_a | paper | ORBIT1 / LEG_a | paper |",
          "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|"]
    P0 = PAPER_T7["LEGENDRE_a"]
    for meth in sc.METHODS:
        r_, p_ = res[meth], PAPER_T7.get(meth)
        pg, po, p2, rg, ro = (("-",) * 5 if p_ is None else
                              (p_["GRID"], p_["ORBIT1"], f"{p_['ORBIT2'] / 3600:.1f}",
                               f"{p_['GRID'] / P0['GRID']:.3f}", f"{p_['ORBIT1'] / P0['ORBIT1']:.2f}"))
        md.append(
            f"| {meth} | {r_['GRID']['cpu']:.1f} | {pg} | {r_['ORBIT1']['cpu']:.1f} | "
            f"{po} | {r_['ORBIT2']['cpu'] / 3600:.2f} | {p2} | "
            f"{r_['GRID']['cpu'] / res['LEGENDRE_a']['GRID']['cpu']:.3f} | {rg} | "
            f"{r_['ORBIT1']['cpu'] / res['LEGENDRE_a']['ORBIT1']['cpu']:.2f} | {ro} |")
    md += ["", "Paper: Pentium IV 1.3 GHz, Fortran90/77 (LEGENDRE_b: Matlab + DLL), binary I/O "
               "included. ORBIT2 = ORBIT1 x 2e6 / 64,440 (paper: x 31.04). GOTTLIEB is not "
               "in the paper (added for comparison)."]
    ss.write_output(outdir / "table7.md", "\n".join(md) + "\n", tag="Table")
    print("\n".join(md))


# ============================================================
# FIGURES
# ============================================================
def fig_ground_track(name="fig04_ground_track"):
    lam, phi, r, t = sc.orbit_points(N_POINTS["ORBIT1"], 31.0)
    lon, lat = np.degrees(lam), np.degrees(phi)
    fig, (ax, axh) = plt.subplots(1, 2, figsize=(ss.TXT_W, 2.5), layout="constrained",
                                  gridspec_kw={"width_ratios": [3.2, 1.0]}, sharey=True)
    sca = ax.scatter(lon, lat, c=t / 86400.0, s=0.15, lw=0, cmap="cividis", rasterized=True)
    cb = fig.colorbar(sca, ax=ax, pad=0.01, aspect=25)
    cb.set_label("Time [d]")
    cb.ax.tick_params(labelsize=7)
    ax.set_xlim(-180, 180)
    ax.set_ylim(-90, 90)
    ax.set_xticks(np.arange(-180, 181, 60))
    ax.set_yticks(np.arange(-90, 91, 30))
    ax.set_xlabel(r"Longitude $\lambda$ [deg]")
    ax.set_ylabel(r"Latitude $\varphi$ [deg]")
    ss.panel_title(ax, "a", f"ORBIT1 ground track ({lon.size:,} points, $\\Delta t$ = 31 s)")
    bins = np.arange(-90, 91, 5)
    lon_g, lat_g = sc.grid_axes()
    hg, _ = np.histogram(np.repeat(lat_g, lon_g.size), bins)
    ho, _ = np.histogram(lat, bins)
    axh.stairs(100 * ho / ho.sum(), bins, orientation="horizontal", fill=True,
               color=ss._tint(ss.METHOD_COLOR["PINES"], 0.35), label="ORBIT1")
    axh.stairs(100 * ho / ho.sum(), bins, orientation="horizontal", color=ss.METHOD_COLOR["PINES"],
               lw=0.8)
    axh.stairs(100 * hg / hg.sum(), bins, orientation="horizontal", color="0.25", lw=0.9, ls="--",
               label="GRID")
    axh.set_xlabel("Points per 5 deg [\\%]" if ss.USE_TEX else "Points per 5 deg [%]")
    axh.legend(loc="center right", fontsize=6)
    ss.panel_title(axh, "b", "Latitude distribution")
    ss.save_figure(fig, name)


def _dumbbell(ax, res, scen):
    ref_o, ref_p = res["LEGENDRE_a"][scen]["cpu"], PAPER_T7["LEGENDRE_a"][scen]
    meths = list(sc.METHODS)
    y = np.arange(len(meths))[::-1]
    vals = []
    for yi, meth in zip(y, meths):
        ours = res[meth][scen]["cpu"] / ref_o
        vals.append(ours)
        if meth in PAPER_T7:  # GOTTLIEB is not in the paper
            pap = PAPER_T7[meth][scen] / ref_p
            vals.append(pap)
            ax.plot([ours, pap], [yi, yi], color="0.75", lw=1.0, zorder=1)
            ax.plot(pap, yi, marker=ss.METHOD_MARKER[meth], mfc="white", mec=ss.C_PAPER, ms=5,
                    ls="none", zorder=2)
        ax.plot(ours, yi, marker=ss.METHOD_MARKER[meth], color=ss.METHOD_COLOR[meth], ms=5,
                ls="none", zorder=3)
        tr = mtransforms.blended_transform_factory(ax.transAxes, ax.transData)
        t = res[meth][scen]["cpu"]
        txt = f"{t:.2f} s" if t < 1 else (f"{t:.1f} s" if t < 100 else f"{t:.0f} s")
        ax.text(0.99, yi + 0.22, txt, transform=tr, ha="right", va="bottom", fontsize=6,
                color=ss.METHOD_COLOR[meth])
    ax.axvline(1.0, color="0.4", lw=0.6, ls=":", zorder=0)
    ax.set_xscale("log")
    lo, hi = min(vals), max(vals)
    ax.set_xlim(lo / 2.5, hi * 12)
    ax.set_yticks(y)
    ax.set_yticklabels([ss.method_label(m) for m in meths])
    ax.set_ylim(-0.6, len(meths) - 0.3)
    ax.grid(axis="y", visible=False)
    ax.set_xlabel(r"CPU time / LEGENDRE\_a [-]" if ss.USE_TEX else "CPU time / LEGENDRE_a [-]")


def fig_efficiency(res, scal, meta, name="fig_efficiency"):
    fig, axs = plt.subplots(1, 3, figsize=(ss.TXT_W, 2.55), layout="constrained",
                            gridspec_kw={"width_ratios": [1, 1, 1.05]})
    _dumbbell(axs[0], res, "GRID")
    ss.panel_title(axs[0], "a", f"GRID ({N_POINTS['GRID']:,} points)")
    _dumbbell(axs[1], res, "ORBIT1")
    axs[1].set_yticklabels([])
    ss.panel_title(axs[1], "b", f"ORBIT1 ({N_POINTS['ORBIT1']:,} points)")
    ax = axs[2]
    Ns = np.asarray(scal["N"], float)
    for meth in sc.METHODS:
        ax.plot(Ns, 1e3 * np.asarray(scal["cpu_per_point"][meth]), marker=ss.METHOD_MARKER[meth],
                ms=3, color=ss.METHOD_COLOR[meth], lw=1.0, label=ss.method_label(meth))
    t_ref = 1e3 * scal["cpu_per_point"]["LEGENDRE_a"][-1]
    ax.plot(Ns, t_ref * (Ns / Ns[-1]) ** 2, color="0.5", lw=0.8, ls="--")
    i_t = max(1, len(Ns) // 4)
    ax.text(Ns[i_t] * 1.08, t_ref * (Ns[i_t] / Ns[-1]) ** 2 * 0.55, r"$\propto N^2$",
            fontsize=7, color="0.4", ha="left", va="top")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ticks = [t for t in (20, 30, 50, 100, 200, 360) if Ns[0] <= t <= Ns[-1]]
    ax.set_xticks(ticks)
    ax.set_xticklabels([str(t) for t in ticks])
    ax.xaxis.set_minor_formatter(plt.NullFormatter())
    ax.set_xlabel(r"Maximum degree $N$ [-]")
    ax.set_ylabel("CPU time per point [ms]")
    ss.panel_title(ax, "c", "Cost per point vs degree")
    handles = [Line2D([], [], marker="o", color="0.15", ls="none", ms=4.5,
                      label=f"This emulation ({meta.get('backend_label', 'NumPy')}, N = {meta['N']})"),
               Line2D([], [], marker="o", mfc="white", mec=ss.C_PAPER, ls="none", ms=4.5,
                      label="Paper, Table 7 (Fortran, Pentium IV)")]
    handles += [Line2D([], [], marker=ss.METHOD_MARKER[m], color=ss.METHOD_COLOR[m], ms=3.5,
                       lw=1.0, label=ss.method_label(m)) for m in sc.METHODS]
    fig.legend(handles=handles, loc="outside upper center", ncol=4, fontsize=6.5)
    ss.save_figure(fig, name)


# ============================================================
# MAIN
# ============================================================
def machine_string(backend="numpy"):
    import platform

    if backend == "rust":
        return f"{platform.system()} {platform.release()}, {sr.toolchain_string()}, " \
               f"scalar kernels (no fast-math), single thread"
    return f"{platform.machine()} {platform.system()} {platform.release()}, Python " \
           f"{platform.python_version()}, NumPy {np.__version__}, single thread"


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--N", type=int, default=360, help="maximum degree (paper: 360)")
    ap.add_argument("--gfc", type=str, default=None, help="ICGEM .gfc model (default: synthetic)")
    ap.add_argument("--subset", type=int, default=6000, help="points timed per scenario")
    ap.add_argument("--nb", type=int, default=600, help="points timed for LEGENDRE_b")
    ap.add_argument("--full", action="store_true", help="time every point (no extrapolation)")
    ap.add_argument("--quick", action="store_true", help="small smoke-test run")
    ap.add_argument("--plot-only", action="store_true", help="redraw from efficiency.json")
    ap.add_argument("--backend", choices=("numpy", "rust"), default="numpy",
                    help="numpy (vectorised) or rust (compiled scalar kernels, sph_rs/)")
    ap.add_argument("--reps", type=int, default=None,
                    help="repetitions, best is kept (default: 1 numpy, 2 rust)")
    ap.add_argument("--out", type=str, default=None,
                    help="output folder (default ./results, or ./results_rust)")
    a = ap.parse_args(argv)
    rust = a.backend == "rust"
    a.out = a.out or str(HERE / ("results_rust" if rust else "results"))
    a.reps = a.reps or (2 if rust else 1)
    if a.quick:
        a.subset, a.nb = 1024, 100
    if rust and not a.quick:
        a.full = True  # compiled: every point is timed
    if rust and not a.plot_only:
        sr.ensure_built()
    outdir = ss.set_figure_dir(a.out)
    ss.apply_style()
    jpath = outdir / "efficiency.json"

    if a.plot_only:
        data = json.loads(jpath.read_text())
    else:
        model = (sc.GravityModel.from_gfc(a.gfc, a.N) if a.gfc else sc.GravityModel.synthetic(a.N))
        timing = ("every point" if a.full else
                  f"process CPU time on {a.subset} evenly spaced points per scenario "
                  f"({a.nb} for LEGENDRE_b), extrapolated; CLENSHAW/GRID per parallel on all points")
        if rust:
            timing = ("CPU time of the synthesis loop inside the binary, " + timing
                      + f", best of {a.reps}")
        meta = {"N": model.N, "model": model.name, "machine": machine_string(a.backend),
                "timing": timing, "backend": a.backend,
                "backend_label": "Rust" if rust else "NumPy"}
        print(f"[efficiency] {meta['model']} | {meta['machine']}")
        print("[efficiency] Table 7 ...")
        res = run_table7(model, a.subset, a.nb, a.full, a.reps, a.backend)
        Ns = [30, 90, 180, 360] if a.quick else [20, 30, 45, 60, 90, 120, 180, 240, 300, 360]
        Ns = [n for n in Ns if n <= model.N]
        print("[efficiency] cost vs degree ...")
        if rust:
            scal = run_scaling(model, Ns, 2000, 2000, reps=3, backend="rust")
        else:
            scal = run_scaling(model, Ns, 256 if a.quick else 512, 16 if a.quick else 48)
        data = {"meta": meta, "table7": res, "scaling": scal}
        ss.write_output(jpath, json.dumps(data, indent=1), tag="Data")
    write_tables(data["table7"], outdir, data["meta"])
    fig_ground_track()
    fig_efficiency(data["table7"], data["scaling"], data["meta"])


if __name__ == "__main__":
    main()
