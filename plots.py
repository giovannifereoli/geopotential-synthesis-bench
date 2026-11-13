"""
plots.py -- the figures, drawn only from the saved result tables (no Rust, no rerun).

    python plots.py              # results/accuracy.csv, timing.csv -> results/fig_*.pdf
    python plots.py my_results   # another results folder (e.g. one written with benchmark.py --out)

benchmark.py calls this at the end of every run. Run it on its own to redraw the figures later.
Text is set with LaTeX when a `latex` binary is found, otherwise with matplotlib's Computer Modern.

x axis = maximum degree N (log scale) everywhere.

One panel per method (blue solid line = acceleration, orange dashed line = gravity gradient for the
accuracy, acceleration + gravity gradient for the CPU time):
  fig_accuracy.pdf      max. relative error over all test points (all three point sets)
  fig_speed.pdf         (a) CPU time per point of LEGENDRE_a in microseconds;
                        (b) CPU time of every other method divided by that of LEGENDRE_a
All methods on one axis, one panel per output (one colour and marker per method), for direct comparison:
  fig_accuracy_all.pdf  max. relative error, acceleration and gravity gradient
  fig_speed_all.pdf     CPU time divided by that of LEGENDRE_a, acceleration and acceleration + gradient;
                        each curve is labelled at its right end
CPU times: the cost of a point does not depend on where the point is, so the three point sets are three
repeat measurements of the same quantity. Curves show their median, shaded bands / bars their range.
Only when the folder holds simd.csv (benchmark.py --phase simd):
  fig_simd.pdf          speed-up from LLVM auto-vectorization, one panel per method
A folder timed with the auto-vectorized build (meta.json) says so in the titles of the CPU-time figures.
"""

import json
import math
import shutil
import statistics
import sys
from pathlib import Path

from benchmark import GATE_ABS, LEVELS, METHODS, read_csv

REF = "LEGENDRE_a"                                    # the speed reference
# Fixed axis ranges for the speed figures, the same for every build and run, so that figures
# of different runs can be compared by eye (widened only if some data falls outside).
RATIO_LIM = (0.58, 2.8)                               # CPU time / CPU time of REF
ABS_LIM = (0.1, 2e4)                                  # CPU time per point of REF [us]
MUTED, LIGHT, BAND = "#52514e", "#e1e0d9", "#ecebe6"
STYLE = {1: dict(color="#2a78d6", ls="-", marker="o", ms=3.2, lw=1.4, label="Acceleration"),
         2: dict(color="#eb6834", ls="--", marker="s", ms=3.0, lw=1.4, label="Acceleration + gravity gradient")}
OVF = dict(ls="none", marker="x", ms=4.5, mew=1.2, color=MUTED)
INK = "#0b0b0b"
# fig_compare: one colour + marker per method (validated categorical palette, fixed order)
CMP = [("LEAR", "#2a78d6", "o"), ("METRIS", "#eb6834", "s"), ("CLENSHAW", "#1baf7a", "^"),
       ("PINES", "#eda100", "D"), ("GOTTLIEB", "#e87ba4", "v"), ("LEGENDRE_b", "#008300", "P"),
       ("HOLMES", "#6250d6", "h"), ("CUNNINGHAM", "#e34948", "*")]


def msz(mk):
    """Marker size; the star reads smaller than the other shapes at equal size."""
    return 6.0 if mk == "*" else 3.8


N_TICKS = (10, 100, 1000)


def setup():
    """LaTeX fonts; returns a function that typesets a method name."""
    import logging
    import matplotlib
    matplotlib.use("Agg")
    logging.getLogger("fontTools").setLevel(logging.ERROR)
    import matplotlib.pyplot as plt
    tex = shutil.which("latex") is not None
    plt.rcParams.update({
        "text.usetex": tex, "font.family": "serif",
        "font.serif": ["Computer Modern Roman"] if tex else ["cmr10", "DejaVu Serif"],
        "mathtext.fontset": "cm", "axes.formatter.use_mathtext": True, "axes.unicode_minus": False,
        "font.size": 8, "axes.titlesize": 8, "axes.labelsize": 8, "legend.fontsize": 7.5,
        "xtick.labelsize": 7, "ytick.labelsize": 7, "axes.linewidth": 0.6,
        "xtick.major.width": 0.6, "ytick.major.width": 0.6, "xtick.minor.width": 0.4,
        "ytick.minor.width": 0.4, "xtick.direction": "in", "ytick.direction": "in",
        "pdf.fonttype": 42})
    if tex:
        return lambda m: r"\texttt{" + m.replace("_", r"\_") + "}"
    return lambda m: m


def log_axis(axis, ticks, fmt=str):
    from matplotlib.ticker import FixedLocator, NullFormatter
    axis.set_major_locator(FixedLocator(ticks))
    axis.set_major_formatter(lambda v, _: fmt(v))
    axis.set_minor_formatter(NullFormatter())


def frame(ax, Ns):
    ax.set_xscale("log")
    ax.set_xlim(Ns[0] / 1.3, Ns[-1] * 1.3)
    log_axis(ax.xaxis, [t for t in N_TICKS if Ns[0] <= t <= Ns[-1]], lambda v: f"{v:.0f}")
    ax.grid(True, which="major", color=LIGHT, lw=0.5)
    ax.tick_params(which="both", top=True, right=True)
    ax.set_axisbelow(True)


def spread(ys, gap):
    """Label positions (axes fractions) moved apart by at least `gap`, order kept."""
    idx = sorted(range(len(ys)), key=ys.__getitem__)
    p = [ys[i] for i in idx]
    for _ in range(500):
        moved = False
        for k in range(1, len(p)):
            d = p[k] - p[k - 1]
            if d < gap - 1e-9:
                p[k - 1] -= (gap - d) / 2
                p[k] += (gap - d) / 2
                moved = True
        if not moved:
            break
    out = [0.0] * len(ys)
    for k, i in enumerate(idx):
        out[i] = p[k]
    return out


def pooled(trows):
    """{(method, N, level): dict(t, t_lo, t_hi, r, r_lo, r_hi)}: CPU time per point [us] and its ratio to
    LEGENDRE_a timed on the same point set, as median and range over the point sets that pass the gate."""
    by = {}
    for r in trows:
        by.setdefault((r["method"], r["N"], r["level"]), {})[r["set"]] = r
    out = {}
    for (m, N, lev), d in by.items():
        ref = by.get((REF, N, lev), {})
        ts = [r["us_per_pt"] for r in d.values() if r["gate_pass"]]
        rs = [r["us_per_pt"] / ref[s]["us_per_pt"] for s, r in d.items() if r["gate_pass"] and s in ref]
        if ts and rs:
            out[(m, N, lev)] = dict(t=statistics.median(ts), t_lo=min(ts), t_hi=max(ts),
                                    r=statistics.median(rs), r_lo=min(rs), r_hi=max(rs))
    return out


def sci(v):
    """2e-12 -> $2 \\times 10^{-12}$"""
    e = math.floor(math.log10(v) + 1e-9)
    return f"${v / 10 ** e:g}" + r" \times 10^{" + f"{e}" + "}$"


def overflow_note(missing):
    """missing: {method: first N not shown} -> 'PINES and GOTTLIEB from N = 1500, CLENSHAW at N = 2000'."""
    by_n = {}
    for m, n in missing.items():
        by_n.setdefault(n, []).append(m)
    return ", ".join(" and ".join(ms) + f" from $N = {n}$" for n, ms in sorted(by_n.items()))


def speed_all(folder, t, Ns, tt, note=""):
    """fig_speed_all.pdf: all methods on one axis per output, CPU time relative to LEGENDRE_a
    (median over the three point sets; vertical bars = their range)."""
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    from matplotlib.transforms import blended_transform_factory
    fig, axs = plt.subplots(1, 2, figsize=(7.16, 4.0), layout="constrained")
    missing = {}
    for ax, lev, title in zip(axs, LEVELS, ("(a) Acceleration", "(b) Acceleration + gravity gradient")):
        frame(ax, Ns)
        ax.set_yscale("log")
        ax.axhspan(1 / 1.05, 1.05, color=BAND, lw=0, zorder=0)
        ax.axhline(1, color=INK, lw=1.1, zorder=2)
        ends, vals = [(REF, None, Ns[-1], 1.0, INK)], [1.0]
        for m, c, mk in CMP:
            q = [(N, t[(m, N, lev)]) for N in Ns if (m, N, lev) in t]
            if not q:
                continue
            pts = [(N, d["r"]) for N, d in q]
            if pts[-1][0] != Ns[-1]:
                missing[m] = Ns[Ns.index(pts[-1][0]) + 1]
            ax.vlines([N for N, _ in q], [d["r_lo"] for _, d in q], [d["r_hi"] for _, d in q],
                      color=c, lw=0.9, alpha=0.55, zorder=2.5)
            ax.plot(*zip(*pts), color=c, marker=mk, ms=msz(mk), lw=1.4, mec="white", mew=0.5, zorder=3)
            ends.append((m, mk, *pts[-1], c))
            vals += [d["r_lo"] for _, d in q] + [d["r_hi"] for _, d in q]
        ylo, yhi = min(RATIO_LIM[0], min(vals) / 1.05), max(RATIO_LIM[1], max(vals) * 1.05)
        ax.set_ylim(ylo, yhi)
        log_axis(ax.yaxis, [v for v in (0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.25, 1.5, 2, 2.5, 3, 4) if ylo <= v <= yhi],
                 lambda v: f"{v:g}")
        # direct labels at the right: marker + name, joined to the end of each curve by a thin leader
        frac = lambda y: (math.log(y) - math.log(ylo)) / (math.log(yhi) - math.log(ylo))
        pos = spread([frac(e[3]) for e in ends], 0.072)
        tr = blended_transform_factory(ax.transData, ax.transAxes)
        xm = Ns[-1] * 1.9
        for (m, mk, N, y, c), p in zip(ends, pos):
            ax.annotate("", xy=(N, y), xytext=(xm, p), textcoords=tr, annotation_clip=False,
                        arrowprops=dict(arrowstyle="-", color=c, lw=0.6, ls=":",
                                        shrinkA=3, shrinkB=0))
            if mk:
                ax.plot([xm], [p], transform=tr, marker=mk, color=c, ms=msz(mk), mec="white", mew=0.5, clip_on=False)
            ax.text(xm * 1.22, p, tt(m), transform=tr, va="center", ha="left", fontsize=7, color=INK)
        ax.set_title(title, loc="left")
        ax.set_xlabel("Maximum degree $N$")
    axs[0].set_ylabel("CPU time / CPU time of " + tt(REF))
    usetex = plt.rcParams["text.usetex"]
    handles = [Patch(fc=BAND, ec="none", label=(r"Within 5\% of " if usetex else "Within 5% of ") + tt(REF)
                     + " (timing noise)"),
               Line2D([], [], color=MUTED, lw=0.9, alpha=0.7,
                      label="Vertical bars: range over the three point sets")]
    if missing:
        handles.append(Line2D([], [], ls="none", label="Curves stop where a method overflows: " + overflow_note(missing)))
    fig.legend(handles=handles, loc="outside lower center", ncol=1, frameon=False, handlelength=1.6,
               labelspacing=0.3)
    fig.suptitle("CPU time per point relative to " + tt(REF) + " (median over the three point sets; below 1 = faster)"
                 + note, fontsize=8)
    fig.savefig(folder / "fig_speed_all.pdf")
    plt.close(fig)


def simd_fig(folder, tt):
    """fig_simd.pdf: CPU time of the scalar build / CPU time of the auto-vectorized build, per method
    (simd.csv: both builds timed side by side on ORBIT1; band = 95% bootstrap interval)."""
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    rows = read_csv(folder / "simd.csv")
    Ns = sorted({r["N"] for r in rows})
    ylo = min(0.8, min(r["ci_lo"] for r in rows) / 1.05)
    yhi = max(1.25, max(r["ci_hi"] for r in rows) * 1.12)
    fig = plt.figure(figsize=(7.16, 5.2), layout="constrained")
    grid = fig.subfigures(1, 1)
    axs = grid.subplots(3, 3, sharex=True, sharey=True)
    for ax, m in zip(axs.flat, METHODS):
        frame(ax, Ns)
        ax.set_yscale("log")
        ax.set_ylim(ylo, yhi)
        log_axis(ax.yaxis, [v for v in (0.8, 0.9, 1, 1.25, 1.5, 2, 2.5, 3) if ylo <= v <= yhi], lambda v: f"{v:g}")
        ax.axhline(1, color=MUTED, lw=0.6)
        for lev in LEVELS:
            q = sorted((r for r in rows if r["method"] == m and r["level"] == lev), key=lambda r: r["N"])
            if q:
                ax.fill_between([r["N"] for r in q], [r["ci_lo"] for r in q], [r["ci_hi"] for r in q],
                                color=STYLE[lev]["color"], alpha=0.2, lw=0)
                ax.plot([r["N"] for r in q], [r["speedup"] for r in q], **STYLE[lev])
        bad = [N for N in Ns if not any(r["method"] == m and r["N"] == N for r in rows)]
        if bad:
            ax.plot(bad, [yhi / 1.09] * len(bad), **OVF)
        ax.set_title(tt(m))
    fig.legend(handles=[Line2D([], [], **STYLE[1]), Line2D([], [], **STYLE[2]),
                        Patch(fc=MUTED, alpha=0.25, ec="none", label=r"95\% interval" if plt.rcParams["text.usetex"]
                              else "95% interval"),
                        Line2D([], [], **{**OVF, "label": "Overflow, not timed"})],
               loc="outside lower center", ncol=4, frameon=False, handlelength=2.6, columnspacing=1.6)
    grid.supxlabel("Maximum degree $N$", fontsize=8)
    grid.supylabel("CPU time without / with auto-vectorization", fontsize=8)
    fig.suptitle("Speed-up from LLVM auto-vectorization (above 1 = faster with it, below 1 = slower; "
                 "both builds timed side by side)", fontsize=8)
    fig.savefig(folder / "fig_simd.pdf")
    plt.close(fig)


def accuracy_all(folder, gate, Ns, tt):
    """fig_accuracy_all.pdf: all methods on one axis per output, max. relative error over all test points."""
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    fig, axs = plt.subplots(1, 2, figsize=(7.16, 3.7), sharey=True, layout="constrained")
    series = [(REF, INK, "X")] + CMP
    missing, handles = {}, []
    for ax, key, title in zip(axs, ("acc", "tensor"), ("(a) Acceleration", "(b) Gravity gradient")):
        frame(ax, Ns)
        ax.set_yscale("log")
        drawn = []                                    # [(points, [names], colour, marker)]
        for m, c, mk in series:
            pts = []
            for N in Ns:
                rr = [r for r in gate if r["N"] == N and r["method"] == m]
                if rr and all(r["finite"] for r in rr):
                    pts.append((N, max(r[key] for r in rr)))
            if not pts:
                continue
            if pts[-1][0] != Ns[-1]:
                missing[m] = Ns[Ns.index(pts[-1][0]) + 1]
            same = [d for d in drawn if d[0] == pts]
            if same:                                  # identical results (LEGENDRE_a and LEGENDRE_b): one curve
                same[0][1].append(m)
                continue
            drawn.append((pts, [m], c, mk))
        for pts, names, c, mk in drawn:
            ax.plot(*zip(*pts), color=c, marker=mk, ms=msz(mk), lw=1.3, mec="white", mew=0.5, zorder=3)
        ax.set_title(title, loc="left")
        ax.set_xlabel("Maximum degree $N$")
        if not handles:
            handles = [Line2D([], [], color=c, marker=mk, ms=msz(mk), lw=1.3, mec="white", mew=0.5,
                              label=", ".join(tt(n) for n in names) + (" (identical)" if len(names) > 1 else ""))
                       for _, names, c, mk in drawn]
    lo, hi = axs[0].get_ylim()
    log_axis(axs[0].yaxis, [v for v in (1e-13, 2e-13, 5e-13, 1e-12, 2e-12, 5e-12, 1e-11) if lo <= v <= hi], sci)
    axs[0].set_ylabel("Max. relative error")
    fig.legend(handles=handles, loc="outside lower center", ncol=3, frameon=False, handlelength=2.4,
               columnspacing=2.0)
    title = r"Max. relative error over all test points (lower = more accurate); all are below the pass threshold $10^{-11}$"
    if missing:
        title += "\nCurves stop where a method overflows: " + overflow_note(missing)
    fig.suptitle(title, fontsize=8)
    fig.savefig(folder / "fig_accuracy_all.pdf")
    plt.close(fig)


def figures(folder):
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    tt = setup()
    folder = Path(folder)
    gate = read_csv(folder / "accuracy.csv")
    trows = read_csv(folder / "timing.csv")
    Ns = sorted({r["N"] for r in gate})
    meta = json.loads((folder / "meta.json").read_text()) if (folder / "meta.json").exists() else {}
    note = "; auto-vectorized build" if meta.get("build", {}).get("name") == "simd" else ""

    # ---- accuracy: worst case over all points of the three sets; overflow = Inf/NaN in any output
    fig = plt.figure(figsize=(7.16, 5.4), layout="constrained")
    grid = fig.subfigures(1, 1)                       # own subfigure, so its x label does not collide with the legend
    axs = grid.subplots(3, 3, sharex=True, sharey=True)
    top = GATE_ABS * 4
    for ax, m in zip(axs.flat, METHODS):
        frame(ax, Ns)
        ax.set_yscale("log")
        ax.set_ylim(1e-13, GATE_ABS * 10)
        ax.axhline(GATE_ABS, color=MUTED, lw=0.8, ls=":")
        for lev, key in zip(LEVELS, ("acc", "tensor")):
            pts = []
            for N in Ns:
                rr = [r for r in gate if r["N"] == N and r["method"] == m]
                if all(r["finite"] for r in rr):
                    pts.append((N, max(r[key] for r in rr)))
            ax.plot(*zip(*pts), **STYLE[lev])
        bad = [N for N in Ns if not all(r["finite"] for r in gate if r["N"] == N and r["method"] == m)]
        if bad:
            ax.plot(bad, [top] * len(bad), **OVF)
        ax.set_title(tt(m))
    fig.legend(handles=[Line2D([], [], **STYLE[1]),
                        Line2D([], [], **{**STYLE[2], "label": "Gravity gradient"}),
                        Line2D([], [], color=MUTED, lw=0.8, ls=":", label=r"Pass threshold $10^{-11}$"),
                        Line2D([], [], **{**OVF, "label": "Overflow (Inf or NaN)"})],
               loc="outside lower center", ncol=4, frameon=False, handlelength=2.2, columnspacing=1.8)
    grid.supxlabel("Maximum degree $N$", fontsize=8)
    grid.supylabel("Max. relative error", fontsize=8)
    fig.savefig(folder / "fig_accuracy.pdf")
    plt.close(fig)

    # ---- CPU time per point: median over the three point sets, band = their range
    t = pooled(trows)
    fig = plt.figure(figsize=(7.16, 4.1), layout="constrained")
    left, right = fig.subfigures(1, 2, width_ratios=(1.15, 3.6), wspace=0.04)

    ax = left.subplots()
    frame(ax, Ns)
    ax.set_yscale("log")
    for lev in LEVELS:
        q = [(N, t[(REF, N, lev)]) for N in Ns if (REF, N, lev) in t]
        ax.fill_between([N for N, _ in q], [d["t_lo"] for _, d in q], [d["t_hi"] for _, d in q],
                        color=STYLE[lev]["color"], alpha=0.2, lw=0)
        ax.plot([N for N, _ in q], [d["t"] for _, d in q], **STYLE[lev])
    lo = min(ABS_LIM[0], min(d["t_lo"] for (m, _, _), d in t.items() if m == REF) / 1.2)
    hi = max(ABS_LIM[1], max(d["t_hi"] for (m, _, _), d in t.items() if m == REF) * 1.2)
    ax.set_ylim(lo, hi)
    log_axis(ax.yaxis, [v for v in (0.1, 1, 10, 100, 1000, 10000) if lo <= v <= hi], lambda v: f"{v:g}")
    ax.set_xlabel("Maximum degree $N$")
    ax.set_ylabel(r"CPU time per point [$\mu$s]")
    ax.set_title("(a) " + tt(REF), loc="left")

    others = [m for m in METHODS if m != REF]
    axs = right.subplots(2, (len(others) + 1) // 2, sharex=True, sharey=True)
    for ax in axs.flat[len(others):]:                 # empty cell: show the x tick labels of the panel above
        ax.axis("off")
        axs[0, -1].xaxis.set_tick_params(labelbottom=True)
    ylo = min(RATIO_LIM[0], min(d["r_lo"] for (m, _, _), d in t.items() if m in others) / 1.06)
    yhi = max(RATIO_LIM[1], max(d["r_hi"] for (m, _, _), d in t.items() if m in others) * 1.12)
    for ax, m in zip(axs.flat, others):
        frame(ax, Ns)
        ax.set_yscale("log")
        ax.set_ylim(ylo, yhi)
        log_axis(ax.yaxis, [v for v in (0.6, 0.8, 1, 1.25, 1.5, 2, 2.5, 3, 4) if ylo <= v <= yhi],
                 lambda v: f"{v:g}")
        ax.axhspan(1 / 1.05, 1.05, color=BAND, lw=0, zorder=0)
        ax.axhline(1, color=MUTED, lw=0.6)
        for lev in LEVELS:
            q = [(N, t[(m, N, lev)]) for N in Ns if (m, N, lev) in t]
            if q:
                ax.fill_between([N for N, _ in q], [d["r_lo"] for _, d in q], [d["r_hi"] for _, d in q],
                                color=STYLE[lev]["color"], alpha=0.2, lw=0)
                ax.plot([N for N, _ in q], [d["r"] for _, d in q], **STYLE[lev])
        bad = [N for N in Ns if any((m, N, lev) not in t for lev in LEVELS)]
        if bad:
            ax.plot(bad, [yhi / 1.09] * len(bad), **OVF)
        ax.set_title(tt(m))
    right.supxlabel("Maximum degree $N$", fontsize=8)
    right.supylabel("CPU time / CPU time of " + tt(REF), fontsize=8)
    right.suptitle("(b) Other methods relative to " + tt(REF) + " (below 1 = faster)" + note, fontsize=8,
                   x=0.02, ha="left")
    fig.legend(handles=[Line2D([], [], **STYLE[1]), Line2D([], [], **STYLE[2]),
                        Patch(fc=BAND, ec="none", label=r"Within 5\% of " + tt(REF) if plt.rcParams["text.usetex"]
                              else "Within 5% of " + REF),
                        Patch(fc=MUTED, alpha=0.25, ec="none", label="Range over the three point sets"),
                        Line2D([], [], **{**OVF, "label": "Overflow, not timed"})],
               loc="outside lower center", ncol=3, frameon=False, handlelength=2.6, columnspacing=1.6)
    fig.savefig(folder / "fig_speed.pdf")
    plt.close(fig)

    speed_all(folder, t, Ns, tt, note)
    accuracy_all(folder, gate, Ns, tt)
    if (folder / "simd.csv").exists():
        simd_fig(folder, tt)


if __name__ == "__main__":
    figures(sys.argv[1] if len(sys.argv) > 1 else Path(__file__).resolve().parent / "results")
