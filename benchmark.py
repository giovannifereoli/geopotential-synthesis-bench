"""
benchmark.py -- accuracy and CPU time of eight spherical-harmonic gravity synthesis methods,
all implemented as scalar Rust kernels in sph_rs/ (built on first use with cargo).

    python benchmark.py                  # gate + timing + analysis -> results/
    python benchmark.py --quick          # smoke test (N = 10, 50, 360; short blocks) -> results_quick/
    python benchmark.py --phase analyze  # recompute tables and figures from the saved raw data in results/
    python plots.py                      # redraw only the figures from the saved tables in results/

Phases
  gate     every method vs a double-double reference (sph_rs ref) at the same points:
           max. relative error per (N, point set), pass/fail per level.
  time     R rounds over all cells (N, level, point set), each round in a new random order;
           per round one fresh process per method, methods in random order; inside each:
           warm-up + B blocks of >= T ms, thread CPU time, a calibration loop before and
           after every block. A/A probe: a second compiled copy of LEGENDRE_a.
  analyze  calibration-normalised medians, bootstrap CIs, winner rule ->
           accuracy.csv, timing.csv, winners.csv, summary.md; figures by plots.py.
"""

import argparse
import csv
import gzip
import json
import os
import platform
import shutil
import subprocess
import tempfile
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
CRATE = HERE / "sph_rs"
TARGET = Path(os.environ.get("SPH_RS_TARGET", Path.home() / ".cache" / "sph_rs" / "target"))
EXE = TARGET / "release" / "sph_rs"

METHODS = ("LEGENDRE_a", "LEGENDRE_b", "PINES", "GOTTLIEB", "METRIS", "CLENSHAW", "LEAR", "HOLMES")
AA = "LEGENDRE_a#2"             # A/A probe: second compiled copy of LEGENDRE_a
N_LIST = (10, 20, 50, 100, 200, 360, 720, 1000, 1500, 2000)
SETS = ("ORBIT1", "UNIFORM", "POLAR")
LEVELS = (1, 2)
GATE_ABS = 1e-11                # pass: max. relative error <= GATE_ABS ...
GATE_REL = 10.0                 # ... and <= GATE_REL x the best method's max. relative error
DELTA_MIN = 0.05                # smallest speed difference that can be called a win
GM, A_REF, J2 = 3.986004415e14, 6378136.3, 1.0826266835531513e-3
W6 = np.array([1, 2, 2, 1, 2, 1.0])[:, None]   # Frobenius weights of xx xy xz yy yz zz


def ncoef(N):
    return (N + 1) * (N + 2) // 2


def k_time(N):
    """Points per timing pass (one pass ~ 5-50 ms at every N)."""
    return int(np.clip(round(4e6 / ncoef(N)), 32, 4096))


# ------------------------------------------------------------------ model and points
def synthetic_model(N, seed=2009):
    """Fully normalised C, S (N+1, N+1): Kaula's rule 1e-5/n^2, EGM96 C20 C22 S22 C30, C00 = 1."""
    rng = np.random.default_rng(seed)
    C, S = np.zeros((N + 1, N + 1)), np.zeros((N + 1, N + 1))
    for n in range(2, N + 1):
        C[n, : n + 1] = rng.normal(0.0, 1e-5 / n**2, n + 1)
        S[n, 1 : n + 1] = rng.normal(0.0, 1e-5 / n**2, n)
    C[0, 0], C[2, 0], C[2, 1], S[2, 1] = 1.0, -4.84165371736e-4, 0.0, 0.0
    C[2, 2], S[2, 2], C[3, 0] = 2.43914352398e-6, -1.40016683654e-6, 9.57254173792e-7
    return C, S


def orbit_xyz(t, sma=6628e3, ecc=0.002, inc=96.5, raan0=-125.0):
    """250 km near-polar orbit with secular J2 drift, body-fixed (no Earth rotation) [m]."""
    i, n0, p = np.radians(inc), np.sqrt(GM / sma**3), sma * (1 - ecc**2)
    fac = 1.5 * J2 * (A_REF / p) ** 2
    raan = np.radians(raan0) - fac * n0 * np.cos(i) * t
    argp = 0.5 * fac * n0 * (5 * np.cos(i) ** 2 - 1) * t
    M = n0 * (1 + 0.5 * fac * np.sqrt(1 - ecc**2) * (3 * np.cos(i) ** 2 - 1)) * t
    E = M.copy()
    for _ in range(8):
        E -= (E - ecc * np.sin(E) - M) / (1 - ecc * np.cos(E))
    th = argp + 2 * np.arctan2(np.sqrt(1 + ecc) * np.sin(E / 2), np.sqrt(1 - ecc) * np.cos(E / 2))
    r = sma * (1 - ecc * np.cos(E))
    return r * np.stack([np.cos(raan) * np.cos(th) - np.sin(raan) * np.sin(th) * np.cos(i),
                         np.sin(raan) * np.cos(th) + np.cos(raan) * np.sin(th) * np.cos(i),
                         np.sin(th) * np.sin(i)])


def point_set(name, K, seed):
    """Body-fixed xyz (3, K) [m]. ORBIT1: K samples of a 23-day 250 km orbit (31 s step);
    UNIFORM: uniform on the sphere, 250-1000 km altitude; POLAR: same altitudes, colatitude
    log-uniform in [1e-8 deg, 1 deg]."""
    if name == "ORBIT1":
        return orbit_xyz(31.0 * np.linspace(0, 64439, K).round())
    rng = np.random.default_rng(1000 * (SETS.index(name) + 1) + seed)
    r = rng.uniform(6628e3, 7378e3, K)
    lam = rng.uniform(-np.pi, np.pi, K)
    if name == "UNIFORM":
        phi = np.arcsin(rng.uniform(-1.0, 1.0, K))
    else:
        phi = rng.choice([-1.0, 1.0], K) * (np.pi / 2 - np.radians(10.0 ** rng.uniform(-8.0, 0.0, K)))
    return r * np.stack([np.cos(phi) * np.cos(lam), np.cos(phi) * np.sin(lam), np.sin(phi)])


# ------------------------------------------------------------------ Rust binary
def build():
    src = [CRATE / "Cargo.toml", *(CRATE / "src").glob("*.rs")]
    if EXE.exists() and max(p.stat().st_mtime for p in src) <= EXE.stat().st_mtime:
        return
    cargo = shutil.which("cargo") or str(Path.home() / ".cargo" / "bin" / "cargo")
    print("[build] cargo build --release", flush=True)
    subprocess.run([cargo, "build", "--release", "--manifest-path", str(CRATE / "Cargo.toml")],
                   env=dict(os.environ, CARGO_TARGET_DIR=str(TARGET)), check=True)


def sph(*args):
    p = subprocess.run([str(EXE), *map(str, args)], capture_output=True, text=True)
    if p.returncode != 0:
        raise RuntimeError(f"sph_rs {args[0]} failed: {p.stderr.strip()}")
    return p


def write_coef(path, C, S):
    """u64 N, f64 GM, f64 a, C[(N+1)^2], S[(N+1)^2]."""
    with open(path, "wb") as fh:
        np.array([C.shape[0] - 1], "<u8").tofile(fh)
        np.array([GM, A_REF], "<f8").tofile(fh)
        np.ascontiguousarray(C, "<f8").tofile(fh)
        np.ascontiguousarray(S, "<f8").tofile(fh)
    return path


def write_xyz(path, xyz):
    """u64 K, x[K], y[K], z[K]."""
    with open(path, "wb") as fh:
        np.array([xyz.shape[1]], "<u8").tofile(fh)
        np.ascontiguousarray(xyz, "<f8").tofile(fh)
    return path


def read_out(path):
    """u64 rows, u64 K, f64 data[rows, K]."""
    with open(path, "rb") as fh:
        rows, k = (int(v) for v in np.fromfile(fh, "<u8", 2))
        return np.fromfile(fh, "<f8", rows * k).reshape(rows, k)


def evaluate(method, coef, pts, N, level, tmp):
    sph("eval", coef, pts, tmp / "out.bin", "--method", method, "--level", level, "--nmax", N)
    return read_out(tmp / "out.bin")


# ------------------------------------------------------------------ accuracy gate
def rel_err(e, nc):
    """Per point |e| / |nc| for the acceleration and (9 rows) the tensor, Frobenius norm."""
    with np.errstate(all="ignore"):
        ea = np.sqrt((e[:3] ** 2).sum(0) / (nc[:3] ** 2).sum(0))
        et = np.sqrt((W6 * e[3:9] ** 2).sum(0) / (W6 * nc[3:9] ** 2).sum(0)) if len(e) == 9 else None
    return ea, et


def run_gate(coef, n_list, k, tmp, threads):
    """One row per (N, set, method): max. relative error over the K points, pass at level 1 / level 2."""
    pts = write_xyz(tmp / "gate.bin", np.concatenate([point_set(s, k, seed=99) for s in SETS], axis=1))
    rows = []
    for N in n_list:
        t0 = time.time()
        sph("ref", coef, pts, tmp / "ref.bin", "--nmax", N, "--threads", threads)
        d = read_out(tmp / "ref.bin")
        hi, lo, nc = d[:9], d[9:18], d[18:27]        # reference = hi + lo, nc = reference - central term
        ba_a, ba_t = rel_err(d[27:36], nc)            # reference: Cartesian (Pines) - spherical form
        t_ref = time.time() - t0
        res = {}
        for m in METHODS:
            o1 = evaluate(m, coef, pts, N, 1, tmp)
            o2 = evaluate(m, coef, pts, N, 2, tmp)
            fin = bool(np.isfinite(o1).all() and np.isfinite(o2).all())
            res[m] = (o1, o2, rel_err((o2 - hi) - lo, nc), fin)
        for i, s in enumerate(SETS):
            sl = slice(i * k, (i + 1) * k)
            v = {}
            for m, (o1, o2, (ea, et), fin) in res.items():
                v[m] = dict(acc=float(ea[sl].max()) if fin else np.inf,
                            tensor=float(et[sl].max()) if fin else np.inf,
                            finite=fin, l1_same=bool(np.array_equal(o1[:, sl], o2[:3, sl])))
            best_a = min(x["acc"] for x in v.values())
            best_t = min(x["tensor"] for x in v.values())
            for m, x in v.items():
                ok1 = x["finite"] and x["l1_same"] and x["acc"] <= min(GATE_ABS, GATE_REL * best_a)
                ok2 = ok1 and x["tensor"] <= min(GATE_ABS, GATE_REL * best_t)
                rows.append(dict(N=N, set=s, method=m, K=k, **x, pass_L1=bool(ok1), pass_L2=bool(ok2),
                                 ref_check_acc=float(ba_a[sl].max()), ref_check_tensor=float(ba_t[sl].max())))
        last = rows[-len(METHODS) * len(SETS):]
        fail = sorted({r["method"] + ("" if r["finite"] else "(overflow)") for r in last
                       if not (r["pass_L1"] and r["pass_L2"])})
        print(f"[gate] N={N:5d}  ref {t_ref:5.1f} s (check {max(ba_a.max(), ba_t.max()):.0e})  "
              f"max rel. err. acc {max(r['acc'] for r in last if r['finite']):.1e} "
              f"grad {max(r['tensor'] for r in last if r['finite']):.1e}  "
              f"fail: {', '.join(fail) or '-'}", flush=True)
    return rows


# ------------------------------------------------------------------ timing
def run_timing(coef, n_list, gate, tmp, blocks, rounds, target_ms, seed, raw):
    """`rounds` passes over all cells, each pass in a new random order of cells and of methods; per pass one
    process per (cell, method) with `blocks` blocks.  A cell's blocks are thus spread over the whole run, so
    a slow change in machine load (other jobs sharing the L2 cache) cannot bias one cell or one method.
    All JSON records of sph_rs bench -> raw (jsonl.gz)."""
    finite = {(r["N"], r["method"]) for r in gate if r["finite"]}
    rng = np.random.default_rng(seed)
    cells0 = [(N, lev, s) for N in n_list for lev in LEVELS for s in SETS]
    t_all, run = time.time(), 0
    with gzip.open(raw, "wt") as fh:
        for rd in range(rounds):
            cells = [cells0[i] for i in rng.permutation(len(cells0))]       # no drift aligned with N
            for ci, (N, lev, s) in enumerate(cells):
                pts = write_xyz(tmp / "time.bin", point_set(s, k_time(N), seed=11))
                meths = [m for m in METHODS + (AA,) if (N, m.split("#")[0]) in finite]
                t0 = time.time()
                for pi in rng.permutation(len(meths)):
                    run += 1
                    p = sph("bench", coef, pts, "--level", lev, "--nmax", N, "--methods", meths[pi],
                            "--blocks", blocks, "--target-ms", target_ms, "--seed", seed + run)
                    for line in p.stdout.splitlines():
                        if line.startswith("{"):
                            fh.write(json.dumps(dict(json.loads(line), N=N, level=lev, set=s, round=rd)) + "\n")
                print(f"[time] round {rd + 1}/{rounds}  {ci + 1:2d}/{len(cells)}  N={N:5d} L{lev} {s:7s}  "
                      f"{len(meths)} processes  {time.time() - t0:5.1f} s  (total {(time.time() - t_all) / 60:5.1f} min)",
                      flush=True)


# ------------------------------------------------------------------ analysis
def boot_ratio(a, b, rng, nboot=4000):
    """median(b) / median(a) and its 95% bootstrap interval (blocks resampled independently)."""
    if a is None or b is None or min(a.size, b.size) < 3:
        return np.nan, np.nan, np.nan
    ra = np.median(a[rng.integers(0, a.size, (nboot, a.size))], axis=1)
    rb = np.median(b[rng.integers(0, b.size, (nboot, b.size))], axis=1)
    lo, hi = np.percentile(rb / ra, [2.5, 97.5])
    return float(np.median(b) / np.median(a)), float(lo), float(hi)


def analyze(recs, gate):
    """Per-block time per point normalised by the calibration loop -> timing and winner rows."""
    K = {(r["N"], r["level"], r["set"]): r["K"] for r in recs if r["t"] == "setup"}
    setup = {(r["N"], r["level"], r["set"], m): s for r in recs if r["t"] == "setup" for m, s in r["setup_s"].items()}
    meas = [r for r in recs if r["t"] == "m"]
    cal_ref = float(np.percentile([0.5 * (r["cal_ns"] + r["cal2_ns"]) for r in meas], 5))
    blk = {}
    for r in meas:
        key = (r["N"], r["level"], r["set"], r["method"])
        blk.setdefault(key, []).append((r["cpu"] / r["reps"] / K[key[:3]], 0.5 * (r["cal_ns"] + r["cal2_ns"])))
    t = {}
    for key, v in blk.items():
        v = np.array(v)
        t[key] = v[:, 0] * cal_ref / v[:, 1]              # all blocks: dropping blocks whose calibration
                                                           # drifted made the tails worse, not better
    passed = {(r["N"], r["set"], r["method"], lev) for r in gate for lev in LEVELS if r[f"pass_L{lev}"]}
    order = {m: i for i, m in enumerate(METHODS + (AA,))}
    rng = np.random.default_rng(2009)
    trows, wrows = [], []
    for key in sorted(t, key=lambda k: (k[1], SETS.index(k[2]), k[0], order[k[3]])):
        N, lev, s, m = key
        med = float(np.median(t[key]))
        lo, hi = np.percentile(np.median(t[key][rng.integers(0, t[key].size, (4000, t[key].size))], axis=1), [2.5, 97.5])
        trows.append(dict(N=N, level=lev, set=s, method=m, K=K[key[:3]], us_per_pt=med * 1e6, ci_lo=lo * 1e6,
                          ci_hi=hi * 1e6, blocks_used=t[key].size,
                          setup_s=setup.get(key, np.nan), gate_pass=(N, s, m.split("#")[0], lev) in passed))
    for N, lev, s in sorted(K, key=lambda k: (k[1], SETS.index(k[2]), k[0])):
        tt = {m: t[(N, lev, s, m)] for m in METHODS + (AA,) if (N, lev, s, m) in t}
        aa = boot_ratio(tt.get("LEGENDRE_a"), tt.get(AA), rng)
        delta = max(DELTA_MIN, 2 * abs(aa[0] - 1)) if np.isfinite(aa[0]) else DELTA_MIN
        cand = sorted((m for m in METHODS if m in tt and (N, s, m, lev) in passed), key=lambda m: np.median(tt[m]))
        if not cand:
            continue
        rat = {m: boot_ratio(tt[cand[0]], tt[m], rng) for m in cand[1:]}
        tied = [cand[0]] + [m for m in cand[1:] if not rat[m][1] > 1 + delta]
        ru = cand[1] if len(cand) > 1 else None
        wrows.append(dict(N=N, level=lev, set=s, fastest=cand[0], verdict="WIN" if len(tied) == 1 else "TIE",
                          tied="+".join(tied), runner_up=ru or "", ratio=rat[ru][0] if ru else np.nan,
                          ratio_lo=rat[ru][1] if ru else np.nan, ratio_hi=rat[ru][2] if ru else np.nan,
                          delta=delta, aa_ratio=aa[0], aa_lo=aa[1], aa_hi=aa[2],
                          failed_gate="+".join(m for m in METHODS if m not in cand)))
    return trows, wrows, cal_ref


# ------------------------------------------------------------------ tables
def write_csv(path, rows):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def read_csv(path):
    rows = list(csv.DictReader(open(path)))
    for r in rows:
        for k, v in r.items():
            if v in ("True", "False"):
                r[k] = v == "True"
            elif k in ("N", "level", "K", "blocks_used"):
                r[k] = int(v)
            else:
                try:
                    r[k] = float(v)
                except ValueError:
                    pass
    return rows


def summary_md(gate, trows, wrows, meta):
    Ns = sorted({r["N"] for r in gate})
    head = "| Method | " + " | ".join(f"N={n}" for n in Ns) + " |\n|---|" + "---:|" * len(Ns)
    out = ["# Results", "",
           f"Generated by `benchmark.py` on {meta.get('machine', {}).get('cpu', '?')}, "
           f"{meta.get('machine', {}).get('rustc', '')}. Synthetic Earth-like model (Kaula), degree N.", "",
           "## Accuracy", "",
           "Max. relative error = largest value, over all test points, of "
           "|computed - reference| / |reference - point-mass term| "
           "(vector norm for the acceleration, Frobenius norm for the gravity gradient). "
           "Dividing by the field minus its point-mass part GM/r^2 measures the error in the part the harmonics add. "
           f"Test points: {len(SETS)} sets x {gate[0]['K']} points ({', '.join(SETS)}). "
           "Reference = double-double synthesis (sph_rs ref). "
           f"Pass: error <= {GATE_ABS:.0e} and <= {GATE_REL:.0f} x the best method; "
           "overflow (Inf/NaN) fails the method at that N.", ""]
    for key, ttl in (("acc", "Acceleration"), ("tensor", "Gravity gradient")):
        out += [f"**{ttl}**", "", head]
        for m in METHODS:
            cells = []
            for N in Ns:
                rr = [r for r in gate if r["N"] == N and r["method"] == m]
                cells.append("overflow" if not all(r["finite"] for r in rr) else f"{max(r[key] for r in rr):.1e}")
            out.append(f"| {m} | " + " | ".join(cells) + " |")
        out.append("")
    chk = max(max(r["ref_check_acc"], r["ref_check_tensor"]) for r in gate)
    out += [f"Reference check (two independent double-double formulations): max difference {chk:.0e}.", ""]
    win = {(r["N"], r["level"], r["set"]): r for r in wrows}
    out += ["## CPU time", "",
            "CPU time per point in microseconds (ORBIT1 points; median of the blocks, normalised by the calibration loop). "
            "**bold** = fastest; tied methods (difference not resolvable, see winners.csv) marked with *. "
            "x = failed the accuracy check, - = overflow (not timed).", ""]
    for lev in LEVELS:
        out += [f"**{'Acceleration' if lev == 1 else 'Acceleration + gravity gradient'}**", "", head]
        for m in METHODS:
            cells = []
            for N in Ns:
                r = next((r for r in trows if (r["N"], r["level"], r["set"], r["method"]) == (N, lev, "ORBIT1", m)), None)
                w = win.get((N, lev, "ORBIT1"))
                if r is None:                             # overflow, or the timing phase has not run yet
                    ovf = not all(g["finite"] for g in gate if g["N"] == N and g["method"] == m)
                    cells.append("-" if ovf else "not timed")
                    continue
                u = r["us_per_pt"]
                v = f"{u:.0f}" if u >= 100 else f"{u:#.3g}"
                if not r["gate_pass"]:
                    v += " x"
                elif w and w["fastest"] == m:
                    v = f"**{v}**"
                elif w and m in w["tied"].split("+"):
                    v += "*"
                cells.append(v)
            out.append(f"| {m} | " + " | ".join(cells) + " |")
        out.append("")
    out += ["## Fastest method per N (all point sets)", "", "| Output | Points | " + " | ".join(f"N={n}" for n in Ns) + " |",
            "|---|---|" + "---|" * len(Ns)]
    for lev in LEVELS:
        for s in SETS:
            cells = []
            for N in Ns:
                w = win.get((N, lev, s))
                cells.append("" if w is None else w["fastest"] if w["verdict"] == "WIN" else "tie: " + w["tied"])
            out.append(f"| {'Acc.' if lev == 1 else 'Acc. + grad.'} | {s} | " + " | ".join(cells) + " |")
    aa = [abs(r["aa_ratio"] - 1) for r in wrows if np.isfinite(r["aa_ratio"])]
    out += ["", f"A/A probe (LEGENDRE_a timed against an identical second copy of itself): median |ratio - 1| = {np.median(aa):.1%}, "
            f"max {np.max(aa):.1%}; a win needs the 95% interval of the time ratio above 1 + max({DELTA_MIN:.0%}, "
            "2 |A/A - 1|).", ""]
    rat = {}
    for r in trows:
        ref = next((q for q in trows if (q["N"], q["level"], q["set"], q["method"]) ==
                    (r["N"], r["level"], r["set"], "LEGENDRE_a")), None)
        if r["method"] in METHODS[1:] and r["gate_pass"] and ref:
            rat.setdefault((r["N"], r["level"], r["method"]), []).append(r["us_per_pt"] / ref["us_per_pt"])
    spr = {k: (max(v) - min(v)) / np.median(v) for k, v in rat.items() if len(v) == len(SETS)}
    if spr:
        kmax = max(spr, key=spr.get)
        out += [f"Point-set check (the cost of a point does not depend on where it is, so the {len(SETS)} point sets are "
                f"repeat measurements): spread (max - min) / median of the time ratio to LEGENDRE_a over the sets = "
                f"{np.median(list(spr.values())):.1%} (median), {spr[kmax]:.1%} (max, {kmax[2]} level {kmax[1]} N={kmax[0]}).", ""]
    return "\n".join(out)


# ------------------------------------------------------------------ driver
def machine_info():
    def run(*cmd):
        try:
            return subprocess.run(cmd, capture_output=True, text=True).stdout.strip()
        except OSError:
            return ""
    rustc = shutil.which("rustc") or str(Path.home() / ".cargo" / "bin" / "rustc")
    rv = dict(line.split(": ", 1) for line in run(rustc, "-vV").splitlines() if ": " in line)
    return dict(cpu=run("sysctl", "-n", "machdep.cpu.brand_string") or platform.processor(),
                perf_cores=run("sysctl", "-n", "hw.perflevel0.physicalcpu"), os=platform.mac_ver()[0] or platform.platform(),
                rustc=f"rustc {rv.get('release', '?')} ({rv.get('host', '?')}), opt-level 3, LTO",
                loadavg=[round(x, 2) for x in os.getloadavg()], date=time.strftime("%Y-%m-%d %H:%M"))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--phase", default="all", choices=("all", "gate", "time", "analyze"))
    ap.add_argument("--out", default=None, help="output folder (default results/)")
    ap.add_argument("--N", type=int, nargs="+", default=list(N_LIST))
    ap.add_argument("--blocks", type=int, default=4, help="timed blocks per process")
    ap.add_argument("--rounds", type=int, default=3, help="passes over all cells (one process per method each)")
    ap.add_argument("--target-ms", type=float, default=100.0)
    ap.add_argument("--k-gate", type=int, default=48, help="gate points per set")
    ap.add_argument("--threads", type=int, default=4, help="threads of the reference")
    ap.add_argument("--seed", type=int, default=2009)
    ap.add_argument("--quick", action="store_true", help="smoke test")
    a = ap.parse_args()
    if a.quick:
        a.N, a.blocks, a.rounds, a.target_ms, a.k_gate = [10, 50, 360], 3, 2, 20.0, 8
    out = Path(a.out or HERE / ("results_quick" if a.quick else "results"))
    out.mkdir(parents=True, exist_ok=True)
    meta_path, raw = out / "meta.json", out / "raw_timing.jsonl.gz"
    meta = json.loads(meta_path.read_text()) if meta_path.exists() else {}
    if a.phase != "analyze":                    # analyze needs only the files saved in out/
        build()
        with tempfile.TemporaryDirectory(prefix="sph_") as td:
            tmp = Path(td)
            coef = write_coef(tmp / "coef.bin", *synthetic_model(max(a.N)))
            if a.phase in ("all", "gate"):
                t0 = time.time()
                gate = run_gate(coef, a.N, a.k_gate, tmp, a.threads)
                write_csv(out / "accuracy.csv", gate)
                meta["machine"] = machine_info()
                meta["gate"] = dict(N=a.N, K_per_set=a.k_gate, ref_threads=a.threads, minutes=(time.time() - t0) / 60)
                meta_path.write_text(json.dumps(meta, indent=1))
            gate = read_csv(out / "accuracy.csv")
            if a.phase in ("all", "time"):
                t0, m0 = time.time(), machine_info()
                run_timing(coef, a.N, gate, tmp, a.blocks, a.rounds, a.target_ms, a.seed, raw)
                meta["timing"] = dict(N=a.N, blocks=a.blocks, rounds=a.rounds, target_ms=a.target_ms, seed=a.seed,
                                      K={N: k_time(N) for N in a.N}, minutes=(time.time() - t0) / 60,
                                      loadavg_start=m0["loadavg"], loadavg_end=machine_info()["loadavg"])
                meta_path.write_text(json.dumps(meta, indent=1))
    gate = read_csv(out / "accuracy.csv")
    with gzip.open(raw, "rt") as fh:
        recs = [json.loads(line) for line in fh]
    trows, wrows, cal_ref = analyze(recs, gate)
    write_csv(out / "timing.csv", trows)
    write_csv(out / "winners.csv", wrows)
    meta["analysis"] = dict(cal_ref_ns=cal_ref, delta_min=DELTA_MIN, gate_abs=GATE_ABS, gate_rel=GATE_REL)
    meta_path.write_text(json.dumps(meta, indent=1))
    (out / "summary.md").write_text(summary_md(gate, trows, wrows, meta))
    from plots import figures
    figures(out)
    print((out / "summary.md").read_text())


if __name__ == "__main__":
    main()
