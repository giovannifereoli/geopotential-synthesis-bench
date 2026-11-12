"""
sph_rust.py -- compiled backend for the five synthesis programmes (crate sph_rs/).

The kernels of sph_core.py are re-written in Rust as scalar, point-by-point loops (the way
the paper's Fortran90 is organised), so that the timings measure the algorithms and not the
NumPy vectorisation strategy. The binary is a stand-alone native executable called through
subprocess with flat binary files: it runs natively (arm64 on Apple silicon) even when this
Python is an x86-64 build under Rosetta, and the synthesis loop is timed inside the binary
(process CPU time, best of `reps`), so that file I/O and process start-up are excluded.

    import sph_rust as sr
    res = sr.run("PINES", model, lam, phi, r, reps=3)   # res["loc"], res["cart"], res["cpu"]
    res = sr.run_grid(model, lam1d, phi1d, r)           # CLENSHAW per parallel
    loc = sr.synthesize("METRIS", model, lam, phi, r)   # drop-in for sph_core.synthesize

The crate is built on first use with `cargo build --release` (no external crates). The build
products go to ~/.cache/sph_rs/target (override with SPH_RS_TARGET) to keep them out of a
cloud-synced project folder.
"""

import json
import os
import platform
import shutil
import subprocess
import tempfile
from pathlib import Path

import numpy as np

import sph_core as sc

CRATE = Path(__file__).resolve().parent / "sph_rs"
TARGET = Path(os.environ.get("SPH_RS_TARGET", Path.home() / ".cache" / "sph_rs" / "target"))
BINARY = TARGET / "release" / ("sph_rs.exe" if os.name == "nt" else "sph_rs")
_BUILT = False
_TMP = None


def _cargo():
    exe = shutil.which("cargo") or str(Path.home() / ".cargo" / "bin" / "cargo")
    if not Path(exe).exists():
        raise RuntimeError("cargo not found: install Rust (https://rustup.rs) for --backend rust")
    return exe


def ensure_built(verbose=True):
    """Build the crate if the binary is missing or older than its sources."""
    global _BUILT
    if _BUILT:
        return BINARY
    src = [CRATE / "Cargo.toml", CRATE / "src" / "main.rs"]
    stale = (not BINARY.exists()) or any(p.stat().st_mtime > BINARY.stat().st_mtime for p in src)
    if stale:
        if verbose:
            print(f"[rust] cargo build --release ({CRATE.name}) ...")
        env = dict(os.environ, CARGO_TARGET_DIR=str(TARGET))
        out = subprocess.run([_cargo(), "build", "--release", "--manifest-path",
                              str(CRATE / "Cargo.toml")], env=env, capture_output=True, text=True)
        if out.returncode != 0:
            raise RuntimeError("cargo build failed:\n" + out.stderr)
    _BUILT = True
    return BINARY


def _tmpdir():
    global _TMP
    if _TMP is None:
        _TMP = tempfile.TemporaryDirectory(prefix="sph_rs_")
    return Path(_TMP.name)


def _f8(x):
    return np.ascontiguousarray(np.atleast_1d(np.asarray(x, dtype="<f8")))


def _model_file(model):
    """Coefficient file of `model` (written once per model object)."""
    path = model._cache.get("rust_file")
    if path is None or not Path(path).exists():
        path = _tmpdir() / f"coef_{id(model):x}.bin"
        with open(path, "wb") as fh:
            np.array([model.N], "<u8").tofile(fh)
            np.array([model.gm, model.a], "<f8").tofile(fh)
            _f8(model.C).tofile(fh)
            _f8(model.S).tofile(fh)
        model._cache["rust_file"] = path
    return path


def _call(method, model, pts_path, K, reps, nmax, grid):
    exe = ensure_built()
    out_path = _tmpdir() / "out.bin"
    cmd = [str(exe), method, str(_model_file(model)), str(pts_path), str(out_path),
           "--reps", str(int(reps))]
    if nmax is not None:
        cmd += ["--nmax", str(int(nmax))]
    if grid:
        cmd += ["--grid"]
    p = subprocess.run(cmd, capture_output=True, text=True)
    if p.returncode != 0:
        raise RuntimeError(f"sph_rs failed ({method}): {p.stderr.strip()}")
    info = json.loads(p.stdout.strip().splitlines()[-1])
    with open(out_path, "rb") as fh:
        rows, k = (int(v) for v in np.fromfile(fh, "<u8", 2))
        data = np.fromfile(fh, "<f8", rows * k).reshape(rows, k).astype(np.float64)
    assert k == K, (k, K)
    info["loc"] = data[:10]
    info["cart"] = data[10:] if rows == 20 else None
    return info


def run(method, model, lam, phi, r, reps=1, nmax=None):
    """
    Synthesis at scattered points with the compiled kernel of `method`.

    Returns a dict: loc (10, K) local-frame quantities (QUANTS order), cart (10, K) or None
    (native Cartesian output of PINES / METRIS), cpu / wall [s] of the synthesis loop (best of
    `reps`), setup [s] (reading the files, recursion tables, METRIS coefficient transforms).
    `nmax` truncates the model to degree nmax inside the binary.
    """
    lam, phi, r = (_f8(x) for x in np.broadcast_arrays(lam, phi, r))
    pts = _tmpdir() / "points.bin"
    with open(pts, "wb") as fh:
        np.array([lam.size], "<u8").tofile(fh)
        for x in (lam, phi, r):
            x.tofile(fh)
    return _call(method, model, pts, lam.size, reps, nmax, grid=False)


def run_grid(model, lam1d, phi1d, r, reps=1, nmax=None):
    """CLENSHAW on a grid at constant r (recursions once per parallel); latitude-major."""
    lam1d, phi1d = _f8(lam1d), _f8(phi1d)
    pts = _tmpdir() / "grid.bin"
    with open(pts, "wb") as fh:
        np.array([lam1d.size, phi1d.size], "<u8").tofile(fh)
        np.array([float(np.ravel(r)[0])], "<f8").tofile(fh)
        lam1d.tofile(fh)
        phi1d.tofile(fh)
    return _call("CLENSHAW", model, pts, lam1d.size * phi1d.size, reps, nmax, grid=True)


def synthesize(method, model, lam, phi, r, cartesian=False):
    """Drop-in for sph_core.synthesize(..., dtype=np.float64) on the compiled backend."""
    res = run(method, model, lam, phi, r)
    if not cartesian:
        return res["loc"]
    if res["cart"] is not None:
        return res["loc"], res["cart"]
    lam_, phi_ = (_f8(x) for x in np.broadcast_arrays(lam, phi, r)[:2])
    return res["loc"], sc.local_to_cart(res["loc"], lam_, phi_)


def clenshaw_grid(model, lam1d, phi1d, r):
    """Drop-in for sph_core.clenshaw_grid(..., dtype=np.float64)."""
    return run_grid(model, lam1d, phi1d, r)["loc"]


def toolchain_string():
    """e.g. 'rustc 1.87.0, aarch64-apple-darwin, opt-level 3 + LTO, no fast-math'."""
    try:
        exe = shutil.which("rustc") or str(Path.home() / ".cargo" / "bin" / "rustc")
        v = subprocess.run([exe, "-vV"], capture_output=True, text=True).stdout
        d = dict(line.split(": ", 1) for line in v.splitlines() if ": " in line)
        return f"rustc {d.get('release', '?')}, {d.get('host', platform.machine())}, opt-level 3 + LTO"
    except Exception:
        return "rustc (unknown version)"


def _selftest(N=30, K=40, seed=1):
    rng = np.random.default_rng(seed)
    model = sc.GravityModel.synthetic(N, seed=7)
    lam = rng.uniform(-np.pi, np.pi, K)
    phi = np.radians(rng.uniform(-85, 85, K))
    r = sc.A_EARTH + rng.uniform(200e3, 500e3, K)
    ref = sc.clenshaw(model, lam, phi, r, dtype=np.longdouble)
    print(f"[{toolchain_string()}]")
    print(f"self-test N={N}, K={K}: min coinciding digits of the Rust kernels vs CLENSHAW*16")
    print("            " + " ".join(f"{q:>5s}" for q in sc.QUANTS))
    for meth in sc.METHODS:
        d = sc.digits(run(meth, model, lam, phi, r)["loc"], ref).min(axis=1)
        print(f"{meth:11s} " + " ".join(f"{v:5.1f}" for v in d))
    lon, lat = np.array([-120.0, 10.0, 77.0]), np.array([-60.0, 5.0, 33.0, 80.0])
    gl, gp, gr = sc.grid_points(300e3, lon, lat)
    a_ = run_grid(model, np.radians(lon), np.radians(lat), gr[0])["loc"]
    b_ = sc.clenshaw(model, gl, gp, gr, dtype=np.longdouble)
    print(f"CLENSHAW grid (Rust) vs CLENSHAW*16 min digits: {sc.digits(a_, b_).min():.1f}")
    cp_ = run("PINES", model, lam, phi, r)["cart"]
    cm_ = run("METRIS", model, lam, phi, r)["cart"]
    print(f"PINES vs METRIS (Rust, Cartesian) min digits: {sc.digits(cp_, cm_).min(axis=1).round(1)}")


if __name__ == "__main__":
    _selftest()
