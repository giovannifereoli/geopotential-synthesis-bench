"""
run_all.py -- every experiment of the Fantino & Casotto (2009) emulation, in order:

    sph_efficiency.py   Table 7 + Fig. 4 + cost vs degree
    sph_precision.py    Figs. 5-10 + precision summary + conditioning analysis
    sph_accuracy.py     Eqs. 31-33 + Fig. 3 (re-uses the precision cache)

    python run_all.py                 # default setup (~10-15 min on one core)
    python run_all.py --quick         # smoke test (~3 min), results in ./quick_results
    python run_all.py --out DIR
    python run_all.py --skip efficiency
    python run_all.py --backend rust  # compiled kernels (sph_rs/), results in ./results_rust
"""

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")  # before numpy is imported: single thread, as in the paper

import argparse
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--quick", action="store_true", help="smoke-test setup of every script")
    ap.add_argument("--out", type=str, default=None,
                    help="output folder (default ./results, or ./quick_results with --quick)")
    ap.add_argument("--N", type=int, default=360, help="maximum degree (paper: 360)")
    ap.add_argument("--gfc", type=str, default=None, help="ICGEM .gfc model (default: synthetic)")
    ap.add_argument("--backend", choices=("numpy", "rust"), default="numpy",
                    help="numpy (vectorised) or rust (compiled scalar kernels; needs cargo)")
    ap.add_argument("--skip", nargs="*", default=(), choices=("efficiency", "precision", "accuracy"))
    a = ap.parse_args(argv)
    sfx = "_rust" if a.backend == "rust" else ""
    out = a.out or str(HERE / (("quick_results" if a.quick else "results") + sfx))
    common = ["--out", out, "--N", str(a.N), "--backend", a.backend] + (["--gfc", a.gfc] if a.gfc else [])
    common += ["--quick"] if a.quick else []

    import sph_accuracy
    import sph_efficiency
    import sph_precision

    steps = (("efficiency", sph_efficiency), ("precision", sph_precision),
             ("accuracy", sph_accuracy))
    t_all = time.perf_counter()
    for name, mod in steps:
        if name in a.skip:
            continue
        t0 = time.perf_counter()
        print(f"\n{'=' * 70}\n[run_all] {name}\n{'=' * 70}")
        mod.main(list(common))
        print(f"[run_all] {name} done in {(time.perf_counter() - t0) / 60:.1f} min")
    print(f"\n[run_all] all done in {(time.perf_counter() - t_all) / 60:.1f} min -> {out}")


if __name__ == "__main__":
    main()
