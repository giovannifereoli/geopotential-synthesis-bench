# Hardware counters (instructions retired, cycles) per point per method, via /usr/bin/time -l.
# Two point counts per run; the difference removes start-up and model loading.
import sys, os, re, subprocess, json
from pathlib import Path
sys.path.insert(0, os.getcwd())
import benchmark as B
out = Path(sys.argv[1]); N = int(sys.argv[2]); out.mkdir(exist_ok=True)
coef = B.write_coef(out / "coef.bin", *B.synthetic_model(N))
Ks = (int(sys.argv[3]), int(sys.argv[4]))
pts = {k: B.write_xyz(out / f"x{k}.bin", B.point_set("ORBIT1", k, seed=11)) for k in Ks}
def run(exe, m, lev, k):
    p = subprocess.run(["/usr/bin/time", "-l", str(exe), "eval", str(coef), str(pts[k]), str(out / "o.bin"),
                        "--method", m, "--level", str(lev), "--nmax", str(N)], capture_output=True, text=True)
    g = lambda key: int(re.search(r"(\d+)\s+" + key, p.stderr).group(1))
    return g("instructions retired"), g("cycles elapsed")
res = []
for build in ("scalar", "simd"):
    exe = B.exe_of(build)
    for lev in (1, 2):
        for m in B.METHODS:
            a = [run(exe, m, lev, k) for k in Ks]
            b = [run(exe, m, lev, k) for k in Ks]  # second pass: keep the lower cycle count
            dk = Ks[1] - Ks[0]
            ins = (min(a[1][0], b[1][0]) - min(a[0][0], b[0][0])) / dk
            cyc = (min(a[1][1], b[1][1]) - min(a[0][1], b[0][1])) / dk
            res.append(dict(build=build, level=lev, method=m, ins=ins, cyc=cyc))
            print(build, lev, m, f"{ins:.4g} ins/pt  {cyc:.4g} cyc/pt  IPC {ins/cyc:.2f}", flush=True)
json.dump(res, open(out / f"hw_N{N}.json", "w"), indent=1)
