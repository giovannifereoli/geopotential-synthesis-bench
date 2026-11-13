//! sph_rs -- scalar Rust kernels for spherical-harmonic gravity synthesis: acceleration (level 1)
//! and acceleration + gravity-gradient tensor (level 2) at scattered body-fixed points.
//! No external crates, IEEE double precision, no fast-math, single thread per kernel.
//!
//!   sph_rs bench <coef.bin> <xyz.bin> [options]          CPU-time benchmark      (src/bench.rs)
//!   sph_rs eval  <coef.bin> <xyz.bin> <out.bin> --method M   evaluate one method (src/bench.rs)
//!   sph_rs ref   <coef.bin> <xyz.bin> <out.bin> [--nmax N]   double-double reference (src/reference.rs)
//!
//! coef.bin : u64 N, f64 gm, f64 a, C[(N+1)^2], S[(N+1)^2]   (row-major [n, m], f64 LE)
//! xyz.bin  : u64 K, x[K], y[K], z[K]                         (body-fixed Cartesian, m)

mod bench;
mod reference;

use std::env;
use std::fs;

// ------------------------------------------------------------------ CPU time
#[repr(C)]
struct Timespec {
    tv_sec: i64,
    tv_nsec: i64,
}
extern "C" {
    fn clock_gettime(clk: i32, ts: *mut Timespec) -> i32;
}

// ------------------------------------------------------------------ binary I/O
fn rd_u64(b: &[u8], off: &mut usize) -> usize {
    let v = u64::from_le_bytes(b[*off..*off + 8].try_into().unwrap());
    *off += 8;
    v as usize
}
fn rd_f64(b: &[u8], off: &mut usize) -> f64 {
    let v = f64::from_le_bytes(b[*off..*off + 8].try_into().unwrap());
    *off += 8;
    v
}
fn rd_vec(b: &[u8], off: &mut usize, n: usize) -> Vec<f64> {
    (0..n).map(|_| rd_f64(b, off)).collect()
}

// ------------------------------------------------------------------ tables
struct Tables {
    l: usize,
    f: Vec<f64>,
    g: Vec<f64>, // [n * l + m]
    h: Vec<f64>,
    k: Vec<f64>,
    j: Vec<f64>,
    r11: Vec<f64>,
    r12: Vec<f64>,
    r22: Vec<f64>,
    alpha: Vec<f64>,
    gamma: Vec<f64>,
    gt: Vec<f64>, // [m * l + n]  (column-major copies for the Clenshaw recursions)
    ht: Vec<f64>,
}

fn tables(nmax: usize) -> Tables {
    let l = nmax + 5;
    let z = || vec![0.0f64; l * l];
    let mut t = Tables {
        l,
        f: vec![0.0; l],
        g: z(), h: z(), k: z(), j: z(), r11: z(), r12: z(), r22: z(), alpha: z(), gamma: z(),
        gt: z(), ht: z(),
    };
    for ni in 1..l {
        let n = ni as f64;
        t.f[ni] = ((2.0 * n + 1.0) * (if ni == 1 { 2.0 } else { 1.0 }) / (2.0 * n)).sqrt();
    }
    for ni in 0..l {
        for mi in 0..=ni.min(l - 1) {
            let (n, m) = (ni as f64, mi as f64);
            let eps = if mi == 0 { 0.5 } else { 1.0 };
            let i = ni * l + mi;
            if ni > mi {
                t.g[i] = ((2.0 * n + 1.0) * (2.0 * n - 1.0) / ((n + m) * (n - m))).sqrt();
            }
            if ni > mi + 1 {
                t.h[i] = ((2.0 * n + 1.0) * (n - m - 1.0) * (n + m - 1.0)
                    / ((2.0 * n - 3.0) * (n + m) * (n - m))).sqrt();
            }
            t.k[i] = (eps * (n - m) * (n + m + 1.0)).sqrt();
            if ni >= mi + 2 {
                t.j[i] = (eps * (n - m) * (n - m - 1.0) * (n + m + 1.0) * (n + m + 2.0)).sqrt();
            }
            t.r11[i] = (eps * (2.0 * n + 1.0) * (n + m + 1.0) * (n + m + 2.0) / (2.0 * n + 3.0)).sqrt();
            t.r12[i] = (eps * (2.0 * n + 1.0) * (n - m) * (n + m + 1.0) * (n + m + 2.0)
                * (n + m + 3.0) / (2.0 * n + 3.0)).sqrt();
            t.r22[i] = (eps * (2.0 * n + 1.0) * (n + m + 1.0) * (n + m + 2.0) * (n + m + 3.0)
                * (n + m + 4.0) / (2.0 * n + 5.0)).sqrt();
            t.alpha[i] = ((n - m + 1.0) * (n + m + 1.0) * (2.0 * n + 1.0) / (2.0 * n + 3.0)).sqrt();
            if mi >= 1 {
                t.gamma[i] = ((if mi == 1 { 2.0 } else { 1.0 }) * (2.0 * n + 1.0) * (n - m + 2.0)
                    * (n - m + 1.0) / (2.0 * n + 3.0)).sqrt();
            }
            t.gt[mi * l + ni] = t.g[i];
            t.ht[mi * l + ni] = t.h[i];
        }
    }
    t
}

// ------------------------------------------------------------------ model
struct Model {
    n: usize,
    n1: usize,
    gm: f64,
    a: f64,
    c: Vec<f64>,  // [n * n1 + m]
    s: Vec<f64>,
    ct: Vec<f64>, // [m * n1 + n]
    st: Vec<f64>,
    t: Tables,
}

fn load_model(path: &str, nmax: Option<usize>) -> Model {
    let b = fs::read(path).expect("cannot read coefficient file");
    let mut off = 0;
    let nf = rd_u64(&b, &mut off);
    let gm = rd_f64(&b, &mut off);
    let a = rd_f64(&b, &mut off);
    let cf = rd_vec(&b, &mut off, (nf + 1) * (nf + 1));
    let sf = rd_vec(&b, &mut off, (nf + 1) * (nf + 1));
    let n = nmax.map_or(nf, |x| x.min(nf));
    let n1 = n + 1;
    let mut c = vec![0.0; n1 * n1];
    let mut s = vec![0.0; n1 * n1];
    let mut ct = vec![0.0; n1 * n1];
    let mut st = vec![0.0; n1 * n1];
    for i in 0..n1 {
        for m in 0..=i {
            c[i * n1 + m] = cf[i * (nf + 1) + m];
            s[i * n1 + m] = if m == 0 { 0.0 } else { sf[i * (nf + 1) + m] };
            ct[m * n1 + i] = c[i * n1 + m];
            st[m * n1 + i] = s[i * n1 + m];
        }
    }
    Model { n, n1, gm, a, c, s, ct, st, t: tables(n) }
}

// ------------------------------------------------------------------ METRIS operators
fn metris_op(er: &[f64], ei: &[f64], kind: char, t: &Tables, l: usize, a: f64) -> (Vec<f64>, Vec<f64>) {
    let lt = t.l;
    let mut or = vec![0.0f64; l * l];
    let mut oi = vec![0.0f64; l * l];
    for n in 0..l - 1 {
        for m in 0..l {
            let (sr, si) = (er[n * l + m], ei[n * l + m]);
            let o = (n + 1) * l;
            if kind == 'z' {
                let al = t.alpha[n * lt + m];
                or[o + m] -= al * sr;
                oi[o + m] -= al * si;
                continue;
            }
            let isx = kind == 'x';
            let sgn = if isx { -1.0 } else { 1.0 };
            // half * src, half = 1/2 (x) or i/2 (y)
            let (hr, hi) = if isx { (0.5 * sr, 0.5 * si) } else { (-0.5 * si, 0.5 * sr) };
            let be = t.r11[n * lt + m];
            if m + 1 < l {
                or[o + m + 1] += sgn * be * hr;
                oi[o + m + 1] += sgn * be * hi;
            }
            if m >= 1 {
                let ga = t.gamma[n * lt + m];
                or[o + m - 1] += ga * hr;
                oi[o + m - 1] += ga * hi;
            }
            if m == 0 {
                // half * conj(src)
                let (cr, ci) = if isx { (0.5 * sr, -0.5 * si) } else { (0.5 * si, 0.5 * sr) };
                or[o + 1] += sgn * be * cr;
                oi[o + 1] += sgn * be * ci;
            }
        }
    }
    for v in or.iter_mut() {
        *v /= a;
    }
    for v in oi.iter_mut() {
        *v /= a;
    }
    (or, oi)
}


fn main() {
    let args: Vec<String> = env::args().collect();
    match args.get(1).map(|s| s.as_str()) {
        Some("bench") | Some("eval") => bench::main(&args),
        Some("ref") => reference::main(&args),
        _ => {
            eprintln!("usage: sph_rs bench|eval|ref ... (see src/main.rs)");
            std::process::exit(2);
        }
    }
}
