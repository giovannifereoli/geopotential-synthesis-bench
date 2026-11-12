//! sph_rs -- compiled (Rust) kernels of the harmonic-synthesis programmes of
//! Fantino & Casotto (2009), J. Geodesy 83:595-619. Same mathematics as sph_core.py, written
//! point by point with scalar loops (as the paper's Fortran90), no external crates, IEEE
//! double precision, no fast-math, single thread.
//!
//!   sph_rs <METHOD> <coef.bin> <points.bin> <out.bin> [--reps R] [--nmax n] [--grid]
//!
//! METHOD: LEGENDRE_a | LEGENDRE_b | PINES | GOTTLIEB | METRIS | CLENSHAW
//!         (--grid: CLENSHAW per parallel)
//! coef.bin   : u64 N, f64 gm, f64 a, C[(N+1)^2], S[(N+1)^2]   (row-major [n, m], f64 LE)
//! points.bin : u64 K, lam[K], phi[K], r[K]                    (rad, rad, m)
//!   (--grid) : u64 nlon, u64 nlat, f64 r, lon[nlon], lat[nlat]
//! out.bin    : u64 rows, u64 K, data[rows * K] row-major; rows = 10 (local frame, QUANTS
//!              order) or 20 (local + Cartesian, CART order) for PINES, GOTTLIEB, METRIS
//! stdout     : one JSON line {"cpu":..,"wall":..,"setup":..,"K":..,"N":..,"reps":..}
//!              cpu / wall = best of R repetitions of the synthesis loop only.

use std::env;
use std::fs;
use std::io::Write;
use std::time::Instant;

// ------------------------------------------------------------------ CPU time
#[repr(C)]
struct Timespec {
    tv_sec: i64,
    tv_nsec: i64,
}
extern "C" {
    fn clock_gettime(clk: i32, ts: *mut Timespec) -> i32;
}
#[cfg(target_os = "macos")]
const CLOCK_PROCESS_CPUTIME_ID: i32 = 12;
#[cfg(not(target_os = "macos"))]
const CLOCK_PROCESS_CPUTIME_ID: i32 = 2;

fn cpu_time() -> f64 {
    let mut ts = Timespec { tv_sec: 0, tv_nsec: 0 };
    unsafe {
        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &mut ts);
    }
    ts.tv_sec as f64 + 1e-9 * ts.tv_nsec as f64
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

// ------------------------------------------------------------------ frames
#[inline]
fn spherical_to_local(d: &[f64; 10], r: f64, phi: f64) -> [f64; 10] {
    let (s, c) = phi.sin_cos();
    let t = s / c;
    let rc = r * c;
    [
        d[0],
        d[1] / rc,
        d[2] / r,
        d[3],
        (d[4] + rc * c * d[3] - s * c * d[2]) / (rc * rc),
        (d[5] + t * d[1]) / (r * rc),
        (d[6] + r * d[3]) / (r * r),
        (d[7] - d[1] / r) / rc,
        (d[8] - d[2] / r) / r,
        d[9],
    ]
}

#[inline]
fn cart_to_local(cart: &[f64; 10], lam: f64, phi: f64) -> [f64; 10] {
    let (sl, cl) = lam.sin_cos();
    let (sp, cp) = phi.sin_cos();
    let rot = [[-sl, cl, 0.0], [-sp * cl, -sp * sl, cp], [cp * cl, cp * sl, sp]];
    let g = [cart[1], cart[2], cart[3]];
    let tt = [[cart[4], cart[5], cart[6]], [cart[5], cart[7], cart[8]], [cart[6], cart[8], cart[9]]];
    let mut gl = [0.0; 3];
    let mut rt = [[0.0; 3]; 3];
    for i in 0..3 {
        for a in 0..3 {
            gl[i] += rot[i][a] * g[a];
            for b in 0..3 {
                rt[i][b] += rot[i][a] * tt[a][b];
            }
        }
    }
    let mut tl = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            for b in 0..3 {
                tl[i][j] += rt[i][b] * rot[j][b];
            }
        }
    }
    [cart[0], gl[0], gl[1], gl[2], tl[0][0], tl[0][1], tl[1][1], tl[0][2], tl[1][2], tl[2][2]]
}

// ------------------------------------------------------------------ 1a. LEGENDRE_a
struct LegA {
    cosm: Vec<f64>,
    sinm: Vec<f64>,
    mm: Vec<f64>,
    mt: Vec<f64>,
    d2c: Vec<f64>,
    k2t: Vec<f64>,
    p: Vec<f64>,
    p1: Vec<f64>,
}

impl LegA {
    fn new(n: usize) -> Self {
        let z = |k: usize| vec![0.0f64; k];
        LegA {
            cosm: z(n + 1), sinm: z(n + 1), mm: (0..=n).map(|m| m as f64).collect(),
            mt: z(n + 1), d2c: z(n + 1), k2t: z(n + 1), p: z(n + 3), p1: z(n + 3),
        }
    }

    fn eval(&mut self, md: &Model, lam: f64, phi: f64, r: f64) -> [f64; 10] {
        let (nmax, n1, t) = (md.n, md.n1, &md.t);
        let l = t.l;
        let (u, c) = phi.sin_cos();
        let tn = u / c;
        let sec2 = 1.0 / (c * c);
        let q = md.a / r;
        for m in 0..=nmax {
            let fm = self.mm[m];
            let (s_, c_) = (fm * lam).sin_cos();
            self.cosm[m] = c_;
            self.sinm[m] = s_;
            let mt = fm * tn;
            self.mt[m] = mt;
            self.d2c[m] = mt * mt - fm * sec2;
            self.k2t[m] = (2.0 * fm + 1.0) * tn;
        }
        self.p.fill(0.0);
        self.p1.fill(0.0);
        let mut acc = [0.0f64; 10];
        let mut w = 1.0f64;
        for n in 0..=nmax {
            // before: p = row n-1, p1 = row n-2. After the swap p1 = row n-1 and p holds
            // row n-2, which is overwritten in place by row n.
            std::mem::swap(&mut self.p, &mut self.p1);
            let e = n + 1;
            {
                let (p, p1) = (&mut self.p, &self.p1);
                if n == 0 {
                    p[0] = 1.0;
                } else {
                    let g = &t.g[n * l..n * l + n];
                    let h = &t.h[n * l..n * l + n];
                    for m in 0..n {
                        // p[m] on the right-hand side is P_{n-2,m} (in-place update)
                        p[m] = g[m] * (u * p1[m]) - h[m] * p[m];
                    }
                    p[n] = t.f[n] * c * p1[n - 1];
                }
            }
            let p = &self.p;
            let kr = &t.k[n * l..n * l + e];
            let jr = &t.j[n * l..n * l + e];
            let cn = &md.c[n * n1..n * n1 + e];
            let sn = &md.s[n * n1..n * n1 + e];
            let (mut px, mut mpy, mut m2px, mut dpx, mut mdpy, mut d2px) = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
            for m in 0..e {
                let (p0, pa, pb) = (p[m], p[m + 1], p[m + 2]);
                let kn = kr[m];
                let dp = kn * pa - self.mt[m] * p0;
                let d2p = self.d2c[m] * p0 - (kn * self.k2t[m]) * pa + jr[m] * pb;
                let (cm, sm) = (self.cosm[m], self.sinm[m]);
                let x = cn[m] * cm + sn[m] * sm;
                let y = sn[m] * cm - cn[m] * sm;
                let fm = self.mm[m];
                px += p0 * x;
                mpy += p0 * fm * y;
                m2px += p0 * fm * fm * x;
                dpx += dp * x;
                mdpy += dp * fm * y;
                d2px += d2p * x;
            }
            let nn1 = n as f64 + 1.0;
            let nn12 = nn1 * (n as f64 + 2.0);
            acc[0] += w * px;
            acc[1] += w * mpy;
            acc[2] += w * dpx;
            acc[3] -= nn1 * w * px;
            acc[4] -= w * m2px;
            acc[5] += w * mdpy;
            acc[6] += w * d2px;
            acc[7] -= nn1 * w * mpy;
            acc[8] -= nn1 * w * dpx;
            acc[9] += nn12 * w * px;
            w *= q;
        }
        let base = md.gm / r;
        let mut d = [0.0; 10];
        for i in 0..10 {
            d[i] = acc[i] * base;
        }
        d[3] /= r;
        d[7] /= r;
        d[8] /= r;
        d[9] /= r * r;
        spherical_to_local(&d, r, phi)
    }
}

// ------------------------------------------------------------------ 1b. LEGENDRE_b
// The same FCR, but organised as the Matlab + DLL programme: for every point the whole
// triangle of fnALFs and of their derivatives is built and stored in freshly allocated
// arrays, and the sums are done afterwards as separate whole-array passes. (The interpreter
// overhead of the original is of course not reproduced by compiled code.)
fn legendre_b(md: &Model, lam: f64, phi: f64, r: f64) -> [f64; 10] {
    let (nmax, n1, t) = (md.n, md.n1, &md.t);
    let l = t.l;
    let w3 = nmax + 3;
    let (u, c) = phi.sin_cos();
    let tn = u / c;
    let mut p = vec![0.0f64; n1 * w3];
    p[0] = 1.0;
    if nmax >= 1 {
        p[w3] = t.g[l] * u;
        p[w3 + 1] = t.f[1] * c;
    }
    for n in 2..=nmax {
        for m in 0..n {
            p[n * w3 + m] = t.g[n * l + m] * u * p[(n - 1) * w3 + m] - t.h[n * l + m] * p[(n - 2) * w3 + m];
        }
        p[n * w3 + n] = t.f[n] * c * p[(n - 1) * w3 + n - 1];
    }
    let mut dp = vec![0.0f64; n1 * n1];
    let mut d2p = vec![0.0f64; n1 * n1];
    for n in 0..n1 {
        for m in 0..n1 {
            let fm = m as f64;
            let (p0, pa, pb) = (p[n * w3 + m], p[n * w3 + m + 1], p[n * w3 + m + 2]);
            let kt = t.k[n * l + m];
            dp[n * n1 + m] = kt * pa - fm * tn * p0;
            d2p[n * n1 + m] = (fm * fm * tn * tn - fm / (c * c)) * p0 - (2.0 * fm + 1.0) * kt * tn * pa
                + t.j[n * l + m] * pb;
        }
    }
    let cm: Vec<f64> = (0..n1).map(|m| (m as f64 * lam).cos()).collect();
    let sm: Vec<f64> = (0..n1).map(|m| (m as f64 * lam).sin()).collect();
    let mut x = vec![0.0f64; n1 * n1];
    let mut y = vec![0.0f64; n1 * n1];
    for n in 0..n1 {
        for m in 0..n1 {
            let (cc, ss) = (md.c[n * n1 + m], md.s[n * n1 + m]);
            x[n * n1 + m] = cc * cm[m] + ss * sm[m];
            y[n * n1 + m] = ss * cm[m] - cc * sm[m];
        }
    }
    let mut rows = vec![[0.0f64; 6]; n1]; // PX, mPY, m2PX, dPX, mdPY, d2PX per degree
    for n in 0..n1 {
        let mut s6 = [0.0f64; 6];
        for m in 0..n1 {
            let fm = m as f64;
            let p0 = p[n * w3 + m];
            let (xx, yy) = (x[n * n1 + m], y[n * n1 + m]);
            s6[0] += p0 * xx;
            s6[1] += p0 * fm * yy;
            s6[2] += p0 * fm * fm * xx;
            s6[3] += dp[n * n1 + m] * xx;
            s6[4] += dp[n * n1 + m] * fm * yy;
            s6[5] += d2p[n * n1 + m] * xx;
        }
        rows[n] = s6;
    }
    let q = md.a / r;
    let mut d = [0.0f64; 10];
    for n in 0..n1 {
        let w = q.powi(n as i32);
        let nn1 = n as f64 + 1.0;
        let [px, mpy, m2px, dpx, mdpy, d2px] = rows[n];
        d[0] += w * px;
        d[1] += w * mpy;
        d[2] += w * dpx;
        d[3] -= nn1 * w * px;
        d[4] -= w * m2px;
        d[5] += w * mdpy;
        d[6] += w * d2px;
        d[7] -= nn1 * w * mpy;
        d[8] -= nn1 * w * dpx;
        d[9] += nn1 * (nn1 + 1.0) * w * px;
    }
    let base = md.gm / r;
    for v in d.iter_mut() {
        *v *= base;
    }
    d[3] /= r;
    d[7] /= r;
    d[8] /= r;
    d[9] /= r * r;
    spherical_to_local(&d, r, phi)
}

// ------------------------------------------------------------------ 2. PINES
struct Pines {
    rm: Vec<f64>, // (s + i t)^m stored at index m + 2
    im: Vec<f64>,
    rows: Vec<Vec<f64>>, // 4 rolling rows of Helmholtz polynomials
}

impl Pines {
    fn new(n: usize) -> Self {
        Pines { rm: vec![0.0; n + 3], im: vec![0.0; n + 3], rows: (0..4).map(|_| vec![0.0; n + 5]).collect() }
    }

    #[inline]
    fn hrow(&mut self, t: &Tables, n: usize, u: f64) {
        let l = t.l;
        let mut an = std::mem::take(&mut self.rows[n % 4]);
        {
            let a1 = &self.rows[(n + 3) % 4];
            let a2 = &self.rows[(n + 2) % 4];
            let g = &t.g[n * l..n * l + n];
            let h = &t.h[n * l..n * l + n];
            for m in 0..n {
                an[m] = g[m] * (u * a1[m]) - h[m] * a2[m];
            }
            an[n] = t.f[n] * a1[n - 1];
        }
        self.rows[n % 4] = an;
    }

    /// Cartesian quantities (CART order).
    fn eval(&mut self, md: &Model, lam: f64, phi: f64, r: f64) -> [f64; 10] {
        let (nmax, n1, t) = (md.n, md.n1, &md.t);
        let l = t.l;
        let (sp, cp) = phi.sin_cos();
        let (sl, cl) = lam.sin_cos();
        let (s, tt, u) = (cp * cl, cp * sl, sp);
        let q = md.a / r;
        self.rm.fill(0.0);
        self.im.fill(0.0);
        self.rm[2] = 1.0;
        for m in 1..=nmax {
            self.rm[m + 2] = s * self.rm[m + 1] - tt * self.im[m + 1];
            self.im[m + 2] = s * self.im[m + 1] + tt * self.rm[m + 1];
        }
        for row in self.rows.iter_mut() {
            row.fill(0.0);
        }
        self.rows[0][0] = 1.0;
        self.hrow(t, 1, u);
        self.hrow(t, 2, u);
        let mut acc = [0.0f64; 14];
        let mut w = 1.0f64;
        for n in 0..=nmax {
            let e = n + 1;
            let a0 = &self.rows[n % 4];
            let b1 = &self.rows[(n + 1) % 4];
            let b2 = &self.rows[(n + 2) % 4];
            let kr = &t.k[n * l..n * l + e];
            let jr = &t.j[n * l..n * l + e];
            let r11 = &t.r11[n * l..n * l + e];
            let r12 = &t.r12[n * l..n * l + e];
            let r22 = &t.r22[n * l..n * l + e];
            let cn = &md.c[n * n1..n * n1 + e];
            let sn = &md.s[n * n1..n * n1 + e];
            let mut o = [0.0f64; 14];
            for m in 0..e {
                let fm = m as f64;
                let (c_, s_) = (cn[m], sn[m]);
                let (rm, im) = (self.rm[m + 2], self.im[m + 2]);
                let (rm1, im1) = (self.rm[m + 1], self.im[m + 1]);
                let (rm2, im2) = (self.rm[m], self.im[m]);
                let d0 = c_ * rm + s_ * im; // D_nm
                let v0 = a0[m];
                let v1 = kr[m] * a0[m + 1];
                let v2 = -r11[m] * b1[m + 1];
                let v3 = jr[m] * a0[m + 2];
                let v4 = -r12[m] * b1[m + 2];
                let v5 = r22[m] * b2[m + 2];
                o[0] += v0 * d0;
                o[1] += v1 * d0;
                o[2] += v2 * d0;
                o[3] += v3 * d0;
                o[4] += v4 * d0;
                o[5] += v5 * d0;
                let e1 = fm * (c_ * rm1 + s_ * im1); // m E_nm
                let f1 = fm * (s_ * rm1 - c_ * im1); // m F_nm
                o[6] += v0 * e1;
                o[7] += v0 * f1;
                o[8] += v1 * e1;
                o[9] += v1 * f1;
                o[10] += v2 * e1;
                o[11] += v2 * f1;
                let mm1 = fm * (fm - 1.0);
                o[12] += v0 * mm1 * (c_ * rm2 + s_ * im2);
                o[13] += v0 * mm1 * (s_ * rm2 - c_ * im2);
            }
            for i in 0..14 {
                acc[i] += w * o[i];
            }
            w *= q;
            if n + 3 <= nmax + 2 {
                self.hrow(t, n + 3, u);
            }
        }
        let mu = md.gm;
        let r1 = mu / (r * r);
        let r2 = r1 / r;
        let v = acc[0] * (mu / r);
        let (a3, a4) = (acc[1] * r1, acc[2] * r1);
        let (a33, a34, a44) = (acc[3] * r2, acc[4] * r2, acc[5] * r2);
        let (a1, a2) = (acc[6] * r1, acc[7] * r1);
        let (a13, a23, a14, a24) = (acc[8] * r2, acc[9] * r2, acc[10] * r2, acc[11] * r2);
        let (a11, a12) = (acc[12] * r2, acc[13] * r2);
        let xh = [s, tt, u];
        let ai4 = [a14, a24, a34];
        let aij = [[a11, a12, a13], [a12, -a11, a23], [a13, a23, a33]];
        let a4r = a4 / r;
        let mut out = [0.0f64; 10];
        out[0] = v;
        out[1] = a1 + s * a4;
        out[2] = a2 + tt * a4;
        out[3] = a3 + u * a4;
        let pairs = [(0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2)];
        for (col, &(i, j)) in pairs.iter().enumerate() {
            let mut vij = aij[i][j] + xh[j] * ai4[i] + xh[i] * ai4[j] + xh[i] * xh[j] * a44;
            if i == j {
                vij += a4r;
            }
            out[4 + col] = vij;
        }
        out
    }
}

// ------------------------------------------------------------------ 2b. GOTTLIEB
// Pines' Helmholtz polynomials in Gottlieb's (1993, NASA CR-188243) normalised form: the
// Cartesian derivatives use d/du Abar_nm = k_nm Abar_n,m+1, d2/du2 Abar_nm = j_nm Abar_n,m+2,
// so only row n is needed (3 rolling rows, two zero columns past the diagonal).
struct Gottlieb {
    rm: Vec<f64>,
    im: Vec<f64>,
    rows: Vec<Vec<f64>>,
}

impl Gottlieb {
    fn new(n: usize) -> Self {
        Gottlieb { rm: vec![0.0; n + 3], im: vec![0.0; n + 3], rows: (0..3).map(|_| vec![0.0; n + 3]).collect() }
    }

    /// Cartesian quantities (CART order).
    fn eval(&mut self, md: &Model, lam: f64, phi: f64, r: f64) -> [f64; 10] {
        let (nmax, n1, t) = (md.n, md.n1, &md.t);
        let l = t.l;
        let (sp, cp) = phi.sin_cos();
        let (sl, cl) = lam.sin_cos();
        let (s, tt, u) = (cp * cl, cp * sl, sp);
        let q = md.a / r;
        self.rm.fill(0.0);
        self.im.fill(0.0);
        self.rm[2] = 1.0;
        for m in 1..=nmax {
            self.rm[m + 2] = s * self.rm[m + 1] - tt * self.im[m + 1];
            self.im[m + 2] = s * self.im[m + 1] + tt * self.rm[m + 1];
        }
        for row in self.rows.iter_mut() {
            row.fill(0.0);
        }
        self.rows[0][0] = 1.0;
        // v gam l j k s t n o h p q r M
        let mut acc = [0.0f64; 14];
        let mut w = 1.0f64;
        for n in 0..=nmax {
            if n >= 1 {
                let mut an = std::mem::take(&mut self.rows[n % 3]);
                {
                    let a1 = &self.rows[(n + 2) % 3];
                    let a2 = &self.rows[(n + 1) % 3];
                    let g = &t.g[n * l..n * l + n];
                    let h = &t.h[n * l..n * l + n];
                    for m in 0..n {
                        an[m] = g[m] * (u * a1[m]) - h[m] * a2[m];
                    }
                    an[n] = t.f[n] * a1[n - 1];
                }
                self.rows[n % 3] = an;
            }
            let e = n + 1;
            let a0 = &self.rows[n % 3];
            let kr = &t.k[n * l..n * l + e];
            let jr = &t.j[n * l..n * l + e];
            let cn = &md.c[n * n1..n * n1 + e];
            let sn = &md.s[n * n1..n * n1 + e];
            let nf = n as f64;
            let mut o = [0.0f64; 14];
            for m in 0..e {
                let fm = m as f64;
                let (c_, s_) = (cn[m], sn[m]);
                let (rm, im) = (self.rm[m + 2], self.im[m + 2]);
                let (rm1, im1) = (self.rm[m + 1], self.im[m + 1]);
                let (rm2, im2) = (self.rm[m], self.im[m]);
                let b0 = c_ * rm + s_ * im; //   C r_m + S i_m
                let b1 = c_ * rm1 + s_ * im1; //  C r_m-1 + S i_m-1
                let a1 = s_ * rm1 - c_ * im1; // -(C i_m-1 - S r_m-1)
                let np1 = nf + fm + 1.0;
                let p0 = a0[m];
                let pb = p0 * b0;
                o[0] += pb;
                o[1] += np1 * pb;
                o[2] += np1 * (np1 + 1.0) * pb;
                let mp = fm * p0;
                let (mb, ma) = (mp * b1, mp * a1);
                o[3] += mb;
                o[4] += ma;
                o[5] += np1 * mb;
                o[6] += np1 * ma;
                let mm1 = (fm - 1.0) * mp;
                o[7] += mm1 * (c_ * rm2 + s_ * im2);
                o[8] += mm1 * (c_ * im2 - s_ * rm2);
                let zg = kr[m] * a0[m + 1];
                let zb = zg * b0;
                o[9] += zb;
                o[10] += np1 * zb;
                o[11] += fm * zg * b1;
                o[12] += fm * zg * a1;
                o[13] += jr[m] * a0[m + 2] * b0;
            }
            for i in 0..14 {
                acc[i] += w * o[i];
            }
            w *= q;
        }
        let [sv, sgam, sl_, sj, sk, ss, st, sn_, so, sh, spp, sq, sr, sm] = acc;
        let lam_ = sgam + u * sh;
        let d1 = ss + u * sq;
        let d2 = st + u * sr;
        let gg = -(spp + sh + u * sm);
        let ff = sl_ + lam_ + u * (spp + sh - gg);
        let mu = md.gm;
        let (m1, m2) = (mu / r, mu / (r * r));
        let m3 = m2 / r;
        [
            m1 * sv,
            m2 * (sj - lam_ * s),
            m2 * (sk - lam_ * tt),
            m2 * (sh - lam_ * u),
            m3 * ((ff * s - 2.0 * d1) * s - lam_ + sn_),
            m3 * ((ff * tt - d2) * s - d1 * tt - so),
            m3 * ((ff * s - d1) * u + gg * s + sq),
            m3 * ((ff * tt - 2.0 * d2) * tt - lam_ - sn_),
            m3 * ((ff * tt - d2) * u + gg * tt + sr),
            m3 * ((ff * u + 2.0 * gg) * u - lam_ + sm),
        ]
    }
}

// ------------------------------------------------------------------ 3. METRIS
struct Metris {
    l: usize,
    er: Vec<f64>, // [(n * l + m) * 10 + q]
    ei: Vec<f64>,
    yr: Vec<Vec<f64>>,
    yi: Vec<Vec<f64>>,
}

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

impl Metris {
    fn new(md: &Model) -> Self {
        let (n, n1, a, t) = (md.n, md.n1, md.a, &md.t);
        let l = n + 3;
        let mut kr = vec![0.0f64; l * l];
        let mut ki = vec![0.0f64; l * l];
        let sc = md.gm / a;
        for i in 0..n1 {
            for m in 0..n1 {
                kr[i * l + m] = sc * md.c[i * n1 + m];
                ki[i * l + m] = -sc * md.s[i * n1 + m];
            }
        }
        let ex = metris_op(&kr, &ki, 'x', t, l, a);
        let ey = metris_op(&kr, &ki, 'y', t, l, a);
        let ez = metris_op(&kr, &ki, 'z', t, l, a);
        let sets = vec![
            metris_op(&ex.0, &ex.1, 'x', t, l, a),
            metris_op(&ex.0, &ex.1, 'y', t, l, a),
            metris_op(&ex.0, &ex.1, 'z', t, l, a),
            metris_op(&ey.0, &ey.1, 'y', t, l, a),
            metris_op(&ey.0, &ey.1, 'z', t, l, a),
            metris_op(&ez.0, &ez.1, 'z', t, l, a),
        ];
        let mut all = vec![(kr, ki), ex, ey, ez];
        all.extend(sets);
        let mut er = vec![0.0f64; l * l * 10];
        let mut ei = vec![0.0f64; l * l * 10];
        for (qi, (a_r, a_i)) in all.iter().enumerate() {
            for i in 0..l * l {
                er[i * 10 + qi] = a_r[i];
                ei[i * 10 + qi] = a_i[i];
            }
        }
        Metris {
            l, er, ei,
            yr: (0..3).map(|_| vec![0.0; l]).collect(),
            yi: (0..3).map(|_| vec![0.0; l]).collect(),
        }
    }

    /// Cartesian quantities (CART order).
    fn eval(&mut self, md: &Model, lam: f64, phi: f64, r: f64) -> [f64; 10] {
        let t = &md.t;
        let (lt, l) = (t.l, self.l);
        let (sp, cp) = phi.sin_cos();
        let (sl, cl) = lam.sin_cos();
        let a = md.a;
        let (xp, yp, zp, rp) = (r * cp * cl / a, r * cp * sl / a, r * sp / a, r / a);
        let r2i = 1.0 / (rp * rp);
        let (zr, zeta_r, zeta_i) = (zp * r2i, xp * r2i, yp * r2i);
        for k in 0..3 {
            self.yr[k].fill(0.0);
            self.yi[k].fill(0.0);
        }
        let mut f = [0.0f64; 10];
        for n in 0..l {
            let mut r0 = std::mem::take(&mut self.yr[n % 3]);
            let mut i0 = std::mem::take(&mut self.yi[n % 3]);
            if n == 0 {
                r0[0] = 1.0 / rp;
            } else {
                let r1 = &self.yr[(n + 2) % 3];
                let i1 = &self.yi[(n + 2) % 3];
                let r2 = &self.yr[(n + 1) % 3];
                let i2 = &self.yi[(n + 1) % 3];
                let g = &t.g[n * lt..n * lt + n];
                let h = &t.h[n * lt..n * lt + n];
                for m in 0..n {
                    r0[m] = g[m] * (zr * r1[m]) - h[m] * (r2i * r2[m]);
                    i0[m] = g[m] * (zr * i1[m]) - h[m] * (r2i * i2[m]);
                }
                let (ar, ai) = (r1[n - 1], i1[n - 1]);
                r0[n] = t.f[n] * (zeta_r * ar - zeta_i * ai);
                i0[n] = t.f[n] * (zeta_r * ai + zeta_i * ar);
            }
            let e = n + 1;
            let er = &self.er[n * l * 10..(n * l + e) * 10];
            let ei = &self.ei[n * l * 10..(n * l + e) * 10];
            // partial sum per degree, then into the total: N+3 additions into f instead of
            // ~N^2/2 (a single running accumulator costs about one digit at N = 360)
            let mut s = [0.0f64; 10];
            for m in 0..e {
                let (yr_, yi_) = (r0[m], i0[m]);
                let cr = &er[m * 10..m * 10 + 10];
                let ci = &ei[m * 10..m * 10 + 10];
                for qi in 0..10 {
                    s[qi] += yr_ * cr[qi] - yi_ * ci[qi];
                }
            }
            for qi in 0..10 {
                f[qi] += s[qi];
            }
            self.yr[n % 3] = r0;
            self.yi[n % 3] = i0;
        }
        f
    }
}

// ------------------------------------------------------------------ 4. CLENSHAW
/// Clenshaw recursions over the degree (n = N .. m) for one order, cosine (C) and sine (S)
/// coefficient families side by side (two lanes): [y, y_r, y_rr, y_u, y_uu, y_ur][p] at n = m.
/// g, h and the coefficients start at degree m.
#[inline(always)]
fn clenshaw_column(g: &[f64], h: &[f64], cc: &[f64], cs: &[f64], m: usize, uq: f64, q2: f64, q: f64) -> [[f64; 2]; 6] {
    let len = cc.len();
    let cs = &cs[..len];
    let (g, h) = (&g[..len + 2], &h[..len + 2]);
    let mut a = [[0.0f64; 2]; 6];
    let mut b = [[0.0f64; 2]; 6];
    for i in (0..len).rev() {
        let gn = g[i + 1];
        let aa = uq * gn;
        let bb = q2 * h[i + 2];
        let qg = q * gn;
        let nn1 = (m + i) as f64 + 1.0;
        let nn2 = nn1 * (nn1 + 1.0);
        let cn = [cc[i], cs[i]];
        let mut y = [[0.0f64; 2]; 6];
        // the terms that do not involve y_{n+1} of the same recursion are grouped, so that
        // the loop-carried dependency is one multiply-add per degree
        for p in 0..2 {
            y[0][p] = aa * a[0][p] + (cn[p] - bb * b[0][p]);
            y[1][p] = aa * a[1][p] + (nn1 * cn[p] - bb * b[1][p]);
            y[2][p] = aa * a[2][p] + (nn2 * cn[p] - bb * b[2][p]);
            y[3][p] = aa * a[3][p] + (qg * a[0][p] - bb * b[3][p]);
            y[4][p] = aa * a[4][p] + ((2.0 * qg) * a[3][p] - bb * b[4][p]);
            y[5][p] = aa * a[5][p] + (qg * a[1][p] - bb * b[5][p]);
        }
        b = a;
        a = y;
    }
    a
}

/// v-terms of Table 5 for one point / parallel: vt[(i * 2 + p) * (N+1) + m],
/// i = v, v_u, v_r, v_uu, v_ur, v_rr; p = cosine (C) / sine (S) part.
fn clenshaw_vterms(md: &Model, u: f64, c: f64, q: f64, vt: &mut [f64]) {
    let (nmax, n1, t) = (md.n, md.n1, &md.t);
    let l = t.l;
    let uq = u * q;
    let q2 = q * q;
    let c2 = c * c;
    let cq = c * q;
    let mut pm = 1.0f64;
    for m in 0..=nmax {
        let gt = &t.gt[m * l..(m + 1) * l];
        let ht = &t.ht[m * l..(m + 1) * l];
        let ct = &md.ct[m * n1..(m + 1) * n1];
        let st = &md.st[m * n1..(m + 1) * n1];
        let y1 = clenshaw_column(&gt[m..nmax + 3], &ht[m..nmax + 3], &ct[m..=nmax], &st[m..=nmax], m, uq, q2, q);
        if m >= 1 {
            pm *= t.f[m] * cq; // Pbar_mm (a/r)^m
        }
        let fm = m as f64;
        let dp = (-fm * u / c2) * pm;
        let d2p = (-fm / c2 + fm * (fm - 2.0) * u * u / (c2 * c2)) * pm;
        for p in 0..2 {
            let (y0, yr, yrr, yu, yuu, yur) = (y1[0][p], y1[1][p], y1[2][p], y1[3][p], y1[4][p], y1[5][p]);
            vt[p * n1 + m] = y0 * pm;
            vt[(2 + p) * n1 + m] = y0 * dp + yu * pm;
            vt[(4 + p) * n1 + m] = yr * pm;
            vt[(6 + p) * n1 + m] = yuu * pm + y0 * d2p + 2.0 * yu * dp;
            vt[(8 + p) * n1 + m] = yr * dp + yur * pm;
            vt[(10 + p) * n1 + m] = yrr * pm;
        }
    }
}

/// Sums over the orders (Table 6) for one longitude, then the local frame.
#[inline]
fn clenshaw_sums(md: &Model, vt: &[f64], cm: &[f64], sm: &[f64], r: f64, phi: f64) -> [f64; 10] {
    let n1 = md.n1;
    let mut sc = [0.0f64; 6];
    let mut ss = [0.0f64; 3];
    let mut sc2 = 0.0f64;
    for m in 0..n1 {
        let (c_, s_) = (cm[m], sm[m]);
        let fm = m as f64;
        let pc0 = vt[m] * c_ + vt[n1 + m] * s_;
        sc[0] += pc0;
        sc2 += fm * fm * pc0;
        for i in 1..6 {
            sc[i] += vt[2 * i * n1 + m] * c_ + vt[(2 * i + 1) * n1 + m] * s_;
        }
        for i in 0..3 {
            ss[i] += fm * (vt[2 * i * n1 + m] * s_ - vt[(2 * i + 1) * n1 + m] * c_);
        }
    }
    let (u, c) = phi.sin_cos();
    let b = md.gm / r;
    let d = [
        b * sc[0],
        -b * ss[0],
        b * c * sc[1],
        -b / r * sc[2],
        -b * sc2,
        -b * c * ss[1],
        b * c * c * sc[3] - b * u * sc[1],
        b / r * ss[2],
        -b * c / r * sc[4],
        b / (r * r) * sc[5],
    ];
    spherical_to_local(&d, r, phi)
}

// ------------------------------------------------------------------ main
fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() < 5 {
        eprintln!("usage: sph_rs <METHOD> <coef.bin> <points.bin> <out.bin> [--reps R] [--nmax n] [--grid]");
        std::process::exit(2);
    }
    let method = args[1].as_str();
    let (mut reps, mut nmax, mut grid) = (1usize, None, false);
    let mut i = 5;
    while i < args.len() {
        match args[i].as_str() {
            "--reps" => { reps = args[i + 1].parse().unwrap(); i += 1; }
            "--nmax" => { nmax = Some(args[i + 1].parse().unwrap()); i += 1; }
            "--grid" => grid = true,
            other => { eprintln!("unknown option {other}"); std::process::exit(2); }
        }
        i += 1;
    }
    let reps = reps.max(1);

    let t_setup = Instant::now();
    let md = load_model(&args[2], nmax);
    let n1 = md.n1;
    let b = fs::read(&args[3]).expect("cannot read point file");
    let mut off = 0;
    let (lam, phi, rad, lon1, lat1, rgrid, k);
    if grid {
        let nlon = rd_u64(&b, &mut off);
        let nlat = rd_u64(&b, &mut off);
        rgrid = rd_f64(&b, &mut off);
        lon1 = rd_vec(&b, &mut off, nlon);
        lat1 = rd_vec(&b, &mut off, nlat);
        k = nlon * nlat;
        lam = vec![];
        phi = vec![];
        rad = vec![];
    } else {
        k = rd_u64(&b, &mut off);
        lam = rd_vec(&b, &mut off, k);
        phi = rd_vec(&b, &mut off, k);
        rad = rd_vec(&b, &mut off, k);
        lon1 = vec![];
        lat1 = vec![];
        rgrid = 0.0;
    }
    let has_cart = !grid && (method == "PINES" || method == "GOTTLIEB" || method == "METRIS");
    let rows = if has_cart { 20 } else { 10 };
    let mut out = vec![0.0f64; rows * k];
    let mut metris = if method == "METRIS" { Some(Metris::new(&md)) } else { None };
    let setup = t_setup.elapsed().as_secs_f64();

    let (mut best_cpu, mut best_wall) = (f64::INFINITY, f64::INFINITY);
    for _ in 0..reps {
        let (c0, w0) = (cpu_time(), Instant::now());
        if grid {
            if method != "CLENSHAW" {
                eprintln!("--grid is only available for CLENSHAW");
                std::process::exit(2);
            }
            let nlon = lon1.len();
            // cos(m lambda), sin(m lambda) for every meridian (once for the whole grid)
            let mut cm = vec![0.0f64; nlon * n1];
            let mut sm = vec![0.0f64; nlon * n1];
            for (jl, &lo) in lon1.iter().enumerate() {
                for m in 0..n1 {
                    let (s_, c_) = (m as f64 * lo).sin_cos();
                    cm[jl * n1 + m] = c_;
                    sm[jl * n1 + m] = s_;
                }
            }
            let mut vt = vec![0.0f64; 12 * n1];
            let q = md.a / rgrid;
            for (il, &la) in lat1.iter().enumerate() {
                let (u, c) = la.sin_cos();
                clenshaw_vterms(&md, u, c, q, &mut vt);
                for jl in 0..nlon {
                    let v = clenshaw_sums(&md, &vt, &cm[jl * n1..(jl + 1) * n1],
                                          &sm[jl * n1..(jl + 1) * n1], rgrid, la);
                    let p = il * nlon + jl;
                    for qi in 0..10 {
                        out[qi * k + p] = v[qi];
                    }
                }
            }
        } else {
            match method {
                "LEGENDRE_a" => {
                    let mut w = LegA::new(md.n);
                    for p in 0..k {
                        let v = w.eval(&md, lam[p], phi[p], rad[p]);
                        for qi in 0..10 {
                            out[qi * k + p] = v[qi];
                        }
                    }
                }
                "LEGENDRE_b" => {
                    for p in 0..k {
                        let v = legendre_b(&md, lam[p], phi[p], rad[p]);
                        for qi in 0..10 {
                            out[qi * k + p] = v[qi];
                        }
                    }
                }
                "PINES" => {
                    let mut w = Pines::new(md.n);
                    for p in 0..k {
                        let cart = w.eval(&md, lam[p], phi[p], rad[p]);
                        let v = cart_to_local(&cart, lam[p], phi[p]);
                        for qi in 0..10 {
                            out[qi * k + p] = v[qi];
                            out[(10 + qi) * k + p] = cart[qi];
                        }
                    }
                }
                "GOTTLIEB" => {
                    let mut w = Gottlieb::new(md.n);
                    for p in 0..k {
                        let cart = w.eval(&md, lam[p], phi[p], rad[p]);
                        let v = cart_to_local(&cart, lam[p], phi[p]);
                        for qi in 0..10 {
                            out[qi * k + p] = v[qi];
                            out[(10 + qi) * k + p] = cart[qi];
                        }
                    }
                }
                "METRIS" => {
                    let w = metris.as_mut().unwrap();
                    for p in 0..k {
                        let cart = w.eval(&md, lam[p], phi[p], rad[p]);
                        let v = cart_to_local(&cart, lam[p], phi[p]);
                        for qi in 0..10 {
                            out[qi * k + p] = v[qi];
                            out[(10 + qi) * k + p] = cart[qi];
                        }
                    }
                }
                "CLENSHAW" => {
                    let mut vt = vec![0.0f64; 12 * n1];
                    let mut cm = vec![0.0f64; n1];
                    let mut sm = vec![0.0f64; n1];
                    for p in 0..k {
                        let (u, c) = phi[p].sin_cos();
                        clenshaw_vterms(&md, u, c, md.a / rad[p], &mut vt);
                        for m in 0..n1 {
                            let (s_, c_) = (m as f64 * lam[p]).sin_cos();
                            cm[m] = c_;
                            sm[m] = s_;
                        }
                        let v = clenshaw_sums(&md, &vt, &cm, &sm, rad[p], phi[p]);
                        for qi in 0..10 {
                            out[qi * k + p] = v[qi];
                        }
                    }
                }
                other => {
                    eprintln!("unknown method {other}");
                    std::process::exit(2);
                }
            }
        }
        best_cpu = best_cpu.min(cpu_time() - c0);
        best_wall = best_wall.min(w0.elapsed().as_secs_f64());
    }

    let mut bytes = Vec::with_capacity(16 + 8 * out.len());
    bytes.extend_from_slice(&(rows as u64).to_le_bytes());
    bytes.extend_from_slice(&(k as u64).to_le_bytes());
    for v in &out {
        bytes.extend_from_slice(&v.to_le_bytes());
    }
    fs::File::create(&args[4]).and_then(|mut fh| fh.write_all(&bytes)).expect("cannot write output");
    println!(
        "{{\"cpu\": {:.9e}, \"wall\": {:.9e}, \"setup\": {:.6e}, \"K\": {}, \"N\": {}, \"reps\": {}}}",
        best_cpu, best_wall, setup, k, md.n, reps
    );
}
