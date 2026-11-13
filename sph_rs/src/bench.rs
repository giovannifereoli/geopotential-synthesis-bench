//! bench.rs -- fair CPU-time benchmark of the synthesis methods, League A (scattered points,
//! scalar code). Design and rules: README.md (timing protocol).
//!
//!   sph_rs bench <coef.bin> <xyz.bin> [--nmax N] [--level 1|2] [--blocks B] [--target-ms T]
//!                [--seed S] [--methods M1,M2,..]
//!   sph_rs eval  <coef.bin> <xyz.bin> <out.bin> --method M [--nmax N] [--level 1|2]
//!
//! xyz.bin : u64 K, x[K], y[K], z[K]  body-fixed Cartesian positions [m]
//! out.bin : u64 rows, u64 K, data[rows * K] row-major, rows = 3 (level 1) or 9 (level 2)
//!
//! One contract for every method:
//!   input  x, y, z (body-fixed);
//!   output level 1 = acceleration (ax, ay, az); level 2 = acceleration + gravity-gradient
//!          tensor (xx, xy, xz, yy, yz, zz), Cartesian body-fixed. The potential is not output.
//!   level chosen at compile time (const generic) so that level 1 does none of the level-2 work;
//!   no Laplace shortcut (all six tensor components are synthesised);
//!   cos(m lambda), sin(m lambda) by the same complex recursion in every spherical method;
//!   (a/r)^n by a running product; partial sums per degree; no allocation in the timed loop;
//!   no hand-written mul_add, no fast-math, no flush-to-zero.
//!
//! The spherical methods are pole-safe (no division by cos(phi)): LEGENDRE_a/b and CLENSHAW
//! recur on Z_nm = Pbar_nm / cos(phi)^min(m,2), so that the 1/cos and 1/cos^2 of the local-frame
//! gradient and tensor cancel analytically (same recursion, ~2 extra multiplies per term);
//! LEAR uses Pbar_nm / cos(phi) at level 1 and adds Pbar_nm / cos(phi)^2 at level 2 (see LEAR).

use crate::{clock_gettime, load_model, metris_op, rd_f64, rd_u64, Model, Tables, Timespec};
use std::fs;
use std::hint::black_box;
use std::io::Write;
use std::time::Instant;

// ------------------------------------------------------------------ timer, QoS
#[cfg(target_os = "macos")]
const CLOCK_THREAD_CPUTIME_ID: i32 = 16;
#[cfg(not(target_os = "macos"))]
const CLOCK_THREAD_CPUTIME_ID: i32 = 3;

fn thread_cpu() -> f64 {
    let mut ts = Timespec { tv_sec: 0, tv_nsec: 0 };
    unsafe {
        clock_gettime(CLOCK_THREAD_CPUTIME_ID, &mut ts);
    }
    ts.tv_sec as f64 + 1e-9 * ts.tv_nsec as f64
}

#[cfg(target_os = "macos")]
extern "C" {
    fn pthread_set_qos_class_self_np(qos: u32, rel: i32) -> i32;
}
/// macOS: QOS_CLASS_USER_INTERACTIVE keeps the thread on a performance core.
#[cfg(target_os = "macos")]
fn set_qos() -> i32 {
    unsafe { pthread_set_qos_class_self_np(0x21, 0) }
}
#[cfg(not(target_os = "macos"))]
fn set_qos() -> i32 {
    -1
}

/// Fixed FMA-free multiply-add chains (8 independent): ns per iteration. Run before every
/// measurement; a drift of the clock (thermal, efficiency core) shows up here.
#[inline(never)]
fn calibrate(iters: usize) -> f64 {
    let m = black_box(0.999_999_9f64);
    let b = black_box(1.0e-9f64);
    let mut v = black_box([1.0f64, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7]);
    let t0 = thread_cpu();
    for _ in 0..iters {
        for x in v.iter_mut() {
            *x = *x * m + b;
        }
        v = black_box(v);
    }
    let dt = thread_cpu() - t0;
    black_box(v);
    dt / iters as f64 * 1e9
}

struct Rng(u64);
impl Rng {
    fn next(&mut self) -> u64 {
        // xorshift64*
        self.0 ^= self.0 >> 12;
        self.0 ^= self.0 << 25;
        self.0 ^= self.0 >> 27;
        self.0.wrapping_mul(0x2545_F491_4F6C_DD1D)
    }
    fn shuffle<T>(&mut self, v: &mut [T]) {
        for i in (1..v.len()).rev() {
            let j = (self.next() % (i as u64 + 1)) as usize;
            v.swap(i, j);
        }
    }
}

// ------------------------------------------------------------------ shared pieces
trait Kernel {
    fn eval<const L2: bool>(&mut self, md: &Model, x: f64, y: f64, z: f64) -> [f64; 9];
}

/// r, u = sin(phi), c = cos(phi), cos(lambda), sin(lambda) from x, y, z (not at the pole itself).
#[inline(always)]
fn sph_geom(x: f64, y: f64, z: f64) -> (f64, f64, f64, f64, f64) {
    let rho2 = x * x + y * y;
    let rho = rho2.sqrt();
    let r = (rho2 + z * z).sqrt();
    (r, z / r, rho / r, x / rho, y / rho)
}

/// cos(m lambda), sin(m lambda), m = 0 .. len-1, by the complex recursion.
#[inline(always)]
fn trig(cm: &mut [f64], sm: &mut [f64], cl: f64, sl: f64) {
    let n = cm.len().min(sm.len());
    cm[0] = 1.0;
    sm[0] = 0.0;
    for m in 1..n {
        let (a, b) = (cm[m - 1], sm[m - 1]);
        cm[m] = a * cl - b * sl;
        sm[m] = b * cl + a * sl;
    }
}

/// Local north-oriented frame (e, n, r) -> Cartesian. t = [ee, en, er, nn, nr, rr].
#[inline(always)]
fn local_to_cart<const L2: bool>(g: [f64; 3], t: &[f64; 6], u: f64, c: f64, cl: f64, sl: f64) -> [f64; 9] {
    let rr = [[-sl, cl, 0.0], [-u * cl, -u * sl, c], [c * cl, c * sl, u]];
    let mut out = [0.0f64; 9];
    for a in 0..3 {
        out[a] = rr[0][a] * g[0] + rr[1][a] * g[1] + rr[2][a] * g[2];
    }
    if L2 {
        let tl = [[t[0], t[1], t[2]], [t[1], t[3], t[4]], [t[2], t[4], t[5]]];
        let mut mm = [[0.0f64; 3]; 3];
        for i in 0..3 {
            for b in 0..3 {
                mm[i][b] = tl[i][0] * rr[0][b] + tl[i][1] * rr[1][b] + tl[i][2] * rr[2][b];
            }
        }
        let pairs = [(0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2)];
        for (k, &(a, b)) in pairs.iter().enumerate() {
            out[3 + k] = rr[0][a] * mm[0][b] + rr[1][a] * mm[1][b] + rr[2][a] * mm[2][b];
        }
    }
    out
}

// ------------------------------------------------------------------ LEGENDRE (pole-safe)
// Per (n, m), with Z = Z_nm, Z1 = Z_n,m+1, Z2 = Z_n,m+2 and cp_m = c^min(m,2),
// cq_m = c^(min(m,2)-1) (cq_0 = 0):
//   P0 = cp_m Z = Pbar, Q0 = cq_m Z = Pbar/c, R0 = Z = Pbar/c^2 (m >= 2), P1 = cp_m+1 Z1, ...
//   dP  = k P1 - m u Q0                                   (= dPbar/dphi)
//   d2P = m(m-1) R0 - m^2 P0 - (2m+1) k u Q1 + j P2       (= d2Pbar/dphi2)
//   en  = m (k Q1 - (m-1) u R0)                           (V_lam,phi + tan V_lam) / c
//   ee  = -m(m-1) R0 - m P0 - u k Q1                      (V_lam,lam - u c V_phi) / c^2
#[inline(always)]
fn leg_term<const L2: bool>(o: &mut [f64; 6], fm: f64, p0: f64, q0: f64, r0: f64, p1: f64, q1: f64,
                            p2: f64, k: f64, j: f64, x: f64, y: f64, u: f64) {
    let dp = k * p1 - fm * u * q0;
    o[0] += p0 * x;
    o[1] += dp * x;
    o[2] += fm * q0 * y;
    if L2 {
        let mm1 = fm * (fm - 1.0);
        let kq1 = k * q1;
        let d2p = mm1 * r0 - fm * fm * p0 - (2.0 * fm + 1.0) * u * kq1 + j * p2;
        o[3] += d2p * x;
        o[4] += fm * (kq1 - (fm - 1.0) * u * r0) * y;
        o[5] += (-(mm1 * r0) - fm * p0 - u * kq1) * x;
    }
}

/// Sums over the order for one degree: [P X, dP X, m Q Y, d2P X, en Y, ee X].
#[inline(always)]
fn leg_row<const L2: bool>(z: &[f64], kr: &[f64], jr: &[f64], cn: &[f64], sn: &[f64], cm: &[f64],
                           sm: &[f64], e: usize, u: f64, c: f64, c2: f64) -> [f64; 6] {
    let z = &z[..e + 2];
    let (kr, jr, cn, sn, cm, sm) = (&kr[..e], &jr[..e], &cn[..e], &sn[..e], &cm[..e], &sm[..e]);
    let mut o = [0.0f64; 6];
    let x0 = cn[0] * cm[0] + sn[0] * sm[0];
    let y0 = sn[0] * cm[0] - cn[0] * sm[0];
    leg_term::<L2>(&mut o, 0.0, z[0], 0.0, z[0], c * z[1], z[1], c2 * z[2], kr[0], jr[0], x0, y0, u);
    if e > 1 {
        let x1 = cn[1] * cm[1] + sn[1] * sm[1];
        let y1 = sn[1] * cm[1] - cn[1] * sm[1];
        leg_term::<L2>(&mut o, 1.0, c * z[1], z[1], z[1], c2 * z[2], c * z[2], c2 * z[3], kr[1], jr[1], x1, y1, u);
    }
    for m in 2..e {
        let x = cn[m] * cm[m] + sn[m] * sm[m];
        let y = sn[m] * cm[m] - cn[m] * sm[m];
        let (z0, z1, z2) = (z[m], z[m + 1], z[m + 2]);
        leg_term::<L2>(&mut o, m as f64, c2 * z0, c * z0, z0, c2 * z1, c * z1, c2 * z2, kr[m], jr[m], x, y, u);
    }
    o
}

#[inline(always)]
fn leg_accumulate<const L2: bool>(acc: &mut [f64; 9], o: &[f64; 6], n: usize, w: f64) {
    let a1 = n as f64 + 1.0;
    acc[0] += w * a1 * o[0];
    acc[1] += w * o[1];
    acc[2] += w * o[2];
    if L2 {
        let a2 = a1 + 1.0;
        acc[3] += w * a1 * a2 * o[0];
        acc[4] += w * o[3];
        acc[5] += w * o[5];
        acc[6] += w * o[4];
        acc[7] += w * a2 * o[2];
        acc[8] += w * a2 * o[1];
    }
}

#[inline(always)]
fn leg_finish<const L2: bool>(acc: &[f64; 9], gm: f64, r: f64, u: f64, c: f64, cl: f64, sl: f64) -> [f64; 9] {
    let b1 = gm / (r * r);
    let g = [b1 * acc[2], b1 * acc[1], -b1 * acc[0]];
    let mut t = [0.0f64; 6];
    if L2 {
        let b2 = b1 / r;
        t = [b2 * (acc[5] - acc[0]), b2 * acc[6], -b2 * acc[7], b2 * (acc[4] - acc[0]), -b2 * acc[8], b2 * acc[3]];
    }
    local_to_cart::<L2>(g, &t, u, c, cl, sl)
}

/// Sectoral step of the Z recursion: Z_nn = f_n c^(1 + min(n-1,2) - min(n,2)) Z_n-1,n-1.
#[inline(always)]
fn zsect(n: usize, c: f64) -> f64 {
    if n <= 2 { 1.0 } else { c }
}

/// LEGENDRE_a: forward column recursion, two rolling rows, sums on the fly.
/// TAG only makes a separate copy of the code for the A/A probe.
pub struct LegA<const TAG: u8> {
    cm: Vec<f64>,
    sm: Vec<f64>,
    p: Vec<f64>,
    p1: Vec<f64>,
}
impl<const TAG: u8> LegA<TAG> {
    fn new(n: usize) -> Self {
        LegA { cm: vec![0.0; n + 1], sm: vec![0.0; n + 1], p: vec![0.0; n + 3], p1: vec![0.0; n + 3] }
    }
}
impl<const TAG: u8> Kernel for LegA<TAG> {
    #[inline(never)]
    fn eval<const L2: bool>(&mut self, md: &Model, x: f64, y: f64, z: f64) -> [f64; 9] {
        let _ = black_box(TAG);
        let (nmax, n1, t) = (md.n, md.n1, &md.t);
        let l = t.l;
        let (r, u, c, cl, sl) = sph_geom(x, y, z);
        let c2 = c * c;
        let q = md.a / r;
        trig(&mut self.cm, &mut self.sm, cl, sl);
        self.p.fill(0.0);
        self.p1.fill(0.0);
        let mut acc = [0.0f64; 9];
        let mut w = 1.0f64;
        for n in 0..=nmax {
            // p held row n-1, p1 row n-2: after the swap p1 = row n-1 and row n overwrites
            // row n-2 in place (columns n+1, n+2 stay zero: sentinels)
            std::mem::swap(&mut self.p, &mut self.p1);
            {
                let (p, p1) = (&mut self.p, &self.p1);
                if n == 0 {
                    p[0] = 1.0;
                } else {
                    let g = &t.g[n * l..n * l + n];
                    let h = &t.h[n * l..n * l + n];
                    let p1n = &p1[..n];
                    let pn = &mut p[..n];
                    for m in 0..n {
                        pn[m] = g[m] * (u * p1n[m]) - h[m] * pn[m];
                    }
                    p[n] = t.f[n] * zsect(n, c) * p1[n - 1];
                }
            }
            let e = n + 1;
            let o = leg_row::<L2>(&self.p, &t.k[n * l..], &t.j[n * l..], &md.c[n * n1..], &md.s[n * n1..],
                                  &self.cm, &self.sm, e, u, c, c2);
            leg_accumulate::<L2>(&mut acc, &o, n, w);
            w *= q;
        }
        leg_finish::<L2>(&acc, md.gm, r, u, c, cl, sl)
    }
}

/// LEGENDRE_b: the same mathematics organised as store-then-sum (the paper's Matlab + DLL
/// programme): the whole triangle of Z_nm is built first (preallocated), then summed.
pub struct LegB {
    cm: Vec<f64>,
    sm: Vec<f64>,
    tri: Vec<f64>, // row n at off(n) = n (n + 5) / 2, length n + 3 (two zero sentinels)
}
#[inline(always)]
fn off(n: usize) -> usize {
    n * (n + 5) / 2
}
impl LegB {
    fn new(n: usize) -> Self {
        LegB { cm: vec![0.0; n + 1], sm: vec![0.0; n + 1], tri: vec![0.0; off(n + 1)] }
    }
}
impl Kernel for LegB {
    #[inline(never)]
    fn eval<const L2: bool>(&mut self, md: &Model, x: f64, y: f64, z: f64) -> [f64; 9] {
        let (nmax, n1, t) = (md.n, md.n1, &md.t);
        let l = t.l;
        let (r, u, c, cl, sl) = sph_geom(x, y, z);
        let c2 = c * c;
        let q = md.a / r;
        trig(&mut self.cm, &mut self.sm, cl, sl);
        let tri = &mut self.tri;
        tri[0] = 1.0;
        tri[1] = 0.0;
        tri[2] = 0.0;
        for n in 1..=nmax {
            let (o, o1) = (off(n), off(n - 1));
            let (lo, hi) = tri.split_at_mut(o);
            let row = &mut hi[..n + 3];
            let r1 = &lo[o1..o1 + n + 2];
            let g = &t.g[n * l..n * l + n];
            if n == 1 {
                row[0] = g[0] * (u * r1[0]);
            } else {
                let o2 = off(n - 2);
                let r2 = &lo[o2..o2 + n + 1];
                let h = &t.h[n * l..n * l + n];
                for m in 0..n {
                    row[m] = g[m] * (u * r1[m]) - h[m] * r2[m];
                }
            }
            row[n] = t.f[n] * zsect(n, c) * r1[n - 1];
            row[n + 1] = 0.0;
            row[n + 2] = 0.0;
        }
        let mut acc = [0.0f64; 9];
        let mut w = 1.0f64;
        for n in 0..=nmax {
            let e = n + 1;
            let o = off(n);
            let s = leg_row::<L2>(&self.tri[o..o + n + 3], &t.k[n * l..], &t.j[n * l..], &md.c[n * n1..],
                                  &md.s[n * n1..], &self.cm, &self.sm, e, u, c, c2);
            leg_accumulate::<L2>(&mut acc, &s, n, w);
            w *= q;
        }
        leg_finish::<L2>(&acc, md.gm, r, u, c, cl, sl)
    }
}

// ------------------------------------------------------------------ PINES
pub struct Pines {
    rm: Vec<f64>, // (s + i t)^m at index m + 2
    im: Vec<f64>,
    rows: [Vec<f64>; 4],
}
impl Pines {
    fn new(n: usize) -> Self {
        Pines { rm: vec![0.0; n + 3], im: vec![0.0; n + 3], rows: std::array::from_fn(|_| vec![0.0; n + 5]) }
    }
}
#[inline(always)]
fn hrow4(rows: &mut [Vec<f64>; 4], t: &Tables, n: usize, u: f64) {
    let l = t.l;
    let mut an = std::mem::take(&mut rows[n % 4]);
    {
        let a1 = &rows[(n + 3) % 4][..n];
        let a2 = &rows[(n + 2) % 4][..n];
        let g = &t.g[n * l..n * l + n];
        let h = &t.h[n * l..n * l + n];
        let anm = &mut an[..n];
        for m in 0..n {
            anm[m] = g[m] * (u * a1[m]) - h[m] * a2[m];
        }
        an[n] = t.f[n] * a1[n - 1];
    }
    rows[n % 4] = an;
}
#[inline(always)]
fn powers(rm: &mut [f64], im: &mut [f64], nmax: usize, s: f64, tt: f64) {
    rm[0] = 0.0;
    im[0] = 0.0;
    rm[1] = 0.0;
    im[1] = 0.0;
    rm[2] = 1.0;
    im[2] = 0.0;
    for m in 1..=nmax {
        rm[m + 2] = s * rm[m + 1] - tt * im[m + 1];
        im[m + 2] = s * im[m + 1] + tt * rm[m + 1];
    }
}
impl Kernel for Pines {
    #[inline(never)]
    fn eval<const L2: bool>(&mut self, md: &Model, x: f64, y: f64, z: f64) -> [f64; 9] {
        let (nmax, n1, t) = (md.n, md.n1, &md.t);
        let l = t.l;
        let r = (x * x + y * y + z * z).sqrt();
        let (s, tt, u) = (x / r, y / r, z / r);
        let q = md.a / r;
        powers(&mut self.rm, &mut self.im, nmax, s, tt);
        for row in self.rows.iter_mut() {
            row.fill(0.0);
        }
        self.rows[0][0] = 1.0;
        hrow4(&mut self.rows, t, 1, u);
        if L2 {
            hrow4(&mut self.rows, t, 2, u);
        }
        let mut acc = [0.0f64; 13];
        let mut w = 1.0f64;
        for n in 0..=nmax {
            let e = n + 1;
            {
                let a0 = &self.rows[n % 4][..e + 2];
                let b1 = &self.rows[(n + 1) % 4][..e + 2];
                let b2 = &self.rows[(n + 2) % 4][..e + 2];
                let (kr, jr) = (&t.k[n * l..n * l + e], &t.j[n * l..n * l + e]);
                let (r11, r12, r22) = (&t.r11[n * l..n * l + e], &t.r12[n * l..n * l + e], &t.r22[n * l..n * l + e]);
                let (cn, sn) = (&md.c[n * n1..n * n1 + e], &md.s[n * n1..n * n1 + e]);
                let (rm, im) = (&self.rm[..e + 2], &self.im[..e + 2]);
                let mut o = [0.0f64; 13];
                for m in 0..e {
                    let fm = m as f64;
                    let (c_, s_) = (cn[m], sn[m]);
                    let d0 = c_ * rm[m + 2] + s_ * im[m + 2];
                    let e1 = fm * (c_ * rm[m + 1] + s_ * im[m + 1]);
                    let f1 = fm * (s_ * rm[m + 1] - c_ * im[m + 1]);
                    let v0 = a0[m];
                    let v1 = kr[m] * a0[m + 1];
                    let v2 = -r11[m] * b1[m + 1];
                    o[0] += v1 * d0;
                    o[1] += v2 * d0;
                    o[2] += v0 * e1;
                    o[3] += v0 * f1;
                    if L2 {
                        let v3 = jr[m] * a0[m + 2];
                        let v4 = -r12[m] * b1[m + 2];
                        let v5 = r22[m] * b2[m + 2];
                        o[4] += v3 * d0;
                        o[5] += v4 * d0;
                        o[6] += v5 * d0;
                        o[7] += v1 * e1;
                        o[8] += v1 * f1;
                        o[9] += v2 * e1;
                        o[10] += v2 * f1;
                        let mm1 = fm * (fm - 1.0);
                        o[11] += v0 * mm1 * (c_ * rm[m] + s_ * im[m]);
                        o[12] += v0 * mm1 * (s_ * rm[m] - c_ * im[m]);
                    }
                }
                let no = if L2 { 13 } else { 4 };
                for i in 0..no {
                    acc[i] += w * o[i];
                }
            }
            w *= q;
            if !L2 && n + 2 <= nmax + 1 {
                hrow4(&mut self.rows, t, n + 2, u);
            }
            if L2 && n + 3 <= nmax + 2 {
                hrow4(&mut self.rows, t, n + 3, u);
            }
        }
        let r1 = md.gm / (r * r);
        let (a3, a4, a1, a2) = (acc[0] * r1, acc[1] * r1, acc[2] * r1, acc[3] * r1);
        let mut out = [a1 + s * a4, a2 + tt * a4, a3 + u * a4, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        if L2 {
            let r2 = r1 / r;
            let (a33, a34, a44) = (acc[4] * r2, acc[5] * r2, acc[6] * r2);
            let (a13, a23, a14, a24) = (acc[7] * r2, acc[8] * r2, acc[9] * r2, acc[10] * r2);
            let (a11, a12) = (acc[11] * r2, acc[12] * r2);
            let xh = [s, tt, u];
            let ai4 = [a14, a24, a34];
            let aij = [[a11, a12, a13], [a12, -a11, a23], [a13, a23, a33]];
            let a4r = a4 / r;
            let pairs = [(0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2)];
            for (col, &(i, j)) in pairs.iter().enumerate() {
                let mut vij = aij[i][j] + xh[j] * ai4[i] + xh[i] * ai4[j] + xh[i] * xh[j] * a44;
                if i == j {
                    vij += a4r;
                }
                out[3 + col] = vij;
            }
        }
        out
    }
}

// ------------------------------------------------------------------ GOTTLIEB
pub struct Gottlieb {
    rm: Vec<f64>,
    im: Vec<f64>,
    rows: [Vec<f64>; 3],
}
impl Gottlieb {
    fn new(n: usize) -> Self {
        Gottlieb { rm: vec![0.0; n + 3], im: vec![0.0; n + 3], rows: std::array::from_fn(|_| vec![0.0; n + 3]) }
    }
}
impl Kernel for Gottlieb {
    #[inline(never)]
    fn eval<const L2: bool>(&mut self, md: &Model, x: f64, y: f64, z: f64) -> [f64; 9] {
        let (nmax, n1, t) = (md.n, md.n1, &md.t);
        let l = t.l;
        let r = (x * x + y * y + z * z).sqrt();
        let (s, tt, u) = (x / r, y / r, z / r);
        let q = md.a / r;
        powers(&mut self.rm, &mut self.im, nmax, s, tt);
        for row in self.rows.iter_mut() {
            row.fill(0.0);
        }
        self.rows[0][0] = 1.0;
        // L1: gam j k h ; L2: + l s t n o p q r M
        let mut acc = [0.0f64; 13];
        let mut w = 1.0f64;
        for n in 0..=nmax {
            if n >= 1 {
                let mut an = std::mem::take(&mut self.rows[n % 3]);
                {
                    let a1 = &self.rows[(n + 2) % 3][..n];
                    let a2 = &self.rows[(n + 1) % 3][..n];
                    let g = &t.g[n * l..n * l + n];
                    let h = &t.h[n * l..n * l + n];
                    let anm = &mut an[..n];
                    for m in 0..n {
                        anm[m] = g[m] * (u * a1[m]) - h[m] * a2[m];
                    }
                    an[n] = t.f[n] * a1[n - 1];
                }
                self.rows[n % 3] = an;
            }
            let e = n + 1;
            let a0 = &self.rows[n % 3][..e + 2];
            let (kr, jr) = (&t.k[n * l..n * l + e], &t.j[n * l..n * l + e]);
            let (cn, sn) = (&md.c[n * n1..n * n1 + e], &md.s[n * n1..n * n1 + e]);
            let (rm, im) = (&self.rm[..e + 2], &self.im[..e + 2]);
            let nf = n as f64;
            let mut o = [0.0f64; 13];
            for m in 0..e {
                let fm = m as f64;
                let (c_, s_) = (cn[m], sn[m]);
                let b0 = c_ * rm[m + 2] + s_ * im[m + 2];
                let b1 = c_ * rm[m + 1] + s_ * im[m + 1];
                let a1 = s_ * rm[m + 1] - c_ * im[m + 1];
                let np1 = nf + fm + 1.0;
                let p0 = a0[m];
                let pb = p0 * b0;
                let mp = fm * p0;
                let (mb, ma) = (mp * b1, mp * a1);
                let zg = kr[m] * a0[m + 1];
                let zb = zg * b0;
                o[0] += np1 * pb;
                o[1] += mb;
                o[2] += ma;
                o[3] += zb;
                if L2 {
                    o[4] += np1 * (np1 + 1.0) * pb;
                    o[5] += np1 * mb;
                    o[6] += np1 * ma;
                    let mm1 = (fm - 1.0) * mp;
                    o[7] += mm1 * (c_ * rm[m] + s_ * im[m]);
                    o[8] += mm1 * (c_ * im[m] - s_ * rm[m]);
                    o[9] += np1 * zb;
                    o[10] += fm * zg * b1;
                    o[11] += fm * zg * a1;
                    o[12] += jr[m] * a0[m + 2] * b0;
                }
            }
            let no = if L2 { 13 } else { 4 };
            for i in 0..no {
                acc[i] += w * o[i];
            }
            w *= q;
        }
        let (sgam, sj, sk, sh) = (acc[0], acc[1], acc[2], acc[3]);
        let lam_ = sgam + u * sh;
        let m2 = md.gm / (r * r);
        let mut out = [m2 * (sj - lam_ * s), m2 * (sk - lam_ * tt), m2 * (sh - lam_ * u), 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        if L2 {
            let (sl_, ss, st, sn_, so, spp, sq, sr, sm) =
                (acc[4], acc[5], acc[6], acc[7], acc[8], acc[9], acc[10], acc[11], acc[12]);
            let d1 = ss + u * sq;
            let d2 = st + u * sr;
            let gg = -(spp + sh + u * sm);
            let ff = sl_ + lam_ + u * (spp + sh - gg);
            let m3 = m2 / r;
            out[3] = m3 * ((ff * s - 2.0 * d1) * s - lam_ + sn_);
            out[4] = m3 * ((ff * tt - d2) * s - d1 * tt - so);
            out[5] = m3 * ((ff * s - d1) * u + gg * s + sq);
            out[6] = m3 * ((ff * tt - 2.0 * d2) * tt - lam_ - sn_);
            out[7] = m3 * ((ff * tt - d2) * u + gg * tt + sr);
            out[8] = m3 * ((ff * u + 2.0 * gg) * u - lam_ + sm);
        }
        out
    }
}

// ------------------------------------------------------------------ METRIS
/// Solid spherical harmonics; the derivative operators act on the coefficients at setup:
/// level 1 stores 3 complex coefficient sets (gradient, degrees <= N+1), level 2 stores 9
/// (gradient + tensor, degrees <= N+2), packed triangular, interleaved per (n, m).
pub struct Metris {
    st: usize,
    nd: usize,
    er: Vec<f64>, // [(n (n+1) / 2 + m) * st + q]
    ei: Vec<f64>,
    yr: [Vec<f64>; 3],
    yi: [Vec<f64>; 3],
}
fn metris_pack(er: &mut [f64], ei: &mut [f64], src: &(Vec<f64>, Vec<f64>), qi: usize, st: usize, nd: usize, l: usize) {
    for n in 0..=nd {
        for m in 0..=n {
            let i = (n * (n + 1) / 2 + m) * st + qi;
            er[i] = src.0[n * l + m];
            ei[i] = src.1[n * l + m];
        }
    }
}
impl Metris {
    fn new(md: &Model, l2: bool) -> Self {
        let (n, n1, a, t) = (md.n, md.n1, md.a, &md.t);
        let l = n + 3;
        let st = if l2 { 9 } else { 3 };
        let nd = n + 1 + l2 as usize;
        let mut kr = vec![0.0f64; l * l];
        let mut ki = vec![0.0f64; l * l];
        let sc = md.gm / a;
        for i in 0..n1 {
            for m in 0..=i {
                kr[i * l + m] = sc * md.c[i * n1 + m];
                ki[i * l + m] = -sc * md.s[i * n1 + m];
            }
        }
        let ntri = (nd + 1) * (nd + 2) / 2;
        let mut er = vec![0.0f64; ntri * st];
        let mut ei = vec![0.0f64; ntri * st];
        let ex = metris_op(&kr, &ki, 'x', t, l, a);
        let ey = metris_op(&kr, &ki, 'y', t, l, a);
        let ez = metris_op(&kr, &ki, 'z', t, l, a);
        drop((kr, ki));
        metris_pack(&mut er, &mut ei, &ex, 0, st, nd, l);
        metris_pack(&mut er, &mut ei, &ey, 1, st, nd, l);
        metris_pack(&mut er, &mut ei, &ez, 2, st, nd, l);
        if l2 {
            let ops: [(&(Vec<f64>, Vec<f64>), char); 6] =
                [(&ex, 'x'), (&ex, 'y'), (&ex, 'z'), (&ey, 'y'), (&ey, 'z'), (&ez, 'z')];
            for (k, (src, kind)) in ops.iter().enumerate() {
                let d = metris_op(&src.0, &src.1, *kind, t, l, a);
                metris_pack(&mut er, &mut ei, &d, 3 + k, st, nd, l);
            }
        }
        Metris {
            st, nd, er, ei,
            yr: std::array::from_fn(|_| vec![0.0; nd + 1]),
            yi: std::array::from_fn(|_| vec![0.0; nd + 1]),
        }
    }
}
impl Kernel for Metris {
    #[inline(never)]
    fn eval<const L2: bool>(&mut self, md: &Model, x: f64, y: f64, z: f64) -> [f64; 9] {
        let st: usize = if L2 { 9 } else { 3 };
        assert!(st == self.st, "METRIS built for the other level");
        let t = &md.t;
        let lt = t.l;
        let nd = self.nd;
        let a = md.a;
        let (xp, yp, zp) = (x / a, y / a, z / a);
        let rp2 = xp * xp + yp * yp + zp * zp;
        let rp = rp2.sqrt();
        let r2i = 1.0 / rp2;
        let (zr, zeta_r, zeta_i) = (zp * r2i, xp * r2i, yp * r2i);
        for k in 0..3 {
            self.yr[k].fill(0.0);
            self.yi[k].fill(0.0);
        }
        let mut f = [0.0f64; 9];
        for n in 0..=nd {
            let mut r0 = std::mem::take(&mut self.yr[n % 3]);
            let mut i0 = std::mem::take(&mut self.yi[n % 3]);
            if n == 0 {
                r0[0] = 1.0 / rp;
            } else {
                let r1 = &self.yr[(n + 2) % 3][..n];
                let i1 = &self.yi[(n + 2) % 3][..n];
                let r2 = &self.yr[(n + 1) % 3][..n];
                let i2 = &self.yi[(n + 1) % 3][..n];
                let g = &t.g[n * lt..n * lt + n];
                let h = &t.h[n * lt..n * lt + n];
                {
                    let (r0n, i0n) = (&mut r0[..n], &mut i0[..n]);
                    for m in 0..n {
                        r0n[m] = g[m] * (zr * r1[m]) - h[m] * (r2i * r2[m]);
                        i0n[m] = g[m] * (zr * i1[m]) - h[m] * (r2i * i2[m]);
                    }
                }
                let (ar, ai) = (r1[n - 1], i1[n - 1]);
                r0[n] = t.f[n] * (zeta_r * ar - zeta_i * ai);
                i0[n] = t.f[n] * (zeta_r * ai + zeta_i * ar);
            }
            let e = n + 1;
            let base = n * (n + 1) / 2 * st;
            let er = &self.er[base..base + e * st];
            let ei = &self.ei[base..base + e * st];
            let mut s = [0.0f64; 9];
            for m in 0..e {
                let (yr_, yi_) = (r0[m], i0[m]);
                let cr = &er[m * st..m * st + st];
                let ci = &ei[m * st..m * st + st];
                for qi in 0..st {
                    s[qi] += yr_ * cr[qi] - yi_ * ci[qi];
                }
            }
            for qi in 0..st {
                f[qi] += s[qi];
            }
            self.yr[n % 3] = r0;
            self.yi[n % 3] = i0;
        }
        f
    }
}

// ------------------------------------------------------------------ CLENSHAW (scattered)
/// Clenshaw over the degree for one order m, cosine and sine families as two lanes.
/// Lanes: 0 y, 1 y_r (weights n+1), 2 y_u; level 2 adds 3 y_rr ((n+1)(n+2)), 4 y_uu, 5 y_ur.
#[inline(always)]
fn clen_column<const L2: bool>(g: &[f64], h: &[f64], cc: &[f64], cs: &[f64], m: usize, uq: f64, q2: f64,
                               q: f64) -> [[f64; 2]; 6] {
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
        let cn = [cc[i], cs[i]];
        let mut y = [[0.0f64; 2]; 6];
        for p in 0..2 {
            y[0][p] = aa * a[0][p] + (cn[p] - bb * b[0][p]);
            y[1][p] = aa * a[1][p] + (nn1 * cn[p] - bb * b[1][p]);
            y[2][p] = aa * a[2][p] + (qg * a[0][p] - bb * b[2][p]);
            if L2 {
                let nn2 = nn1 * (nn1 + 1.0);
                y[3][p] = aa * a[3][p] + (nn2 * cn[p] - bb * b[3][p]);
                y[4][p] = aa * a[4][p] + ((2.0 * qg) * a[2][p] - bb * b[4][p]);
                y[5][p] = aa * a[5][p] + (qg * a[1][p] - bb * b[5][p]);
            }
        }
        b = a;
        a = y;
    }
    a
}

/// CLENSHAW for scattered points: per order m, Clenshaw over n, seeded by
/// zs = (a/r)^m Z_mm; then the pole-safe local-frame terms (as in LEGENDRE_a).
pub struct Clen;
impl Kernel for Clen {
    #[inline(never)]
    fn eval<const L2: bool>(&mut self, md: &Model, x: f64, y: f64, z: f64) -> [f64; 9] {
        let (nmax, n1, t) = (md.n, md.n1, &md.t);
        let l = t.l;
        let (r, u, c, cl, sl) = sph_geom(x, y, z);
        let c2 = c * c;
        let q = md.a / r;
        let (uq, q2) = (u * q, q * q);
        let mut zs = 1.0f64;
        let (mut cmv, mut smv) = (1.0f64, 0.0f64);
        let mut acc = [0.0f64; 9];
        for m in 0..=nmax {
            if m >= 1 {
                zs *= t.f[m] * q * zsect(m, c);
                let (a_, b_) = (cmv, smv);
                cmv = a_ * cl - b_ * sl;
                smv = b_ * cl + a_ * sl;
            }
            let yv = clen_column::<L2>(&t.gt[m * l + m..m * l + nmax + 3], &t.ht[m * l + m..m * l + nmax + 3],
                                       &md.ct[m * n1 + m..m * n1 + n1], &md.st[m * n1 + m..m * n1 + n1],
                                       m, uq, q2, q);
            let fm = m as f64;
            let (cpm, cqm) = match m {
                0 => (1.0, 0.0),
                1 => (c, 1.0),
                _ => (c2, c),
            };
            let (sp, sq, sr) = (zs * cpm, zs * cqm, zs);
            let xx = |i: usize| yv[i][0] * cmv + yv[i][1] * smv;
            let yy = |i: usize| yv[i][1] * cmv - yv[i][0] * smv;
            let (x0, xr, xu) = (xx(0), xx(1), xx(2));
            let y0 = yy(0);
            acc[0] += sp * xr;
            acc[1] += c * sp * xu - fm * u * sq * x0;
            acc[2] += fm * sq * y0;
            if L2 {
                let (xrr, xuu, xur) = (xx(3), xx(4), xx(5));
                let (yr1, yu) = (yy(1), yy(2));
                let mm1 = fm * (fm - 1.0);
                acc[3] += sp * xrr;
                acc[4] += c * sp * (xur + xu) - fm * u * sq * (xr + x0);
                acc[5] += fm * sq * (yr1 + y0);
                acc[6] += (mm1 * sr - fm * fm * sp) * x0 - (2.0 * fm + 1.0) * u * sp * xu + c2 * sp * xuu - sp * xr;
                acc[7] += (-(mm1 * sr) - fm * sp) * x0 - u * sp * xu - sp * xr;
                acc[8] += fm * sp * yu - mm1 * u * sr * y0;
            }
        }
        let b1 = md.gm / (r * r);
        let g = [b1 * acc[2], b1 * acc[1], -b1 * acc[0]];
        let mut tl = [0.0f64; 6];
        if L2 {
            let b2 = b1 / r;
            tl = [b2 * acc[7], b2 * acc[8], -b2 * acc[5], b2 * acc[6], -b2 * acc[4], b2 * acc[3]];
        }
        local_to_cart::<L2>(g, &tl, u, c, cl, sl)
    }
}

// ------------------------------------------------------------------ LEAR (Eckman-normalised)
// Lear (1986) in the fully normalised form of Eckman, Brown & Adamo (2016), NASA/TP-2016-218604,
// app. C.7 (learnorm.m). Variables: s = sin(phi), c = cos(phi);
//   pn  = Pbar_n0,  ppn  = dPbar_n0/ds;
//   pnm = Pbar_nm / c (m >= 1),  ppnm = dPbar_nm/dphi   -- no 1/c anywhere, pole-safe.
// Recursions (learnorm.m, constants of each line folded into setup tables as for every method):
//   pn(n)   = ((2n-1) s norm1 pn(n-1) - (n-1) norm2 pn(n-2)) / n
//   ppn(n)  = norm1 (s ppn(n-1) + n pn(n-1))
//   pnm(n,n)= (2n-1) norm11 c pnm(n-1,n-1),   ppnm(n,n) = -n s pnm(n,n)
//   e3      = norm1m pnm(n-1,m)
//   pnm(n,m)= ((2n-1) s e3 - (n+m-1) norm2m pnm(n-2,m)) / (n-m),   ppnm(n,m) = -n s pnm + (n+m) e3
// Level 1 (Lear): a_r = mu/r^2 [-sum (n+1) rb C_n0 pn - c sum (n+1) rb sum_m e4 pnm],
//   a_e = mu/r^2 sum rb sum_m m e5 pnm,  a_n = mu/r^2 [c sum rb C_n0 ppn + sum rb sum_m e4 ppnm],
//   e4 = C cos(m lam) + S sin(m lam),  e5 = S cos(m lam) - C sin(m lam).
//
// Level 2 is our extension (Lear gives the acceleration only). Local-frame terms per (n, m), times
// mu rb / r^3, with P = Pbar_nm and ' = d/dphi:
//   ee = -m^2 P/c^2 - (s/c) P' - (n+1) P     nn = P'' - (n+1) P       rr = (n+1)(n+2) P
//   en = m e5 (P'/c + s P/c^2)               er = -(n+2) m e5 P/c     nr = -(n+2) P'
// er, nr, rr and the (n+1) P terms reuse the level-1 sums. The 1/c and 1/c^2 of ee, en, nn are
// removed in Lear's own style, per kind of column (no division by c, no Legendre ODE):
//   m = 0 : one more s-derivative, as Lear's ppn:  pppn(n) = norm1 (s pppn(n-1) + (n+1) ppn(n-1));
//           ee + (n+1)P = -s ppn,  P'' = c^2 pppn - s ppn.
//   m = 1 : pnm = Q(s) is a polynomial; Q_s, Q_ss by differentiating its column recursion;
//           ee + (n+1)P = -c (Q + s Q_s),  en = c Q_s e5,  P'' = c (c^2 Q_ss - 3 s Q_s - Q).
//   m >= 2: q = Pbar_nm / c^2 by the same column recursion as pnm, seed q(n,n) = (2n-1) norm11
//           pnm(n-1,n-1). Lear's derivative identity divided by c^2 gives G = P'/c, its
//           phi-derivative divided by c gives P'':
//           G = -n s q + (n+m) norm1m q(n-1,m),   P'' = (1-n) s G - n c^2 q + (n+m) norm1m G(n-1,m);
//           ee + (n+1)P = -m^2 q - s G,  en = m e5 (G + s q).
pub struct Lear {
    l: usize,
    cm: Vec<f64>,
    sm: Vec<f64>,
    p: Vec<f64>,
    p1: Vec<f64>,
    q: Vec<f64>,   // level 2: q = Pbar/c^2, row n
    q1: Vec<f64>,  // level 2: row n-1
    gb: Vec<f64>,  // level 2: G = P'/c of row n-1, overwritten by row n
    za: Vec<f64>,  // (2n-1) norm1 / n
    zb: Vec<f64>,  // (n-1) norm2 / n
    zn1: Vec<f64>, // norm1
    sec: Vec<f64>, // (2n-1) norm11 (n >= 2)
    n1m: Vec<f64>, // [n * l + m] norm1m
    ta: Vec<f64>,  // (2n-1) / (n-m)
    tb: Vec<f64>,  // (n+m-1) norm2m / (n-m)
}
impl Lear {
    fn new(n: usize) -> Self {
        let l = n + 2;
        let mut k = Lear {
            l, cm: vec![0.0; n + 1], sm: vec![0.0; n + 1], p: vec![0.0; n + 2], p1: vec![0.0; n + 2],
            q: vec![0.0; n + 2], q1: vec![0.0; n + 2], gb: vec![0.0; n + 2],
            za: vec![0.0; l], zb: vec![0.0; l], zn1: vec![0.0; l], sec: vec![0.0; l],
            n1m: vec![0.0; l * l], ta: vec![0.0; l * l], tb: vec![0.0; l * l],
        };
        for ni in 1..=n {
            let fnn = ni as f64;
            let norm1 = ((2.0 * fnn + 1.0) / (2.0 * fnn - 1.0)).sqrt();
            k.zn1[ni] = norm1;
            k.za[ni] = (2.0 * fnn - 1.0) * norm1 / fnn;
            if ni >= 2 {
                let norm2 = ((2.0 * fnn + 1.0) / (2.0 * fnn - 3.0)).sqrt();
                k.zb[ni] = (fnn - 1.0) * norm2 / fnn;
                k.sec[ni] = (2.0 * fnn - 1.0) * ((2.0 * fnn + 1.0) / (2.0 * fnn)).sqrt() / (2.0 * fnn - 1.0);
            }
            for mi in 1..ni {
                let fm = mi as f64;
                let i = ni * l + mi;
                k.n1m[i] = ((fnn - fm) * (2.0 * fnn + 1.0) / ((fnn + fm) * (2.0 * fnn - 1.0))).sqrt();
                k.ta[i] = (2.0 * fnn - 1.0) / (fnn - fm);
                if ni >= 3 {
                    let norm2m = ((fnn - fm) * (fnn - fm - 1.0) * (2.0 * fnn + 1.0)
                        / ((fnn + fm) * (fnn + fm - 1.0) * (2.0 * fnn - 3.0))).sqrt();
                    k.tb[i] = (fnn + fm - 1.0) * norm2m / (fnn - fm);
                }
            }
        }
        k
    }
}
impl Kernel for Lear {
    #[inline(never)]
    fn eval<const L2: bool>(&mut self, md: &Model, x: f64, y: f64, z: f64) -> [f64; 9] {
        let (nmax, n1, l) = (md.n, md.n1, self.l);
        let (r, s, c, cl, sl) = sph_geom(x, y, z);
        let qr = md.a / r;
        let c2 = c * c;
        trig(&mut self.cm, &mut self.sm, cl, sl);
        self.p.fill(0.0);
        self.p1.fill(0.0);
        if L2 {
            self.q.fill(0.0);
            self.q1.fill(0.0);
            self.gb.fill(0.0);
        }
        let (mut pn, mut pn1, mut ppn) = (1.0f64, 0.0f64, 0.0f64); // pn(n), pn(n-1), ppn(n)
        let mut pppn = 0.0f64; // level 2: d2 Pbar_n0 / ds2
        // level 2, m = 1 column: Q_s and Q_ss at rows n-1 (a) and n-2 (b)
        let (mut d1a, mut d1b, mut dd1a, mut dd1b) = (0.0f64, 0.0f64, 0.0f64, 0.0f64);
        let (mut zr, mut zn, mut t1, mut te, mut t3) = (0.0f64, 0.0f64, 0.0f64, 0.0f64, 0.0f64);
        let (mut aee, mut ann, mut aen, mut aer, mut anr, mut arr) = (0.0f64, 0.0f64, 0.0f64, 0.0f64, 0.0f64, 0.0f64);
        let mut w = 1.0f64;
        for n in 0..=nmax {
            let fnn = n as f64;
            let (mut e1, mut e2, mut e3s) = (0.0f64, 0.0f64, 0.0f64);
            let (mut xee, mut xnn, mut xen) = (0.0f64, 0.0f64, 0.0f64); // level 2, m >= 1
            if n >= 1 {
                // zonal
                if L2 {
                    pppn = self.zn1[n] * (s * pppn + (fnn + 1.0) * ppn);
                }
                ppn = self.zn1[n] * (s * ppn + fnn * pn);
                let pnew = self.za[n] * (s * pn) - self.zb[n] * pn1;
                pn1 = pn;
                pn = pnew;
                // tesseral + sectorial: p1 = row n-1, row n overwrites row n-2 in place
                std::mem::swap(&mut self.p, &mut self.p1);
                if L2 {
                    std::mem::swap(&mut self.q, &mut self.q1);
                }
                let (p, p1, qv, qv1, gb) = (&mut self.p, &self.p1, &mut self.q, &self.q1, &mut self.gb);
                let (cn, sn) = (&md.c[n * n1..n * n1 + n + 1], &md.s[n * n1..n * n1 + n + 1]);
                let (cm, sm) = (&self.cm[..n + 1], &self.sm[..n + 1]);
                let (n1m, ta, tb) = (&self.n1m[n * l..n * l + n], &self.ta[n * l..n * l + n], &self.tb[n * l..n * l + n]);
                let ns = -fnn * s;
                let (a_s, nc2) = ((1.0 - fnn) * s, fnn * c2);
                if n >= 2 {
                    // m = 1
                    let e3 = n1m[1] * p1[1];
                    let v = ta[1] * (s * e3) - tb[1] * p[1];
                    let pp = ns * v + (fnn + 1.0) * e3;
                    p[1] = v;
                    let e4 = cn[1] * cm[1] + sn[1] * sm[1];
                    let e5 = sn[1] * cm[1] - cn[1] * sm[1];
                    e1 += e4 * v;
                    e2 += e5 * v;
                    e3s += e4 * pp;
                    if L2 {
                        let d = ta[1] * (e3 + s * (n1m[1] * d1a)) - tb[1] * d1b;
                        let dd = ta[1] * (2.0 * (n1m[1] * d1a) + s * (n1m[1] * dd1a)) - tb[1] * dd1b;
                        (d1b, d1a, dd1b, dd1a) = (d1a, d, dd1a, dd);
                        xee -= e4 * (c * (v + s * d));
                        xnn += e4 * (c * (c2 * dd - 3.0 * s * d - v));
                        xen += e5 * (c * d);
                    }
                }
                for m in 2..n {
                    let fm = m as f64;
                    let e3 = n1m[m] * p1[m];
                    let v = ta[m] * (s * e3) - tb[m] * p[m];
                    let pp = ns * v + (fnn + fm) * e3;
                    p[m] = v;
                    let e4 = cn[m] * cm[m] + sn[m] * sm[m];
                    let e5 = sn[m] * cm[m] - cn[m] * sm[m];
                    e1 += e4 * v;
                    e2 += fm * e5 * v;
                    e3s += e4 * pp;
                    if L2 {
                        let e3q = n1m[m] * qv1[m];
                        let vq = ta[m] * (s * e3q) - tb[m] * qv[m];
                        qv[m] = vq;
                        let g = ns * vq + (fnn + fm) * e3q;
                        let gp = gb[m];
                        gb[m] = g;
                        xee -= e4 * (fm * fm * vq + s * g);
                        xnn += e4 * (a_s * g - nc2 * vq + (fnn + fm) * (n1m[m] * gp));
                        xen += fm * e5 * (g + s * vq);
                    }
                }
                let v = if n == 1 { 3.0f64.sqrt() } else { self.sec[n] * c * p1[n - 1] };
                p[n] = v;
                let e4 = cn[n] * cm[n] + sn[n] * sm[n];
                let e5 = sn[n] * cm[n] - cn[n] * sm[n];
                e1 += e4 * v;
                e2 += fnn * e5 * v;
                e3s += e4 * (ns * v);
                if L2 {
                    if n == 1 {
                        // Q = sqrt(3), Q_s = Q_ss = 0
                        xee -= e4 * (c * v);
                        xnn -= e4 * (c * v);
                    } else {
                        let vq = self.sec[n] * p1[n - 1];
                        qv[n] = vq;
                        let g = ns * vq; // G(n-1, n) = 0
                        gb[n] = g;
                        xee -= e4 * (fnn * fnn * vq + s * g);
                        xnn += e4 * (a_s * g - nc2 * vq);
                        xen += fnn * e5 * (g + s * vq);
                    }
                }
            }
            let c0 = md.c[n * n1];
            let a1 = fnn + 1.0;
            zr += w * (a1 * c0 * pn);
            zn += w * (c0 * ppn);
            t1 += w * (a1 * e1);
            te += w * e2;
            t3 += w * e3s;
            if L2 {
                let (a2, sp) = (fnn + 2.0, s * ppn);
                aee += w * (xee - c0 * sp);
                ann += w * (xnn + c0 * (c2 * pppn - sp));
                aen += w * xen;
                aer += w * (a2 * e2);
                anr += w * (a2 * (c * (c0 * ppn) + e3s));
                arr += w * (a1 * a2 * (c0 * pn + c * e1));
            }
            w *= qr;
        }
        let b1 = md.gm / (r * r);
        let ar = -b1 * (zr + c * t1);
        let ae = b1 * te;
        let an = b1 * (c * zn + t3);
        let e5 = ar * c - an * s;
        let (gx, gy, gz) = (e5 * cl - ae * sl, e5 * sl + ae * cl, ar * s + an * c);
        if !L2 {
            return [gx, gy, gz, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        }
        let b2 = b1 / r;
        let p0 = zr + c * t1; // sum rb (n+1) P e4
        let t = [b2 * (aee - p0), b2 * aen, -b2 * aer, b2 * (ann - p0), -b2 * anr, b2 * arr];
        let mut out = local_to_cart::<true>([ae, an, ar], &t, s, c, cl, sl);
        (out[0], out[1], out[2]) = (gx, gy, gz); // level-1 rows exactly as at level 1
        out
    }
}

// ------------------------------------------------------------------ driver
pub const ALL: [&str; 8] = [
    "LEGENDRE_a", "LEGENDRE_a#2", "LEGENDRE_b", "PINES", "GOTTLIEB", "METRIS", "CLENSHAW", "LEAR",
];

enum Kern {
    A(LegA<0>),
    A2(LegA<1>),
    B(LegB),
    P(Pines),
    G(Gottlieb),
    M(Metris),
    C(Clen),
    L(Lear),
}

fn make(name: &str, md: &Model, l2: bool) -> Kern {
    match name {
        "LEGENDRE_a" => Kern::A(LegA::new(md.n)),
        "LEGENDRE_a#2" => Kern::A2(LegA::new(md.n)),
        "LEGENDRE_b" => Kern::B(LegB::new(md.n)),
        "PINES" => Kern::P(Pines::new(md.n)),
        "GOTTLIEB" => Kern::G(Gottlieb::new(md.n)),
        "METRIS" => Kern::M(Metris::new(md, l2)),
        "CLENSHAW" => Kern::C(Clen),
        "LEAR" => Kern::L(Lear::new(md.n)),
        other => {
            eprintln!("unknown method {other}");
            std::process::exit(2);
        }
    }
}

/// One pass over the points; the outputs go through black_box into a checksum.
#[inline(never)]
fn run_pass<K: Kernel, const L2: bool>(k: &mut K, md: &Model, pts: &[[f64; 3]]) -> f64 {
    let nout = if L2 { 9 } else { 3 };
    let mut sink = 0.0f64;
    for p in pts {
        let p = black_box(*p);
        let o = black_box(k.eval::<L2>(md, p[0], p[1], p[2]));
        for i in 0..nout {
            sink += o[i];
        }
    }
    sink
}

impl Kern {
    fn pass(&mut self, md: &Model, pts: &[[f64; 3]], l2: bool) -> f64 {
        macro_rules! go {
            ($k:expr) => {
                if l2 { run_pass::<_, true>($k, md, pts) } else { run_pass::<_, false>($k, md, pts) }
            };
        }
        match self {
            Kern::A(k) => go!(k),
            Kern::A2(k) => go!(k),
            Kern::B(k) => go!(k),
            Kern::P(k) => go!(k),
            Kern::G(k) => go!(k),
            Kern::M(k) => go!(k),
            Kern::C(k) => go!(k),
            Kern::L(k) => go!(k),
        }
    }
    fn one(&mut self, md: &Model, p: &[f64; 3], l2: bool) -> [f64; 9] {
        macro_rules! go {
            ($k:expr) => {
                if l2 { $k.eval::<true>(md, p[0], p[1], p[2]) } else { $k.eval::<false>(md, p[0], p[1], p[2]) }
            };
        }
        match self {
            Kern::A(k) => go!(k),
            Kern::A2(k) => go!(k),
            Kern::B(k) => go!(k),
            Kern::P(k) => go!(k),
            Kern::G(k) => go!(k),
            Kern::M(k) => go!(k),
            Kern::C(k) => go!(k),
            Kern::L(k) => go!(k),
        }
    }
}

pub(crate) fn read_xyz(path: &str) -> Vec<[f64; 3]> {
    let b = fs::read(path).expect("cannot read point file");
    let mut off = 0;
    let k = rd_u64(&b, &mut off);
    let mut v = vec![[0.0f64; 3]; k];
    for c in 0..3 {
        for p in v.iter_mut() {
            p[c] = rd_f64(&b, &mut off);
        }
    }
    v
}

struct Opts {
    nmax: Option<usize>,
    level: usize,
    blocks: usize,
    target_ms: f64,
    seed: u64,
    methods: Vec<String>,
    method: String,
    cal_iters: usize,
}

fn parse_opts(args: &[String]) -> Opts {
    let mut o = Opts {
        nmax: None, level: 2, blocks: 11, target_ms: 100.0, seed: 2009,
        methods: ALL.iter().map(|s| s.to_string()).collect(), method: String::new(), cal_iters: 2_000_000,
    };
    let mut i = 0;
    while i < args.len() {
        let v = args.get(i + 1).cloned().unwrap_or_default();
        match args[i].as_str() {
            "--nmax" => o.nmax = Some(v.parse().unwrap()),
            "--level" => o.level = v.parse().unwrap(),
            "--blocks" => o.blocks = v.parse().unwrap(),
            "--target-ms" => o.target_ms = v.parse().unwrap(),
            "--seed" => o.seed = v.parse().unwrap(),
            "--methods" => o.methods = v.split(',').map(|s| s.to_string()).collect(),
            "--method" => o.method = v,
            "--cal-iters" => o.cal_iters = v.parse().unwrap(),
            other => {
                eprintln!("unknown option {other}");
                std::process::exit(2);
            }
        }
        i += 2;
    }
    assert!(o.level == 1 || o.level == 2, "--level must be 1 or 2");
    o
}

pub fn main(args: &[String]) {
    match args[1].as_str() {
        "eval" => main_eval(args),
        _ => main_bench(args),
    }
}

/// sph_rs eval <coef.bin> <xyz.bin> <out.bin> --method M [--nmax N] [--level L]
fn main_eval(args: &[String]) {
    if args.len() < 5 {
        eprintln!("usage: sph_rs eval <coef.bin> <xyz.bin> <out.bin> --method M [--nmax N] [--level 1|2]");
        std::process::exit(2);
    }
    let o = parse_opts(&args[5..]);
    let l2 = o.level == 2;
    let md = load_model(&args[2], o.nmax);
    let pts = read_xyz(&args[3]);
    let mut k = make(&o.method, &md, l2);
    let rows = if l2 { 9 } else { 3 };
    let kk = pts.len();
    let mut out = vec![0.0f64; rows * kk];
    for (p, x) in pts.iter().enumerate() {
        let v = k.one(&md, x, l2);
        for qi in 0..rows {
            out[qi * kk + p] = v[qi];
        }
    }
    let mut bytes = Vec::with_capacity(16 + 8 * out.len());
    bytes.extend_from_slice(&(rows as u64).to_le_bytes());
    bytes.extend_from_slice(&(kk as u64).to_le_bytes());
    for v in &out {
        bytes.extend_from_slice(&v.to_le_bytes());
    }
    fs::File::create(&args[4]).and_then(|mut fh| fh.write_all(&bytes)).expect("cannot write output");
}

/// sph_rs bench <coef.bin> <xyz.bin> [options]: one JSON object per line on stdout.
fn main_bench(args: &[String]) {
    if args.len() < 4 {
        eprintln!("usage: sph_rs bench <coef.bin> <xyz.bin> [--nmax N] [--level 1|2] [--blocks B] \
                   [--target-ms T] [--seed S] [--methods M1,M2,..]");
        std::process::exit(2);
    }
    let o = parse_opts(&args[4..]);
    let l2 = o.level == 2;
    let qos = set_qos();
    let t0 = Instant::now();
    let md = load_model(&args[2], o.nmax);
    let load_s = t0.elapsed().as_secs_f64();
    let pts = read_xyz(&args[3]);
    let mut kern: Vec<Kern> = Vec::new();
    let mut setup = Vec::new();
    for name in &o.methods {
        let t1 = Instant::now();
        kern.push(make(name, &md, l2));
        setup.push(t1.elapsed().as_secs_f64());
    }
    let nm = kern.len();
    let setup_js: Vec<String> = o.methods.iter().zip(&setup).map(|(m, s)| format!("\"{m}\": {s:.6e}")).collect();
    println!("{{\"t\": \"setup\", \"N\": {}, \"level\": {}, \"K\": {}, \"blocks\": {}, \"target_ms\": {}, \
              \"qos\": {}, \"clock\": \"thread_cputime\", \"load_s\": {:.6e}, \"setup_s\": {{{}}}}}",
             md.n, o.level, pts.len(), o.blocks, o.target_ms, qos, load_s, setup_js.join(", "));
    let mut rng = Rng(o.seed.wrapping_mul(0x9E37_79B9_7F4A_7C15) | 1);
    let mut checksum = vec![0.0f64; nm];
    // warm-up block (discarded): two passes per method, the second sets the repetitions
    let mut order: Vec<usize> = (0..nm).collect();
    rng.shuffle(&mut order);
    let mut reps = vec![1usize; nm];
    for &mi in &order {
        checksum[mi] += kern[mi].pass(&md, &pts, l2);
        let c0 = thread_cpu();
        checksum[mi] += kern[mi].pass(&md, &pts, l2);
        let dt = (thread_cpu() - c0).max(1e-7);
        reps[mi] = ((o.target_ms * 1e-3 / dt).ceil() as usize).clamp(1, 1_000_000);
        println!("{{\"t\": \"pilot\", \"method\": \"{}\", \"cpu_pass\": {:.6e}, \"reps\": {}}}", o.methods[mi], dt, reps[mi]);
    }
    std::io::stdout().flush().ok();
    for blk in 0..o.blocks {
        rng.shuffle(&mut order);
        for (pos, &mi) in order.iter().enumerate() {
            let cal = calibrate(o.cal_iters);
            let (c0, w0) = (thread_cpu(), Instant::now());
            for _ in 0..reps[mi] {
                checksum[mi] += kern[mi].pass(&md, &pts, l2);
            }
            let cpu = thread_cpu() - c0;
            let wall = w0.elapsed().as_secs_f64();
            let cal2 = calibrate(o.cal_iters);
            println!("{{\"t\": \"m\", \"blk\": {}, \"pos\": {}, \"method\": \"{}\", \"cpu\": {:.9e}, \"wall\": {:.9e}, \
                      \"reps\": {}, \"cal_ns\": {:.6e}, \"cal2_ns\": {:.6e}}}", blk, pos, o.methods[mi], cpu, wall, reps[mi], cal, cal2);
        }
        std::io::stdout().flush().ok();
    }
    let cs: Vec<String> = o.methods.iter().zip(&checksum).map(|(m, s)| format!("\"{m}\": {s:.17e}")).collect();
    println!("{{\"t\": \"end\", \"checksum\": {{{}}}}}", cs.join(", "));
}
