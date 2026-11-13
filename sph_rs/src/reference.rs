//! reference.rs -- accuracy reference: acceleration + gravity-gradient tensor in double-double
//! arithmetic (~32 significant digits) with an extended exponent range, by two independent
//! formulations that serve as a cross-check of each other:
//!
//!   A  spherical: pole-safe Legendre column recursion on Z_nm = Pbar_nm / cos(phi)^min(m,2),
//!      sums per degree in the local (east, north, radial) frame, rotated to Cartesian;
//!   B  Cartesian: Pines' derived Legendre functions A_nm(u) and the complex powers
//!      (x + i y)^m / r^m, tensor assembled directly in Cartesian coordinates.
//!
//! The Legendre values and the powers are carried as X-numbers (Fukushima 2012, J Geod 86:271):
//! a double-double mantissa times 2^(480 i), so that no recursion under- or overflows at any
//! degree or latitude. Values are converted to plain double-double only when they enter a sum.
//!
//!   sph_rs ref <coef.bin> <xyz.bin> <out.bin> [--nmax N] [--threads T]
//!
//! out.bin : u64 rows = 36, u64 K, data[36 * K] row-major (f64):
//!   rows  0.. 8  A, high part     (ax ay az, xx xy xz yy yz zz)
//!   rows  9..17  A, low part      (reference = high + low)
//!   rows 18..26  A minus the central term (-GM x / r^3 and its gradient), rounded
//!   rows 27..35  B minus A, rounded (cross-check)

use crate::bench::read_xyz;
use crate::{load_model, Model};
use std::fs;
use std::io::Write;
use std::ops::{Add, Div, Mul, Neg, Sub};
use std::time::Instant;

// ------------------------------------------------------------------ double-double
#[derive(Clone, Copy, Debug, Default)]
pub struct DD {
    hi: f64,
    lo: f64,
}
#[inline(always)]
fn two_sum(a: f64, b: f64) -> (f64, f64) {
    let s = a + b;
    let bb = s - a;
    (s, (a - (s - bb)) + (b - bb))
}
#[inline(always)]
fn qts(a: f64, b: f64) -> DD {
    let s = a + b;
    DD { hi: s, lo: b - (s - a) }
}
const Z: DD = DD { hi: 0.0, lo: 0.0 };
impl DD {
    #[inline(always)]
    fn f(x: f64) -> DD {
        DD { hi: x, lo: 0.0 }
    }
    /// multiply by an exact power of two
    #[inline(always)]
    fn scale(self, p: f64) -> DD {
        DD { hi: self.hi * p, lo: self.lo * p }
    }
    #[inline(always)]
    fn mulf(self, b: f64) -> DD {
        let p = self.hi * b;
        let e = self.hi.mul_add(b, -p) + self.lo * b;
        qts(p, e)
    }
    fn sqrt(self) -> DD {
        if self.hi <= 0.0 {
            return Z;
        }
        let q = self.hi.sqrt();
        let r = self - DD::f(q) * DD::f(q);
        DD::f(q) + DD::f(r.hi / (2.0 * q))
    }
}
impl Add for DD {
    type Output = DD;
    #[inline(always)]
    fn add(self, b: DD) -> DD {
        let (s, e) = two_sum(self.hi, b.hi);
        let (t, f) = two_sum(self.lo, b.lo);
        let d = qts(s, e + t);
        qts(d.hi, d.lo + f)
    }
}
impl Neg for DD {
    type Output = DD;
    #[inline(always)]
    fn neg(self) -> DD {
        DD { hi: -self.hi, lo: -self.lo }
    }
}
impl Sub for DD {
    type Output = DD;
    #[inline(always)]
    fn sub(self, b: DD) -> DD {
        self + (-b)
    }
}
impl Mul for DD {
    type Output = DD;
    #[inline(always)]
    fn mul(self, b: DD) -> DD {
        let p = self.hi * b.hi;
        let e = self.hi.mul_add(b.hi, -p) + (self.hi * b.lo + self.lo * b.hi);
        qts(p, e)
    }
}
impl Div for DD {
    type Output = DD;
    fn div(self, b: DD) -> DD {
        let q1 = self.hi / b.hi;
        let r = self - b.mulf(q1);
        let q2 = r.hi / b.hi;
        let r = r - b.mulf(q2);
        let q3 = r.hi / b.hi;
        qts(q1, q2) + DD::f(q3)
    }
}
impl Mul<f64> for DD {
    type Output = DD;
    #[inline(always)]
    fn mul(self, b: f64) -> DD {
        self.mulf(b)
    }
}

// ------------------------------------------------------------------ X-numbers
const BIG: f64 = f64::from_bits(((1023 + 480) as u64) << 52); // 2^480, the exponent step
const BIGI: f64 = f64::from_bits(((1023 - 480) as u64) << 52); // 2^-480
const WIN: f64 = f64::from_bits(((1023 + 240) as u64) << 52); // 2^240, normalisation window
const WINI: f64 = f64::from_bits(((1023 - 240) as u64) << 52); // 2^-240

/// value = v * 2^(480 e), |v.hi| in [2^-240, 2^240) unless zero. The window is as wide as the
/// step, so the representation is unique: two values whose exponents differ by 2 or more differ
/// in magnitude by more than 2^480 and the smaller one cannot affect a double-double sum.
#[derive(Clone, Copy, Default)]
struct X {
    v: DD,
    e: i32,
}
impl X {
    #[inline(always)]
    fn new(v: DD, e: i32) -> X {
        let mut x = X { v, e };
        if v.hi == 0.0 {
            return X { v: Z, e: 0 };
        }
        while x.v.hi.abs() >= WIN {
            x.v = x.v.scale(BIGI);
            x.e += 1;
        }
        while x.v.hi.abs() < WINI {
            x.v = x.v.scale(BIG);
            x.e -= 1;
        }
        x
    }
    #[inline(always)]
    fn zero(&self) -> bool {
        self.v.hi == 0.0
    }
    /// to plain double-double (0 below 2^-720 ~ 1e-217, negligible: every sum holds the O(1) degree-0 term)
    #[inline(always)]
    fn dd(self) -> DD {
        match self.e {
            0 => self.v,
            -1 => self.v.scale(BIGI),
            1 => self.v.scale(BIG),
            e if e < -1 => Z,
            _ => DD::f(f64::INFINITY),
        }
    }
    /// a x + b y
    #[inline(always)]
    fn lin2(a: DD, x: X, b: DD, y: X) -> X {
        if x.zero() {
            return X::new(b * y.v, y.e);
        }
        if y.zero() {
            return X::new(a * x.v, x.e);
        }
        let e = x.e.max(y.e);
        let al = |v: DD, ev: i32| match e - ev {
            0 => v,
            1 => v.scale(BIGI),
            _ => Z,
        };
        X::new(al(a * x.v, x.e) + al(b * y.v, y.e), e)
    }
}
/// (product of an X mantissa and a DD) * 2^(480 e) as plain double-double
#[inline(always)]
fn xdd(v: DD, e: i32) -> DD {
    X::new(v, e).dd()
}

// ------------------------------------------------------------------ tables (triangle, DD)
struct Tab {
    f: Vec<DD>,
    g: Vec<DD>,
    h: Vec<DD>,
    k: Vec<DD>,
    j: Vec<DD>,
    r11: Vec<DD>,
    r12: Vec<DD>,
    r22: Vec<DD>,
}
#[inline(always)]
fn tri(n: usize, m: usize) -> usize {
    n * (n + 1) / 2 + m
}
fn prod(v: &[f64]) -> DD {
    v.iter().fold(DD::f(1.0), |a, &b| a * DD::f(b))
}
fn sq(num: &[f64], den: &[f64]) -> DD {
    (prod(num) / prod(den)).sqrt()
}
fn tab(nt: usize) -> Tab {
    let sz = tri(nt + 1, 0);
    let mut t = Tab {
        f: vec![Z; nt + 1], g: vec![Z; sz], h: vec![Z; sz], k: vec![Z; sz], j: vec![Z; sz],
        r11: vec![Z; sz], r12: vec![Z; sz], r22: vec![Z; sz],
    };
    for ni in 1..=nt {
        let n = ni as f64;
        t.f[ni] = sq(&[2.0 * n + 1.0, if ni == 1 { 2.0 } else { 1.0 }], &[2.0 * n]);
    }
    for ni in 0..=nt {
        for mi in 0..=ni {
            let (n, m) = (ni as f64, mi as f64);
            let eps = if mi == 0 { 0.5 } else { 1.0 };
            let i = tri(ni, mi);
            if ni > mi {
                t.g[i] = sq(&[2.0 * n + 1.0, 2.0 * n - 1.0], &[n + m, n - m]);
            }
            if ni > mi + 1 {
                t.h[i] = sq(&[2.0 * n + 1.0, n - m - 1.0, n + m - 1.0], &[2.0 * n - 3.0, n + m, n - m]);
            }
            t.k[i] = sq(&[eps, n - m, n + m + 1.0], &[1.0]);
            if ni >= mi + 2 {
                t.j[i] = sq(&[eps, n - m, n - m - 1.0, n + m + 1.0, n + m + 2.0], &[1.0]);
            }
            t.r11[i] = sq(&[eps, 2.0 * n + 1.0, n + m + 1.0, n + m + 2.0], &[2.0 * n + 3.0]);
            t.r12[i] = sq(&[eps, 2.0 * n + 1.0, n - m, n + m + 1.0, n + m + 2.0, n + m + 3.0], &[2.0 * n + 3.0]);
            t.r22[i] = sq(&[eps, 2.0 * n + 1.0, n + m + 1.0, n + m + 2.0, n + m + 3.0, n + m + 4.0], &[2.0 * n + 5.0]);
        }
    }
    t
}

// ------------------------------------------------------------------ geometry
struct Geo {
    x: [DD; 3],
    r: DD,
    u: DD,
    c: DD,
    cl: DD,
    sl: DD,
}
fn geo(p: &[f64; 3]) -> Geo {
    let x = [DD::f(p[0]), DD::f(p[1]), DD::f(p[2])];
    let rho2 = x[0] * x[0] + x[1] * x[1];
    let rho = rho2.sqrt();
    let r = (rho2 + x[2] * x[2]).sqrt();
    Geo { x, r, u: x[2] / r, c: rho / r, cl: x[0] / rho, sl: x[1] / rho }
}

/// central term -GM x / r^3 and its gradient, scaled by C00
fn central(md: &Model, g: &Geo) -> [DD; 9] {
    let gm = DD::f(md.gm) * DD::f(md.c[0]);
    let r2 = g.r * g.r;
    let r3 = r2 * g.r;
    let r5 = r3 * r2;
    let mut o = [Z; 9];
    for i in 0..3 {
        o[i] = -(gm * g.x[i] / r3);
    }
    let pairs = [(0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2)];
    for (k, &(i, j)) in pairs.iter().enumerate() {
        let mut v = DD::f(3.0) * g.x[i] * g.x[j] / r5;
        if i == j {
            v = v - DD::f(1.0) / r3;
        }
        o[3 + k] = gm * v;
    }
    o
}

// ------------------------------------------------------------------ A: spherical, pole-safe
fn form_a(md: &Model, t: &Tab, g: &Geo) -> [DD; 9] {
    let nmax = md.n;
    let (u, c) = (g.u, g.c);
    let c2 = c * c;
    let q = DD::f(md.a) / g.r;
    let mut cm = vec![Z; nmax + 1];
    let mut sm = vec![Z; nmax + 1];
    cm[0] = DD::f(1.0);
    for m in 1..=nmax {
        cm[m] = cm[m - 1] * g.cl - sm[m - 1] * g.sl;
        sm[m] = sm[m - 1] * g.cl + cm[m - 1] * g.sl;
    }
    let mut p1 = vec![X::default(); nmax + 3]; // row n-1
    let mut p2 = vec![X::default(); nmax + 3]; // row n-2
    let mut zd = vec![Z; nmax + 3];
    let mut acc = [Z; 9];
    let mut w = DD::f(1.0);
    for n in 0..=nmax {
        let mut cur = vec![X::default(); nmax + 3];
        if n == 0 {
            cur[0] = X::new(DD::f(1.0), 0);
        } else {
            for m in 0..n {
                let i = tri(n, m);
                cur[m] = X::lin2(t.g[i] * u, p1[m], -t.h[i], p2[m]);
            }
            let zs = if n <= 2 { DD::f(1.0) } else { c };
            let pv = p1[n - 1];
            cur[n] = X::new(pv.v * (t.f[n] * zs), pv.e);
        }
        for m in 0..n + 3 {
            zd[m] = cur[m].dd();
        }
        // sums over the order: [P X, dP X, m Q Y, d2P X, en Y, ee X]
        let mut o = [Z; 6];
        for m in 0..=n {
            let fm = DD::f(m as f64);
            let i = tri(n, m);
            let (cn, sn) = (DD::f(md.c[n * md.n1 + m]), DD::f(md.s[n * md.n1 + m]));
            let xx = cn * cm[m] + sn * sm[m];
            let yy = sn * cm[m] - cn * sm[m];
            let (p0, q0, r0, pp1, q1, pp2) = match m {
                0 => (zd[0], Z, zd[0], c * zd[1], zd[1], c2 * zd[2]),
                1 => (c * zd[1], zd[1], zd[1], c2 * zd[2], c * zd[2], c2 * zd[3]),
                _ => (c2 * zd[m], c * zd[m], zd[m], c2 * zd[m + 1], c * zd[m + 1], c2 * zd[m + 2]),
            };
            let (k, j) = (t.k[i], t.j[i]);
            let dp = k * pp1 - fm * u * q0;
            o[0] = o[0] + p0 * xx;
            o[1] = o[1] + dp * xx;
            o[2] = o[2] + fm * q0 * yy;
            let mm1 = fm * DD::f(m as f64 - 1.0);
            let kq1 = k * q1;
            let d2p = mm1 * r0 - fm * fm * p0 - DD::f(2.0 * m as f64 + 1.0) * u * kq1 + j * pp2;
            o[3] = o[3] + d2p * xx;
            o[4] = o[4] + fm * (kq1 - DD::f(m as f64 - 1.0) * u * r0) * yy;
            o[5] = o[5] + (-(mm1 * r0) - fm * p0 - u * kq1) * xx;
        }
        let a1 = DD::f(n as f64 + 1.0);
        let a2 = DD::f(n as f64 + 2.0);
        acc[0] = acc[0] + w * a1 * o[0];
        acc[1] = acc[1] + w * o[1];
        acc[2] = acc[2] + w * o[2];
        acc[3] = acc[3] + w * a1 * a2 * o[0];
        acc[4] = acc[4] + w * o[3];
        acc[5] = acc[5] + w * o[5];
        acc[6] = acc[6] + w * o[4];
        acc[7] = acc[7] + w * a2 * o[2];
        acc[8] = acc[8] + w * a2 * o[1];
        w = w * q;
        p2 = std::mem::replace(&mut p1, cur);
    }
    let b1 = DD::f(md.gm) / (g.r * g.r);
    let b2 = b1 / g.r;
    let gl = [b1 * acc[2], b1 * acc[1], -(b1 * acc[0])];
    let tl6 = [b2 * (acc[5] - acc[0]), b2 * acc[6], -(b2 * acc[7]), b2 * (acc[4] - acc[0]), -(b2 * acc[8]), b2 * acc[3]];
    // local (e, n, r) -> Cartesian
    let (u, c, cl, sl) = (g.u, g.c, g.cl, g.sl);
    let rr = [[-sl, cl, Z], [-(u * cl), -(u * sl), c], [c * cl, c * sl, u]];
    let mut out = [Z; 9];
    for a in 0..3 {
        out[a] = rr[0][a] * gl[0] + rr[1][a] * gl[1] + rr[2][a] * gl[2];
    }
    let tl = [[tl6[0], tl6[1], tl6[2]], [tl6[1], tl6[3], tl6[4]], [tl6[2], tl6[4], tl6[5]]];
    let mut mm = [[Z; 3]; 3];
    for i in 0..3 {
        for b in 0..3 {
            mm[i][b] = tl[i][0] * rr[0][b] + tl[i][1] * rr[1][b] + tl[i][2] * rr[2][b];
        }
    }
    let pairs = [(0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2)];
    for (k, &(a, b)) in pairs.iter().enumerate() {
        out[3 + k] = rr[0][a] * mm[0][b] + rr[1][a] * mm[1][b] + rr[2][a] * mm[2][b];
    }
    out
}

// ------------------------------------------------------------------ B: Pines, Cartesian
fn form_b(md: &Model, t: &Tab, g: &Geo) -> [DD; 9] {
    let nmax = md.n;
    let (s, tt, u) = (g.x[0] / g.r, g.x[1] / g.r, g.u);
    let q = DD::f(md.a) / g.r;
    // (s + i t)^m, shared exponent; index m + 2 (indices 0, 1 = zero)
    let mut pr = vec![Z; nmax + 3];
    let mut pi = vec![Z; nmax + 3];
    let mut pe = vec![0i32; nmax + 3];
    pr[2] = DD::f(1.0);
    for m in 1..=nmax {
        let (a, b, e) = (pr[m + 1], pi[m + 1], pe[m + 1]);
        let (re, im) = (a * s - b * tt, b * s + a * tt);
        let big = re.hi.abs().max(im.hi.abs());
        let xr = X::new(DD::f(big), e); // exponent decision on the larger part
        let sc = match xr.e - e {
            0 => 1.0,
            -1 => BIG,
            1 => BIGI,
            _ => 1.0,
        };
        pr[m + 2] = re.scale(sc);
        pi[m + 2] = im.scale(sc);
        pe[m + 2] = xr.e;
    }
    let w_ = nmax + 5;
    let mut rows: Vec<Vec<X>> = (0..4).map(|_| vec![X::default(); w_]).collect();
    let hrow = |rows: &mut Vec<Vec<X>>, n: usize| {
        let mut an = vec![X::default(); w_];
        for m in 0..n {
            let i = tri(n, m);
            an[m] = X::lin2(t.g[i] * u, rows[(n + 3) % 4][m], -t.h[i], rows[(n + 2) % 4][m]);
        }
        let pv = rows[(n + 3) % 4][n - 1];
        an[n] = X::new(pv.v * t.f[n], pv.e);
        rows[n % 4] = an;
    };
    rows[0][0] = X::new(DD::f(1.0), 0);
    hrow(&mut rows, 1);
    hrow(&mut rows, 2);
    let mut acc = [Z; 13];
    let mut w = DD::f(1.0);
    for n in 0..=nmax {
        let mut o = [Z; 13];
        {
            let a0 = &rows[n % 4];
            let b1 = &rows[(n + 1) % 4];
            let b2 = &rows[(n + 2) % 4];
            for m in 0..=n {
                let i = tri(n, m);
                let fm = DD::f(m as f64);
                let (c_, s_) = (DD::f(md.c[n * md.n1 + m]), DD::f(md.s[n * md.n1 + m]));
                let (d0, de) = (c_ * pr[m + 2] + s_ * pi[m + 2], pe[m + 2]);
                let (e1, f1, ee) = (fm * (c_ * pr[m + 1] + s_ * pi[m + 1]), fm * (s_ * pr[m + 1] - c_ * pi[m + 1]), pe[m + 1]);
                let mm1 = fm * DD::f(m as f64 - 1.0);
                let (g1, h1, ge) = (mm1 * (c_ * pr[m] + s_ * pi[m]), mm1 * (s_ * pr[m] - c_ * pi[m]), pe[m]);
                let v0 = a0[m];
                let v1 = (t.k[i] * a0[m + 1].v, a0[m + 1].e);
                let v2 = (-(t.r11[i] * b1[m + 1].v), b1[m + 1].e);
                let v3 = (t.j[i] * a0[m + 2].v, a0[m + 2].e);
                let v4 = (-(t.r12[i] * b1[m + 2].v), b1[m + 2].e);
                let v5 = (t.r22[i] * b2[m + 2].v, b2[m + 2].e);
                o[0] = o[0] + xdd(v1.0 * d0, v1.1 + de);
                o[1] = o[1] + xdd(v2.0 * d0, v2.1 + de);
                o[2] = o[2] + xdd(v0.v * e1, v0.e + ee);
                o[3] = o[3] + xdd(v0.v * f1, v0.e + ee);
                o[4] = o[4] + xdd(v3.0 * d0, v3.1 + de);
                o[5] = o[5] + xdd(v4.0 * d0, v4.1 + de);
                o[6] = o[6] + xdd(v5.0 * d0, v5.1 + de);
                o[7] = o[7] + xdd(v1.0 * e1, v1.1 + ee);
                o[8] = o[8] + xdd(v1.0 * f1, v1.1 + ee);
                o[9] = o[9] + xdd(v2.0 * e1, v2.1 + ee);
                o[10] = o[10] + xdd(v2.0 * f1, v2.1 + ee);
                o[11] = o[11] + xdd(v0.v * g1, v0.e + ge);
                o[12] = o[12] + xdd(v0.v * h1, v0.e + ge);
            }
        }
        for i in 0..13 {
            acc[i] = acc[i] + w * o[i];
        }
        w = w * q;
        if n + 3 <= nmax + 2 {
            hrow(&mut rows, n + 3);
        }
    }
    let r1 = DD::f(md.gm) / (g.r * g.r);
    let (a3, a4, a1, a2) = (acc[0] * r1, acc[1] * r1, acc[2] * r1, acc[3] * r1);
    let mut out = [a1 + s * a4, a2 + tt * a4, a3 + u * a4, Z, Z, Z, Z, Z, Z];
    let r2 = r1 / g.r;
    let (a33, a34, a44) = (acc[4] * r2, acc[5] * r2, acc[6] * r2);
    let (a13, a23, a14, a24) = (acc[7] * r2, acc[8] * r2, acc[9] * r2, acc[10] * r2);
    let (a11, a12) = (acc[11] * r2, acc[12] * r2);
    let xh = [s, tt, u];
    let ai4 = [a14, a24, a34];
    let aij = [[a11, a12, a13], [a12, -a11, a23], [a13, a23, a33]];
    let a4r = a4 / g.r;
    let pairs = [(0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2)];
    for (col, &(i, j)) in pairs.iter().enumerate() {
        let mut v = aij[i][j] + xh[j] * ai4[i] + xh[i] * ai4[j] + xh[i] * xh[j] * a44;
        if i == j {
            v = v + a4r;
        }
        out[3 + col] = v;
    }
    out
}

// ------------------------------------------------------------------ driver
pub fn main(args: &[String]) {
    if args.len() < 5 {
        eprintln!("usage: sph_rs ref <coef.bin> <xyz.bin> <out.bin> [--nmax N] [--threads T]");
        std::process::exit(2);
    }
    let (mut nmax, mut threads) = (None, 4usize);
    let mut i = 5;
    while i < args.len() {
        let v = args.get(i + 1).cloned().unwrap_or_default();
        match args[i].as_str() {
            "--nmax" => nmax = Some(v.parse().unwrap()),
            "--threads" => threads = v.parse().unwrap(),
            o => {
                eprintln!("unknown option {o}");
                std::process::exit(2);
            }
        }
        i += 2;
    }
    let t0 = Instant::now();
    let md = load_model(&args[2], nmax);
    let pts = read_xyz(&args[3]);
    let t = tab(md.n + 2);
    let t_tab = t0.elapsed().as_secs_f64();
    let kk = pts.len();
    let mut res: Vec<[[DD; 9]; 3]> = vec![[[Z; 9]; 3]; kk];
    let chunk = (kk + threads - 1) / threads.max(1);
    std::thread::scope(|sc| {
        for (ci, part) in res.chunks_mut(chunk.max(1)).enumerate() {
            let (md, t, pts) = (&md, &t, &pts);
            sc.spawn(move || {
                for (k, slot) in part.iter_mut().enumerate() {
                    let g = geo(&pts[ci * chunk + k]);
                    *slot = [form_a(md, t, &g), form_b(md, t, &g), central(md, &g)];
                }
            });
        }
    });
    let rows = 36;
    let mut out = vec![0.0f64; rows * kk];
    for (p, [a, b, cen]) in res.iter().enumerate() {
        for q in 0..9 {
            out[q * kk + p] = a[q].hi;
            out[(9 + q) * kk + p] = a[q].lo;
            out[(18 + q) * kk + p] = (a[q] - cen[q]).hi;
            out[(27 + q) * kk + p] = (b[q] - a[q]).hi;
        }
    }
    let mut bytes = Vec::with_capacity(16 + 8 * out.len());
    bytes.extend_from_slice(&(rows as u64).to_le_bytes());
    bytes.extend_from_slice(&(kk as u64).to_le_bytes());
    for v in &out {
        bytes.extend_from_slice(&v.to_le_bytes());
    }
    fs::File::create(&args[4]).and_then(|mut fh| fh.write_all(&bytes)).expect("cannot write output");
    eprintln!("ref: N={} K={} tables {:.2}s total {:.2}s", md.n, kk, t_tab, t0.elapsed().as_secs_f64());
}
