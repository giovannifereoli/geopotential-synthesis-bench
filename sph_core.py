"""
sph_core.py -- the five spherical-harmonic synthesis programmes of

    G. Fantino & S. Casotto (2009), "Methods of harmonic synthesis for global geopotential
    models and their first-, second- and third-order gradients", J. Geodesy 83:595-619

re-implemented in NumPy so that their precision, accuracy and relative speed can be studied.

    LEGENDRE_a  forward column recursion (FCR) on fully normalised associated Legendre
                functions (fnALFs), spherical coordinates, vectorised over the points
                (compiled Fortran90 in the paper)
    LEGENDRE_b  the same algorithm called point by point from an interpreted loop, with
                whole-triangle array operations per point (PoliMi DLL under Matlab)
    PINES       FCR on fully normalised Helmholtz polynomials, Cartesian (Pines 1973)
    METRIS      recursion on fully normalised solid spherical harmonics, the derivatives
                obtained by acting on the coefficients (Metris et al. 1999), Cartesian
    CLENSHAW    Clenshaw summation over the degree on fnALFs, either for scattered points
                (vectorised over points and orders) or per parallel on a grid

Every routine is dtype-generic: dtype=np.float64 is the paper's REAL*8, and
dtype=np.longdouble (80-bit extended, ~19 digits on x86-64) stands in for REAL*16.

All routines return the ten quantities of the paper in the local north-oriented frame
(e = east, n = north, r = radial), in the order QUANTS:

    V, V_lambda, V_phi, V_r, V_lambda lambda, V_lambda phi, V_phi phi, V_lambda r, V_phi r, V_rr

PINES and METRIS also return the Cartesian quantities (order CART) with cartesian=True.

Run `python sph_core.py` for a self-test (all five programmes agree at N = 30, Laplace).
"""

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")  # single thread, as the paper's Pentium IV

import time
from pathlib import Path

import numpy as np

# ============================================================
# CONSTANTS AND LABELS
# ============================================================
GM_EARTH = 3.986004415e14  # [m^3 s^-2]  (EGM96)
A_EARTH = 6378136.3  # [m]
J2_EARTH = 1.0826266835531513e-3  # = -sqrt(5) * C20bar

METHODS = ("LEGENDRE_a", "LEGENDRE_b", "PINES", "METRIS", "CLENSHAW")
QUANTS = ("V", "Vl", "Vp", "Vr", "Vll", "Vlp", "Vpp", "Vlr", "Vpr", "Vrr")
CART = ("V", "Vx", "Vy", "Vz", "Vxx", "Vxy", "Vxz", "Vyy", "Vyz", "Vzz")
QUANT_TEX = {
    "V": r"$V$", "Vl": r"$V_{\lambda}$", "Vp": r"$V_{\varphi}$", "Vr": r"$V_{r}$",
    "Vll": r"$V_{\lambda\lambda}$", "Vlp": r"$V_{\lambda\varphi}$",
    "Vpp": r"$V_{\varphi\varphi}$", "Vlr": r"$V_{\lambda r}$", "Vpr": r"$V_{\varphi r}$",
    "Vrr": r"$V_{rr}$",
    "Vx": r"$V_{x}$", "Vy": r"$V_{y}$", "Vz": r"$V_{z}$", "Vxx": r"$V_{xx}$",
    "Vxy": r"$V_{xy}$", "Vxz": r"$V_{xz}$", "Vyy": r"$V_{yy}$", "Vyz": r"$V_{yz}$",
    "Vzz": r"$V_{zz}$",
}
QUANT_UNIT = {"V": r"m$^2$\,s$^{-2}$", "g": r"m\,s$^{-2}$", "T": r"s$^{-2}$"}

# default number of points per vectorised block (memory ~ chunk * N * 8 B * #arrays)
DEFAULT_CHUNK = {"LEGENDRE_a": 256, "PINES": 256, "METRIS": 1024, "CLENSHAW": 256}


# ============================================================
# RECURSION COEFFICIENTS (cached per N and dtype)
# ============================================================
_TABLES = {}


def tables(N, dtype=np.float64):
    """
    Recursion / ratio coefficients as (N+5, N+5) arrays indexed [n, m] (zero outside their domain).

    f[n]      sectorial  P_nn = f_n cos(phi) P_n-1,n-1           sqrt((2n+1)(1+d_1n)/(2n))
    g, h      column     P_nm = g_nm sin(phi) P_n-1,m - h_nm P_n-2,m
    k, j      d/dphi     P_nm/P_n,m+1 and P_nm/P_n,m+2 normalisation ratios
    R11,R12,R22  Pines   Abar_nm normalisation ratios to A_n+1,m+1 / A_n+1,m+2 / A_n+2,m+2
    alpha,beta,gamma     Metris d/dz, (d/dx, d/dy) raising and lowering coefficients
    msr2      the second coefficient of Pines' modified sub-sectorial recursion
    """
    key = (int(N), np.dtype(dtype).name)
    T = _TABLES.get(key)
    if T is not None:
        return T
    L = int(N) + 5
    dt = np.dtype(dtype).type
    i = np.arange(L)
    n = np.repeat(i[:, None], L, 1).astype(dtype)
    m = np.repeat(i[None, :], L, 0).astype(dtype)
    eps = np.where(m == 0, dt(0.5), dt(1.0))  # (2 - delta_0m) / 2
    one, two = dt(1), dt(2)

    def root(x, mask):
        return np.sqrt(np.where(mask, x, dt(0)))

    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        nn = np.arange(L).astype(dtype)
        f = np.zeros(L, dtype)
        f[1:] = np.sqrt((two * nn[1:] + one) * (one + (nn[1:] == 1)) / (two * nn[1:]))
        T = {
            "f": f,
            "g": root((2 * n + 1) * (2 * n - 1) / ((n + m) * (n - m)), n > m),
            "h": root((2 * n + 1) * (n - m - 1) * (n + m - 1) / ((2 * n - 3) * (n + m) * (n - m)),
                      n > m + 1),
            "k": root(eps * (n - m) * (n + m + 1), n >= m),
            "j": root(eps * (n - m) * (n - m - 1) * (n + m + 1) * (n + m + 2), n >= m + 2),
            "R11": root(eps * (2 * n + 1) * (n + m + 1) * (n + m + 2) / (2 * n + 3), n >= m),
            "R12": root(eps * (2 * n + 1) * (n - m) * (n + m + 1) * (n + m + 2) * (n + m + 3)
                        / (2 * n + 3), n >= m),
            "R22": root(eps * (2 * n + 1) * (n + m + 1) * (n + m + 2) * (n + m + 3) * (n + m + 4)
                        / (2 * n + 5), n >= m),
            "alpha": root((n - m + 1) * (n + m + 1) * (2 * n + 1) / (2 * n + 3), n >= m),
            "gamma": root(np.where(m == 1, two, one) * (2 * n + 1) * (n - m + 2) * (n - m + 1)
                          / (2 * n + 3), (n >= m) & (m >= 1)),
            "msr2": root(eps * (2 * n + 1) * (n - m) * (n - m - 1) / (2 * n - 1), n >= m + 2),
        }
    T["beta"] = T["R11"]  # same normalisation ratio (Metris raising operator)
    _TABLES[key] = T
    return T


# ============================================================
# GRAVITY MODEL
# ============================================================
class GravityModel:
    """Fully normalised Stokes coefficients C[n, m], S[n, m] (0 <= m <= n <= N), GM, a."""

    def __init__(self, C, S, gm=GM_EARTH, a=A_EARTH, name="model"):
        C = np.tril(np.array(C, dtype=np.float64))
        S = np.tril(np.array(S, dtype=np.float64))
        S[:, 0] = 0.0
        self.C, self.S = C, S
        self.gm, self.a, self.name = float(gm), float(a), name
        self._cache = {}

    @property
    def N(self):
        return self.C.shape[0] - 1

    def truncate(self, N):
        return GravityModel(self.C[: N + 1, : N + 1], self.S[: N + 1, : N + 1], self.gm,
                            self.a, f"{self.name} (N={N})")

    def coef(self, dtype=np.float64):
        """(C, S) cast to dtype (exact upcast of the double coefficients)."""
        key = ("coef", np.dtype(dtype).name)
        if key not in self._cache:
            self._cache[key] = (self.C.astype(dtype), self.S.astype(dtype))
        return self._cache[key]

    @classmethod
    def synthetic(cls, N=360, seed=2009):
        """
        Earth-like synthetic model: Kaula's rule sigma_n = 1e-5 / n^2 for every coefficient,
        with the EGM96 values of C20, C22, S22, C30 (the real low-degree field), C00 = 1.
        """
        rng = np.random.default_rng(seed)
        C = np.zeros((N + 1, N + 1))
        S = np.zeros((N + 1, N + 1))
        for n in range(2, N + 1):
            sig = 1e-5 / n**2
            C[n, : n + 1] = rng.normal(0.0, sig, n + 1)
            S[n, 1 : n + 1] = rng.normal(0.0, sig, n)
        C[0, 0] = 1.0
        low = {(2, 0): (-4.84165371736e-4, 0.0), (2, 1): (0.0, 0.0),
               (2, 2): (2.43914352398e-6, -1.40016683654e-6), (3, 0): (9.57254173792e-7, 0.0)}
        for (n, m), (c, s) in low.items():
            if n <= N:
                C[n, m], S[n, m] = c, s
        return cls(C, S, GM_EARTH, A_EARTH, f"synthetic Kaula (N={N}, seed={seed})")

    @classmethod
    def from_gfc(cls, path, N=None):
        """Read an ICGEM .gfc file (e.g. EGM96.gfc, EGM2008.gfc) up to degree N."""
        path = Path(path)
        gm, a, nmax_file, rows = GM_EARTH, A_EARTH, None, []
        with path.open("r", errors="replace") as fh:
            in_head = True
            for line in fh:
                tok = line.split()
                if not tok:
                    continue
                if in_head:
                    key = tok[0].lower()
                    if key == "earth_gravity_constant":
                        gm = float(tok[1].replace("D", "E").replace("d", "e"))
                    elif key == "radius":
                        a = float(tok[1].replace("D", "E").replace("d", "e"))
                    elif key == "max_degree":
                        nmax_file = int(tok[1])
                    elif key == "end_of_head":
                        in_head = False
                    continue
                if tok[0] not in ("gfc", "gfct"):
                    continue
                n, m = int(tok[1]), int(tok[2])
                if N is not None and n > N:
                    continue
                rows.append((n, m, float(tok[3].replace("D", "E")), float(tok[4].replace("D", "E"))))
        nmax = max(r[0] for r in rows) if N is None else int(N)
        if nmax_file is not None and N is None:
            nmax = min(nmax, nmax_file)
        C = np.zeros((nmax + 1, nmax + 1))
        S = np.zeros((nmax + 1, nmax + 1))
        for n, m, c, s in rows:
            if n <= nmax:
                C[n, m], S[n, m] = c, s
        return cls(C, S, gm, a, f"{path.stem} (N={nmax})")


# ============================================================
# COORDINATES AND FRAMES
# ============================================================
def sph_to_cart(lam, phi, r):
    c = np.cos(phi)
    return r * c * np.cos(lam), r * c * np.sin(lam), r * np.sin(phi)


def local_frame(lam, phi):
    """(K, 3, 3) rotation, rows = east, north, radial unit vectors in Cartesian components."""
    sl, cl, sp, cp = np.sin(lam), np.cos(lam), np.sin(phi), np.cos(phi)
    R = np.empty(np.shape(lam) + (3, 3), dtype=np.result_type(lam, phi))
    R[..., 0, 0], R[..., 0, 1], R[..., 0, 2] = -sl, cl, 0.0
    R[..., 1, 0], R[..., 1, 1], R[..., 1, 2] = -sp * cl, -sp * sl, cp
    R[..., 2, 0], R[..., 2, 1], R[..., 2, 2] = cp * cl, cp * sl, sp
    return R


def cart_to_local(cart, lam, phi):
    """Cartesian (10, K) in CART order -> local (10, K) in QUANTS order (Tables 1-2)."""
    R = local_frame(lam, phi)
    g = np.moveaxis(cart[1:4], 0, -1)
    gl = np.einsum("kij,kj->ki", R, g)
    xx, xy, xz, yy, yz, zz = cart[4:10]
    T = np.stack([np.stack([xx, xy, xz], -1), np.stack([xy, yy, yz], -1),
                  np.stack([xz, yz, zz], -1)], -2)
    Tl = np.einsum("kia,kab,kjb->kij", R, T, R)
    return np.stack([cart[0], gl[:, 0], gl[:, 1], gl[:, 2], Tl[:, 0, 0], Tl[:, 0, 1],
                     Tl[:, 1, 1], Tl[:, 0, 2], Tl[:, 1, 2], Tl[:, 2, 2]])


def spherical_to_local(D, r, phi):
    """
    Partial derivatives w.r.t. (lambda, phi, r), order
    [V, d_l, d_p, d_r, d_ll, d_lp, d_pp, d_lr, d_pr, d_rr], -> local frame (Table 3).
    """
    V, dl, dp, dr, dll, dlp, dpp, dlr, dpr, drr = D
    c, s = np.cos(phi), np.sin(phi)
    t = s / c
    rc = r * c
    out = np.empty_like(D)
    out[0] = V
    out[1] = dl / rc
    out[2] = dp / r
    out[3] = dr
    out[4] = (dll + rc * c * dr - s * c * dp) / (rc * rc)
    out[5] = (dlp + t * dl) / (r * rc)
    out[6] = (dpp + r * dr) / (r * r)
    out[7] = (dlr - dl / r) / rc
    out[8] = (dpr - dp / r) / r
    out[9] = drr
    return out


def _as(x, dtype):
    return np.ascontiguousarray(np.atleast_1d(np.asarray(x)).astype(dtype))


def _chunks(K, size):
    for s in range(0, K, size):
        yield slice(s, min(K, s + size))


# ============================================================
# 1a. LEGENDRE_a -- FCR on fnALFs, vectorised over the points
# ============================================================
def _legendre_coef(model, dtype):
    """Per degree n: (2(n+1), 3) matrix, columns [C;S], [mS;-mC], [m^2 C; m^2 S]."""
    key = ("legendre", np.dtype(dtype).name)
    if key not in model._cache:
        C, S = model.coef(dtype)
        W = []
        for n in range(model.N + 1):
            m = np.arange(n + 1).astype(dtype)
            c, s = C[n, : n + 1], S[n, : n + 1]
            W.append(np.ascontiguousarray(np.stack(
                [np.concatenate([c, s]), np.concatenate([m * s, -m * c]),
                 np.concatenate([m * m * c, m * m * s])], axis=1)))
        model._cache[key] = W
    return model._cache[key]


def _legendre_block(model, lam, phi, r, dtype):
    N, mu = model.N, model.gm
    T = tables(N, dtype)
    W = _legendre_coef(model, dtype)
    f, g, h, k, j = T["f"], T["g"], T["h"], T["k"], T["j"]
    K = lam.size
    u, c = np.sin(phi), np.cos(phi)
    t = u / c
    sec2 = 1.0 / (c * c)
    q = model.a / r
    mm = np.arange(N + 1).astype(dtype)
    ang = lam[:, None] * mm[None, :]
    cosm, sinm = np.cos(ang), np.sin(ang)
    mt = mm[None, :] * t[:, None]  # m tan(phi)
    d2c = mt * mt - mm[None, :] * sec2[:, None]  # m^2 tan^2 - m sec^2
    k2t = (2 * mm + 1)[None, :] * t[:, None]  # (2m+1) tan(phi)
    ucol = u[:, None]

    bufs = [np.zeros((K, N + 3), dtype) for _ in range(3)]
    acc = np.zeros((10, K), dtype)  # lumped sums without mu/r
    w = np.ones(K, dtype)  # (a/r)^n
    for n in range(N + 1):
        P, P1, P2 = bufs[n % 3], bufs[(n - 1) % 3], bufs[(n - 2) % 3]
        if n == 0:
            P[:, 0] = 1.0
        else:
            P[:, :n] = g[n, :n] * (ucol * P1[:, :n]) - h[n, :n] * P2[:, :n]
            P[:, n] = f[n] * c * P1[:, n - 1]
        e = n + 1
        P0, Pa, Pb = P[:, :e], P[:, 1 : e + 1], P[:, 2 : e + 2]
        kn = k[n, :e]
        dP = kn * Pa - mt[:, :e] * P0
        d2P = d2c[:, :e] * P0 - (kn * k2t[:, :e]) * Pa + j[n, :e] * Pb
        cm, sm = cosm[:, :e], sinm[:, :e]
        Wn = W[n]
        s0 = np.concatenate([P0 * cm, P0 * sm], axis=1) @ Wn  # PX, mPY, m2PX
        s1 = np.concatenate([dP * cm, dP * sm], axis=1) @ Wn[:, :2]  # dPX, mdPY
        s2 = np.concatenate([d2P * cm, d2P * sm], axis=1) @ Wn[:, 0]  # d2PX
        PX, mPY, m2PX = s0[:, 0], s0[:, 1], s0[:, 2]
        dPX, mdPY = s1[:, 0], s1[:, 1]
        n1, n12 = n + 1.0, (n + 1.0) * (n + 2.0)
        acc[0] += w * PX
        acc[1] += w * mPY
        acc[2] += w * dPX
        acc[3] -= n1 * w * PX
        acc[4] -= w * m2PX
        acc[5] += w * mdPY
        acc[6] += w * s2
        acc[7] -= n1 * w * mPY
        acc[8] -= n1 * w * dPX
        acc[9] += n12 * w * PX
        w = w * q
    base = mu / r
    D = acc * base
    D[[3, 7, 8]] /= r
    D[9] /= r * r
    return spherical_to_local(D, r, phi)


def legendre_a(model, lam, phi, r, dtype=np.float64, chunk=None):
    """LEGENDRE_a: FCR on fnALFs, all points of a block processed together (Section 3)."""
    lam, phi, r = (_as(x, dtype) for x in np.broadcast_arrays(lam, phi, r))
    chunk = chunk or DEFAULT_CHUNK["LEGENDRE_a"]
    out = np.empty((10, lam.size), dtype)
    for sl in _chunks(lam.size, chunk):
        out[:, sl] = _legendre_block(model, lam[sl], phi[sl], r[sl], dtype)
    return out


# ============================================================
# 1b. LEGENDRE_b -- same algorithm, point by point (interpreted driver + compiled kernels)
# ============================================================
def _legendre_b_point(lam, phi, r, C, S, T, N, mu, a, Mrow, Ncol):
    f, g, h = T["f"], T["g"], T["h"]
    Kt, J = T["k"][: N + 1, : N + 1], T["j"][: N + 1, : N + 1]
    u, c = np.sin(phi), np.cos(phi)
    t = u / c
    P = np.zeros((N + 1, N + 3), C.dtype)
    P[0, 0] = 1.0
    if N >= 1:
        P[1, 0] = g[1, 0] * u
        P[1, 1] = f[1] * c
    for n in range(2, N + 1):
        P[n, :n] = g[n, :n] * u * P[n - 1, :n] - h[n, :n] * P[n - 2, :n]
        P[n, n] = f[n] * c * P[n - 1, n - 1]
    P0, Pa, Pb = P[:, : N + 1], P[:, 1 : N + 2], P[:, 2 : N + 3]
    dP = Kt * Pa - Mrow * t * P0
    d2P = (Mrow * Mrow * t * t - Mrow / (c * c)) * P0 - (2 * Mrow + 1) * Kt * t * Pa + J * Pb
    cm, sm = np.cos(Mrow * lam), np.sin(Mrow * lam)
    X = C * cm + S * sm
    Y = S * cm - C * sm
    PX = (P0 * X).sum(1)
    mPY = (P0 * Mrow * Y).sum(1)
    m2PX = (P0 * Mrow * Mrow * X).sum(1)
    dPX = (dP * X).sum(1)
    mdPY = (dP * Mrow * Y).sum(1)
    d2PX = (d2P * X).sum(1)
    w = (a / r) ** Ncol
    n1 = Ncol + 1
    base = mu / r
    D = np.array([
        (w * PX).sum(), (w * mPY).sum(), (w * dPX).sum(), -(n1 * w * PX).sum() / r,
        -(w * m2PX).sum(), (w * mdPY).sum(), (w * d2PX).sum(), -(n1 * w * mPY).sum() / r,
        -(n1 * w * dPX).sum() / r, (n1 * (n1 + 1) * w * PX).sum() / (r * r)]) * base
    return D


def legendre_b(model, lam, phi, r, dtype=np.float64, chunk=None):
    """LEGENDRE_b: FCR on fnALFs, one call per observation point (Matlab + DLL style)."""
    lam, phi, r = (_as(x, dtype) for x in np.broadcast_arrays(lam, phi, r))
    N = model.N
    T = tables(N, dtype)
    C, S = model.coef(dtype)
    Mrow = np.arange(N + 1).astype(dtype)
    Ncol = Mrow.copy()
    D = np.empty((10, lam.size), dtype)
    for i in range(lam.size):
        D[:, i] = _legendre_b_point(lam[i], phi[i], r[i], C, S, T, N, model.gm, model.a,
                                    Mrow, Ncol)
    return spherical_to_local(D, r, phi)


# ============================================================
# 2. PINES -- FCR on fully normalised Helmholtz polynomials, Cartesian
# ============================================================
def helmholtz_fcr(u, N, dtype=np.float64):
    """Fully normalised Helmholtz polynomials Abar_nm(u), (K, N+1, N+1), via the FCR."""
    T = tables(N, dtype)
    u = _as(u, dtype)
    A = np.zeros((u.size, N + 1, N + 1), dtype)
    A[:, 0, 0] = 1.0
    for n in range(1, N + 1):
        A[:, n, :n] = T["g"][n, :n] * u[:, None] * A[:, n - 1, :n] - T["h"][n, :n] * A[:, n - 2, :n]
        A[:, n, n] = T["f"][n] * A[:, n - 1, n - 1]
    return A


def helmholtz_msr(u, N, dtype=np.float64):
    """
    Fully normalised Helmholtz polynomials via Pines' original *modified sub-sectorial
    recursion* (MSR): sectorial seed, then downward in m along each row,
        (n-m) A_nm = u A_n,m+1 - A_n-1,m+1     (unnormalised)
    which is numerically unstable for high degree (Fig. 3 of the paper).
    """
    T = tables(N, dtype)
    u = _as(u, dtype)
    k, s2, f = T["k"], T["msr2"], T["f"]
    A = np.zeros((u.size, N + 1, N + 2), dtype)
    A[:, 0, 0] = 1.0
    for n in range(1, N + 1):
        A[:, n, n] = f[n] * A[:, n - 1, n - 1]
        for m in range(n - 1, -1, -1):
            A[:, n, m] = (k[n, m] * u * A[:, n, m + 1] - s2[n, m] * A[:, n - 1, m + 1]) / (n - m)
    return A[:, :, : N + 1]


def legendre_fcr(phi, N, dtype=np.float64):
    """fnALFs Pbar_nm(sin phi), (K, N+1, N+1), via the standard FCR (as in LEGENDRE_a)."""
    T = tables(N, dtype)
    phi = _as(phi, dtype)
    u, c = np.sin(phi), np.cos(phi)
    P = np.zeros((phi.size, N + 1, N + 1), dtype)
    P[:, 0, 0] = 1.0
    for n in range(1, N + 1):
        P[:, n, :n] = T["g"][n, :n] * u[:, None] * P[:, n - 1, :n] - T["h"][n, :n] * P[:, n - 2, :n]
        P[:, n, n] = T["f"][n] * c * P[:, n - 1, n - 1]
    return P


def dlegendre_dphi(P, phi, dtype=np.float64):
    """dPbar_nm/dphi = k_nm Pbar_n,m+1 - m tan(phi) Pbar_nm from a (K, N+1, N+1) triangle."""
    K, N1, _ = P.shape
    T = tables(N1 - 1, dtype)
    t = np.tan(_as(phi, dtype))
    Pa = np.zeros_like(P)
    Pa[:, :, :-1] = P[:, :, 1:]
    m = np.arange(N1).astype(dtype)
    return T["k"][None, :N1, :N1] * Pa - m[None, None, :] * t[:, None, None] * P


def _pines_coef(model, dtype):
    key = ("pines", np.dtype(dtype).name)
    if key not in model._cache:
        C, S = model.coef(dtype)
        T = tables(model.N, dtype)
        out = []
        for n in range(model.N + 1):
            e = n + 1
            m = np.arange(e).astype(dtype)
            c, s = C[n, :e], S[n, :e]
            Rv = np.stack([np.ones(e, dtype), T["k"][n, :e], -T["R11"][n, :e], T["j"][n, :e],
                           -T["R12"][n, :e], T["R22"][n, :e]])
            R1 = Rv[:3] * m
            mm1 = m * (m - 1)
            out.append((
                np.ascontiguousarray((Rv * c)[..., None]), np.ascontiguousarray((Rv * s)[..., None]),
                np.ascontiguousarray(np.stack([R1 * c, R1 * s], -1)),
                np.ascontiguousarray(np.stack([R1 * s, -R1 * c], -1)),
                np.ascontiguousarray(np.stack([mm1 * c, mm1 * s], -1)),
                np.ascontiguousarray(np.stack([mm1 * s, -mm1 * c], -1)),
            ))
        model._cache[key] = out
    return model._cache[key]


def _pines_block(model, lam, phi, r, dtype):
    N, mu = model.N, model.gm
    T = tables(N, dtype)
    W = _pines_coef(model, dtype)
    f, g, h = T["f"], T["g"], T["h"]
    K = lam.size
    cp = np.cos(phi)
    s, t, u = cp * np.cos(lam), cp * np.sin(lam), np.sin(phi)  # direction cosines
    q = model.a / r
    # r_m + i i_m = (s + i t)^m, stored at column m+2 (columns 0, 1 = m-2, m-1 < 0 -> 0)
    RM = np.zeros((K, N + 3), dtype)
    IM = np.zeros((K, N + 3), dtype)
    RM[:, 2] = 1.0
    for m in range(1, N + 1):
        RM[:, m + 2] = s * RM[:, m + 1] - t * IM[:, m + 1]
        IM[:, m + 2] = s * IM[:, m + 1] + t * RM[:, m + 1]
    # Helmholtz rows n, n+1, n+2 in 4 rolling buffers
    L = N + 5
    A = [np.zeros((K, L), dtype) for _ in range(4)]
    A[0][:, 0] = 1.0
    ucol = u[:, None]

    def hrow(n):
        An, A1, A2 = A[n % 4], A[(n - 1) % 4], A[(n - 2) % 4]
        An[:, :n] = g[n, :n] * (ucol * A1[:, :n]) - h[n, :n] * A2[:, :n]
        An[:, n] = f[n] * A1[:, n - 1]

    hrow(1)
    hrow(2)
    acc = np.zeros((14, K), dtype)
    w = np.ones(K, dtype)
    for n in range(N + 1):
        e = n + 1
        A0, B1, B2 = A[n % 4], A[(n + 1) % 4], A[(n + 2) % 4]
        var = np.stack([A0[:, :e], A0[:, 1 : e + 1], B1[:, 1 : e + 1], A0[:, 2 : e + 2],
                        B1[:, 2 : e + 2], B2[:, 2 : e + 2]])  # (6, K, e)
        W0r, W0i, W1r, W1i, W2r, W2i = W[n]
        rm, im = RM[:, 2 : e + 2], IM[:, 2 : e + 2]
        rm1, im1 = RM[:, 1 : e + 1], IM[:, 1 : e + 1]
        rm2, im2 = RM[:, :e], IM[:, :e]
        o6 = (var * rm) @ W0r + (var * im) @ W0i  # V, a3, a4, a33, a34, a44
        v3 = var[:3]
        o3 = (v3 * rm1) @ W1r + (v3 * im1) @ W1i  # (a1, a2), (a13, a23), (a14, a24)
        o1 = (A0[:, :e] * rm2) @ W2r + (A0[:, :e] * im2) @ W2i  # a11, a12
        acc[0:6] += w * o6[:, :, 0]
        acc[6:12] += w * np.moveaxis(o3, 2, 1).reshape(6, K)
        acc[12:14] += w * o1.T
        w = w * q
        if n + 3 <= N + 2:
            hrow(n + 3)
    r1 = mu / (r * r)
    r2 = r1 / r
    V = acc[0] * (mu / r)
    a3, a4 = acc[1] * r1, acc[2] * r1
    a33, a34, a44 = acc[3] * r2, acc[4] * r2, acc[5] * r2
    a1, a2 = acc[6] * r1, acc[7] * r1
    a13, a23, a14, a24 = acc[8] * r2, acc[9] * r2, acc[10] * r2, acc[11] * r2
    a11, a12 = acc[12] * r2, acc[13] * r2
    xh = (s, t, u)
    ai4 = (a14, a24, a34)
    aij = ((a11, a12, a13), (a12, -a11, a23), (a13, a23, a33))
    a4r = a4 / r
    out = np.empty((10, K), dtype)
    out[0] = V
    out[1], out[2], out[3] = a1 + s * a4, a2 + t * a4, a3 + u * a4
    col = 4
    for i_, j_ in ((0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2)):
        # a44 built from Abar_n+2,m+2 is d(a4)/dr - (1/r) sum x_k d(a4)/ds_k + a4/r, i.e. the
        # combination (a44_Pines - a4/r) of the chain rule already (identity
        # (n+m+2) A_n+1,m+1 + u A_n+1,m+2 = A_n+2,m+2 - A_n+1,m+1)
        Vij = aij[i_][j_] + xh[j_] * ai4[i_] + xh[i_] * ai4[j_] + xh[i_] * xh[j_] * a44
        if i_ == j_:
            Vij = Vij + a4r
        out[col] = Vij
        col += 1
    return out


def pines(model, lam, phi, r, dtype=np.float64, chunk=None, cartesian=False):
    """PINES: FCR on fully normalised Helmholtz polynomials, Cartesian output (Section 4)."""
    lam, phi, r = (_as(x, dtype) for x in np.broadcast_arrays(lam, phi, r))
    chunk = chunk or DEFAULT_CHUNK["PINES"]
    cart = np.empty((10, lam.size), dtype)
    for sl in _chunks(lam.size, chunk):
        cart[:, sl] = _pines_block(model, lam[sl], phi[sl], r[sl], dtype)
    loc = cart_to_local(cart, lam, phi)
    return (loc, cart) if cartesian else loc


# ============================================================
# 3. METRIS -- solid spherical harmonics, derivatives acting on the coefficients
# ============================================================
def _metris_op(E, kind, T):
    """Apply d/dx, d/dy or d/dz (in units of 1/a) to a functional Re sum E_nm Ybar_nm."""
    L = E.shape[0]
    al, be, ga = T["alpha"][:L, :L], T["beta"][:L, :L], T["gamma"][:L, :L]
    out = np.zeros_like(E)
    src = E[:-1]
    if kind == "z":
        out[1:, :] -= al[:-1, :] * src
    else:
        half = 0.5 if kind == "x" else 0.5j
        sgn = -1.0 if kind == "x" else 1.0
        out[1:, 1:] += (sgn * half) * be[:-1, :-1] * src[:, :-1]  # m -> m+1
        out[1:, :-1] += half * ga[:-1, 1:] * src[:, 1:]  # m -> m-1 (m >= 1)
        out[1:, 1] += (sgn * half) * be[:-1, 0] * np.conj(src[:, 0])  # m = 0 -> 1, real part
    return out


def _metris_coef(model, dtype):
    key = ("metris", np.dtype(dtype).name)
    if key not in model._cache:
        N, a = model.N, model.a
        C, S = model.coef(dtype)
        T = tables(N, dtype)
        L = N + 3
        ctype = np.result_type(dtype, np.complex64)
        EK = np.zeros((L, L), ctype)
        EK[: N + 1, : N + 1] = (model.gm / a) * (C - 1j * S)
        Ex, Ey, Ez = (_metris_op(EK, k, T) / a for k in "xyz")
        Es = [EK, Ex, Ey, Ez, _metris_op(Ex, "x", T) / a, _metris_op(Ex, "y", T) / a,
              _metris_op(Ex, "z", T) / a, _metris_op(Ey, "y", T) / a,
              _metris_op(Ey, "z", T) / a, _metris_op(Ez, "z", T) / a]
        E = np.stack(Es, axis=-1)  # (L, L, 10)
        model._cache[key] = (np.ascontiguousarray(E.real), np.ascontiguousarray(E.imag))
    return model._cache[key]


def _metris_block(model, lam, phi, r, dtype):
    N = model.N
    T = tables(N, dtype)
    Er, Ei = _metris_coef(model, dtype)
    f, g, h = T["f"], T["g"], T["h"]
    K = lam.size
    x, y, z = sph_to_cart(lam, phi, r)
    a = model.a
    xp, yp, zp, rp = x / a, y / a, z / a, r / a
    r2i = 1.0 / (rp * rp)
    zr, zeta_r, zeta_i = zp * r2i, xp * r2i, yp * r2i
    zcol, r2col = zr[:, None], r2i[:, None]
    L = N + 3
    Yr = [np.zeros((K, L), dtype) for _ in range(3)]
    Yi = [np.zeros((K, L), dtype) for _ in range(3)]
    F = np.zeros((K, 10), dtype)
    for n in range(N + 3):
        R0, R1, R2 = Yr[n % 3], Yr[(n - 1) % 3], Yr[(n - 2) % 3]
        I0, I1, I2 = Yi[n % 3], Yi[(n - 1) % 3], Yi[(n - 2) % 3]
        if n == 0:
            R0[:, 0] = 1.0 / rp
        else:
            R0[:, :n] = g[n, :n] * (zcol * R1[:, :n]) - h[n, :n] * (r2col * R2[:, :n])
            I0[:, :n] = g[n, :n] * (zcol * I1[:, :n]) - h[n, :n] * (r2col * I2[:, :n])
            ar, ai = R1[:, n - 1], I1[:, n - 1]
            R0[:, n] = f[n] * (zeta_r * ar - zeta_i * ai)
            I0[:, n] = f[n] * (zeta_r * ai + zeta_i * ar)
        e = n + 1
        F += R0[:, :e] @ Er[n, :e] - I0[:, :e] @ Ei[n, :e]
    return F.T


def metris(model, lam, phi, r, dtype=np.float64, chunk=None, cartesian=False):
    """METRIS: recursion on fnSSHs (Section 5), Cartesian output."""
    lam, phi, r = (_as(x, dtype) for x in np.broadcast_arrays(lam, phi, r))
    chunk = chunk or DEFAULT_CHUNK["METRIS"]
    cart = np.empty((10, lam.size), dtype)
    for sl in _chunks(lam.size, chunk):
        cart[:, sl] = _metris_block(model, lam[sl], phi[sl], r[sl], dtype)
    loc = cart_to_local(cart, lam, phi)
    return (loc, cart) if cartesian else loc


# ============================================================
# 4. CLENSHAW -- Clenshaw summation over the degree (Tables 4-6)
# ============================================================
def _clenshaw_vterms(model, u, c, q, dtype):
    """
    v-terms of Table 5 for K points (or K parallels): (6, 2, K, N+1) for
    [v, v_u, v_r, v_uu, v_ur, v_rr], second axis = cosine (C) / sine (S) part.
    """
    N = model.N
    T = tables(N, dtype)
    C, S = model.coef(dtype)
    g, h, f = T["g"], T["h"], T["f"]
    K = u.size
    coef = np.stack([C, S])  # (2, N+1, N+1)
    uq, q2 = (u * q)[:, None], (q * q)[:, None]
    qc = q[:, None]
    y1 = np.zeros((6, 2, K, N + 1), dtype)  # y_{n+1}
    y2 = np.zeros((6, 2, K, N + 1), dtype)  # y_{n+2}
    yfin = np.empty((6, 2, K, N + 1), dtype)
    for n in range(N, -1, -1):
        e = n + 1
        b1, b2 = y1[..., :e], y2[..., :e]
        gn = g[n + 1, :e]
        Yn = (uq * gn) * b1 - (q2 * h[n + 2, :e]) * b2
        cn = coef[:, n, :e][:, None, :]
        Yn[0] += cn
        Yn[1] += (n + 1.0) * cn
        Yn[2] += (n + 1.0) * (n + 2.0) * cn
        qg = qc * gn
        Yn[3] += qg * b1[0]
        Yn[4] += (2.0 * qg) * b1[3]
        Yn[5] += qg * b1[1]
        yfin[..., n] = Yn[..., n]
        y2, y1 = b1, Yn
    y0, yr, yrr, yu, yuu, yur = yfin
    mm = np.arange(N + 1).astype(dtype)
    cq = (c * q)[:, None]
    pm = np.empty((K, N + 1), dtype)
    pm[:, 0] = 1.0
    pm[:, 1:] = np.cumprod(f[1 : N + 1][None, :] * cq, axis=1)  # Pbar_mm (a/r)^m
    c2 = (c * c)[:, None]
    uc = u[:, None]
    dp = (-mm * uc / c2) * pm
    d2p = (-mm / c2 + mm * (mm - 2) * uc * uc / (c2 * c2)) * pm
    v = y0 * pm
    vu = y0 * dp + yu * pm
    vr = yr * pm
    vuu = yuu * pm + y0 * d2p + 2.0 * yu * dp
    vur = yr * dp + yur * pm
    vrr = yrr * pm
    return np.stack([v, vu, vr, vuu, vur, vrr])


def _clenshaw_assemble(Sc, Ss1, Sc2, r, phi, mu):
    """Table 6. Sc: (6, ...) sums of x1 cos + x2 sin; Ss1: (3, ...) sums of m(x1 sin - x2 cos)
    for v, v_u, v_r; Sc2: sum of m^2 (v1 cos + v2 sin)."""
    c, u = np.cos(phi), np.sin(phi)
    b = mu / r
    D = np.empty((10,) + np.shape(Sc[0]), Sc.dtype)
    D[0] = b * Sc[0]
    D[1] = -b * Ss1[0]
    D[2] = b * c * Sc[1]
    D[3] = -b / r * Sc[2]
    D[4] = -b * Sc2
    D[5] = -b * c * Ss1[1]
    D[6] = b * c * c * Sc[3] - b * u * Sc[1]
    D[7] = b / r * Ss1[2]
    D[8] = -b * c / r * Sc[4]
    D[9] = b / (r * r) * Sc[5]
    return D


def _clenshaw_block(model, lam, phi, r, dtype):
    N = model.N
    u, c = np.sin(phi), np.cos(phi)
    q = model.a / r
    vt = _clenshaw_vterms(model, u, c, q, dtype)  # (6, 2, K, N+1)
    mm = np.arange(N + 1).astype(dtype)
    ang = lam[:, None] * mm[None, :]
    cm, sm = np.cos(ang), np.sin(ang)
    Pc = vt[:, 0] * cm + vt[:, 1] * sm  # (6, K, N+1)
    Ps = vt[:3, 0] * sm - vt[:3, 1] * cm
    Sc = Pc.sum(-1)
    Ss1 = Ps @ mm
    Sc2 = Pc[0] @ (mm * mm)
    D = _clenshaw_assemble(Sc, Ss1, Sc2, r, phi, model.gm)
    return spherical_to_local(D, r, phi)


def clenshaw(model, lam, phi, r, dtype=np.float64, chunk=None):
    """CLENSHAW for scattered points: vectorised over points and orders."""
    lam, phi, r = (_as(x, dtype) for x in np.broadcast_arrays(lam, phi, r))
    chunk = chunk or DEFAULT_CHUNK["CLENSHAW"]
    out = np.empty((10, lam.size), dtype)
    for sl in _chunks(lam.size, chunk):
        out[:, sl] = _clenshaw_block(model, lam[sl], phi[sl], r[sl], dtype)
    return out


def clenshaw_grid(model, lam1d, phi1d, r, dtype=np.float64):
    """
    CLENSHAW on a grid (constant r): Clenshaw recursions once per parallel, then the sums over
    the orders for all meridians as (nlat, N+1) @ (N+1, nlon) products (Section 6).
    Returns (10, nlat * nlon), latitude-major (same order as grid_points()).
    """
    lam1d, phi1d = _as(lam1d, dtype), _as(phi1d, dtype)
    r = dtype(r) if np.ndim(r) == 0 else _as(r, dtype)[0]
    N = model.N
    nlat, nlon = phi1d.size, lam1d.size
    u, c = np.sin(phi1d), np.cos(phi1d)
    q = np.full(nlat, model.a / r, dtype)
    vt = _clenshaw_vterms(model, u, c, q, dtype)  # (6, 2, nlat, N+1)
    mm = np.arange(N + 1).astype(dtype)
    ang = mm[:, None] * lam1d[None, :]
    cm, sm = np.cos(ang), np.sin(ang)  # (N+1, nlon)
    Sc = vt[:, 0] @ cm + vt[:, 1] @ sm  # (6, nlat, nlon)
    vm = vt[:3] * mm
    Ss1 = vm[:, 0] @ sm - vm[:, 1] @ cm  # (3, nlat, nlon)
    Sc2 = (vt[0, 0] * mm * mm) @ cm + (vt[0, 1] * mm * mm) @ sm
    PHI = np.repeat(phi1d, nlon)
    R = np.full(nlat * nlon, r, dtype)
    D = _clenshaw_assemble(Sc.reshape(6, -1), Ss1.reshape(3, -1), Sc2.reshape(-1), R, PHI,
                           model.gm)
    return spherical_to_local(D, R, PHI)


# ============================================================
# DISPATCH
# ============================================================
_FUNCS = {"LEGENDRE_a": legendre_a, "LEGENDRE_b": legendre_b, "PINES": pines,
          "METRIS": metris, "CLENSHAW": clenshaw}


def synthesize(method, model, lam, phi, r, dtype=np.float64, chunk=None, cartesian=False):
    """Local quantities (10, K) of `method`; with cartesian=True also the Cartesian (10, K)
    (native for PINES/METRIS, rotated from the local frame for the spherical codes)."""
    fn = _FUNCS[method]
    if method in ("PINES", "METRIS"):
        return fn(model, lam, phi, r, dtype=dtype, chunk=chunk, cartesian=cartesian)
    loc = fn(model, lam, phi, r, dtype=dtype, chunk=chunk)
    if not cartesian:
        return loc
    lam_, phi_ = (_as(x, dtype) for x in np.broadcast_arrays(lam, phi, r)[:2])
    return loc, local_to_cart(loc, lam_, phi_)


def local_to_cart(loc, lam, phi):
    """Inverse of cart_to_local."""
    R = local_frame(lam, phi)
    gl = np.moveaxis(loc[1:4], 0, -1)
    g = np.einsum("kji,kj->ki", R, gl)
    ll, lp, pp, lr, pr, rr = loc[4], loc[5], loc[6], loc[7], loc[8], loc[9]
    Tl = np.stack([np.stack([ll, lp, lr], -1), np.stack([lp, pp, pr], -1),
                   np.stack([lr, pr, rr], -1)], -2)
    T = np.einsum("kai,kab,kbj->kij", R, Tl, R)
    return np.stack([loc[0], g[:, 0], g[:, 1], g[:, 2], T[:, 0, 0], T[:, 0, 1], T[:, 0, 2],
                     T[:, 1, 1], T[:, 1, 2], T[:, 2, 2]])


# ============================================================
# OBSERVATION SCENARIOS (paper p. 12)
# ============================================================
def grid_axes(step=1.0):
    """GRID: lambda = -179..180 deg, phi = -89..89 deg (step 1 deg -> 360 x 179 = 64,440)."""
    lon = np.arange(-179.0, 180.0 + 1e-9, step)
    lat = np.arange(-89.0, 89.0 + 1e-9, step)
    return lon, lat


def grid_points(h=250e3, lon=None, lat=None, a=A_EARTH):
    """Flattened latitude-major GRID points (lam, phi [rad], r [m]) at height h above a."""
    if lon is None or lat is None:
        lon0, lat0 = grid_axes()
        lon = lon0 if lon is None else lon
        lat = lat0 if lat is None else lat
    LAT, LON = np.meshgrid(np.radians(lat), np.radians(lon), indexing="ij")
    r = np.full(LAT.size, a + h)
    return LON.ravel(), LAT.ravel(), r


def orbit_points(n_points=64440, dt=31.0, sma=6628e3, ecc=0.002, inc=96.5, raan0=-125.0,
                 argp0=0.0, M0=0.0, gm=GM_EARTH, J2=J2_EARTH, Re=A_EARTH):
    """
    ORBIT1/ORBIT2: near-circular, near-polar orbit with the secular J2 drift of the node,
    perigee and mean anomaly; positions in the body-fixed frame without Earth rotation (as in
    Fig. 4 of the paper). Returns lam, phi [rad], r [m], t [s].
    """
    t = dt * np.arange(n_points)
    i = np.radians(inc)
    n0 = np.sqrt(gm / sma**3)
    p = sma * (1 - ecc**2)
    fac = 1.5 * J2 * (Re / p) ** 2
    raan = np.radians(raan0) - fac * n0 * np.cos(i) * t
    argp = np.radians(argp0) + 0.5 * fac * n0 * (5 * np.cos(i) ** 2 - 1) * t
    M = np.radians(M0) + n0 * (1 + 0.5 * fac * np.sqrt(1 - ecc**2) * (3 * np.cos(i) ** 2 - 1)) * t
    E = M.copy()
    for _ in range(8):
        E -= (E - ecc * np.sin(E) - M) / (1 - ecc * np.cos(E))
    nu = 2 * np.arctan2(np.sqrt(1 + ecc) * np.sin(E / 2), np.sqrt(1 - ecc) * np.cos(E / 2))
    rr = sma * (1 - ecc * np.cos(E))
    th = argp + nu
    x = rr * (np.cos(raan) * np.cos(th) - np.sin(raan) * np.sin(th) * np.cos(i))
    y = rr * (np.sin(raan) * np.cos(th) + np.cos(raan) * np.sin(th) * np.cos(i))
    z = rr * np.sin(th) * np.sin(i)
    return np.arctan2(y, x), np.arcsin(z / rr), rr, t


# ============================================================
# HELPERS
# ============================================================
def digits(x, ref, cap=16.0):
    """Number of coinciding decimal digits, -log10(|x - ref| / |ref|), capped at `cap`."""
    x = np.asarray(x, dtype=np.longdouble)
    ref = np.asarray(ref, dtype=np.longdouble)
    with np.errstate(divide="ignore", invalid="ignore"):
        rel = np.abs(x - ref) / np.abs(ref)
        d = -np.log10(rel.astype(np.float64))
    d = np.where(np.isfinite(d), d, cap)
    return np.minimum(d, cap)


def timed(fn, *args, **kw):
    """(result, wall seconds) of fn(*args, **kw) (process CPU time is also reported)."""
    t0, c0 = time.perf_counter(), time.process_time()
    out = fn(*args, **kw)
    return out, time.perf_counter() - t0, time.process_time() - c0


# ============================================================
# SELF-TEST
# ============================================================
def _selftest(N=30, K=40, seed=1):
    rng = np.random.default_rng(seed)
    model = GravityModel.synthetic(N, seed=7)
    lam = rng.uniform(-np.pi, np.pi, K)
    phi = np.radians(rng.uniform(-85, 85, K))
    r = A_EARTH + rng.uniform(200e3, 500e3, K)
    res = {}
    for meth in METHODS:
        res[meth] = synthesize(meth, model, lam, phi, r)
    ref = clenshaw(model, lam, phi, r, dtype=np.longdouble)
    print(f"self-test N={N}, K={K}: min coinciding digits vs CLENSHAW*16")
    print("            " + " ".join(f"{q:>5s}" for q in QUANTS))
    for meth in METHODS:
        d = digits(res[meth], ref).min(axis=1)
        print(f"{meth:11s} " + " ".join(f"{v:5.1f}" for v in d))
    lap = {m: np.abs(res[m][4] + res[m][6] + res[m][9]).max() for m in METHODS}
    print("max |Laplace| [s^-2]: " + ", ".join(f"{m} {v:.1e}" for m, v in lap.items()))
    # grid mode = scattered mode
    lon, lat = np.array([-120.0, 10.0, 77.0]), np.array([-60.0, 5.0, 33.0, 80.0])
    gl, gp, gr = grid_points(300e3, lon, lat)
    a_ = clenshaw_grid(model, np.radians(lon), np.radians(lat), gr[0])
    b_ = legendre_a(model, gl, gp, gr)
    print(f"CLENSHAW grid vs LEGENDRE_a min digits: {digits(a_, b_).min():.1f}")
    # Cartesian cross-check PINES vs METRIS
    _, cp_ = pines(model, lam, phi, r, cartesian=True)
    _, cm_ = metris(model, lam, phi, r, cartesian=True)
    print(f"PINES vs METRIS (Cartesian) min digits: {digits(cp_, cm_).min(axis=1).round(1)}")
    # finite-difference check of the radial derivatives
    dr = 1.0
    vp = legendre_a(model, lam, phi, r + dr)
    vm = legendre_a(model, lam, phi, r - dr)
    fd = (vp[3] - vm[3]) / (2 * dr)
    print(f"FD check V_rr: max rel err {np.max(np.abs(fd - res['LEGENDRE_a'][9]) / np.abs(fd)):.1e}")


if __name__ == "__main__":
    _selftest()
