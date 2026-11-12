"""
Journal figure style for the spherical-harmonic synthesis comparison (Fantino & Casotto 2009).

Copied from MCMCPlots.py (FIGURE STYLE / FIGURE OUTPUT sections) so that the SPH scripts do
not import the MCMC stack (corner, statsmodels, GMM). Same conventions:

    * printed size: single column COL_W, full text width TXT_W [in]
    * LaTeX text at 7-8 pt (mathtext + Computer Modern fallback if latex/dvipng are missing)
    * Okabe-Ito colors, one meaning per color in every figure
    * panel titles "(a) ..." on the left, axis labels "Name [unit]"
    * figure legends above the panels (fig.legend(loc="outside upper center")),
      layout="constrained"
    * one vector PDF per figure, written atomically (robust to a slow OneDrive folder)

    import sph_style as ss
    ss.apply_style()
    ss.set_figure_dir("results")
    fig, ax = plt.subplots(figsize=(ss.COL_W, 2.4), layout="constrained")
    ...
    ss.save_figure(fig, "fig_name")
"""

import io
import os
import shutil
import time
from pathlib import Path

import matplotlib

if os.environ.get("MPLBACKEND") is None and not os.environ.get("DISPLAY") and os.name != "nt":
    # Batch runs (python run_all.py) need no GUI; interactive sessions can still override.
    try:
        import sys

        if not hasattr(sys, "ps1") and "ipykernel" not in sys.modules:
            matplotlib.use("Agg")
    except Exception:
        pass

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np


def _patch_type1_subset():
    r"""
    Keep encoding slot 0 of the usetex Type-1 fonts embedded in PDFs (matplotlib 3.11).

    matplotlib 3.11 subsets the TeX Type-1 fonts it embeds, and Type1Font.subset() overwrites
    encoding slot 0 with .notdef: the minus sign (cmsy slot 0) and \Gamma (cmr slot 0) vanish
    from every usetex PDF. When a font really uses slot 0 it is embedded whole; every other
    font is still subset. Idempotent, and a no-op on versions without Type1Font.subset.
    """
    try:
        from matplotlib._type1font import Type1Font
    except ImportError:
        return
    orig = getattr(Type1Font, "subset", None)
    if orig is None or getattr(orig, "_jf_keeps_slot0", False):
        return

    def subset(self, characters, name_prefix):
        characters = frozenset(characters)
        if 0 in characters and self.prop["Encoding"].get(0, ".notdef") != ".notdef":
            return self
        return orig(self, characters, name_prefix)

    subset._jf_keeps_slot0 = True
    Type1Font.subset = subset


_patch_type1_subset()

# ============================================================
# FIGURE STYLE (journal)
# ============================================================
COL_W = 3.5  # single-column width [in]
TXT_W = 7.16  # full text width (double column) [in]
SHOW_FIGS = False  # plt.show() after each figure is saved (interactive backends only)
FALLBACK_DIR = Path.home() / "figures_fallback"


def _ensure_tex_on_path():
    """
    True if latex + dvipng are reachable. GUI-launched processes (VS Code, Jupyter, Spyder on
    macOS) often do not inherit the login-shell PATH, so the standard TeX locations are added
    before giving up; otherwise matplotlib would silently fall back to mathtext and DejaVu.
    """
    if shutil.which("latex") and shutil.which("dvipng"):
        return True
    import glob

    extra = ["/Library/TeX/texbin", "/opt/homebrew/bin", "/usr/local/bin", "/usr/texbin"]
    extra += sorted(glob.glob("/usr/local/texlive/*/bin/*"), reverse=True)
    path = os.environ.get("PATH", "").split(os.pathsep)
    os.environ["PATH"] = os.pathsep.join(path + [d for d in extra if os.path.isdir(d)])
    return bool(shutil.which("latex") and shutil.which("dvipng"))


USE_TEX = _ensure_tex_on_path() and os.environ.get("SPH_NO_TEX", "0") != "1"
if not USE_TEX:
    print("[figures] latex/dvipng not used: falling back to mathtext (Computer Modern)")
EN_DASH = "--" if USE_TEX else "-"

STYLE = {
    "text.usetex": USE_TEX,
    "text.latex.preamble": r"\usepackage{amsmath}\usepackage{amssymb}\usepackage{mathrsfs}",
    "mathtext.fontset": "cm",
    "axes.formatter.use_mathtext": True,
    "font.family": "serif",
    "font.serif": ["Computer Modern Roman", "cmr10", "DejaVu Serif"],
    "font.size": 8,
    "axes.labelsize": 8,
    "axes.titlesize": 8,
    "legend.fontsize": 7,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "axes.linewidth": 0.6,
    "axes.edgecolor": "0.15",
    "axes.labelcolor": "0.1",
    "axes.titlepad": 4.0,
    "axes.labelpad": 2.5,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "axes.axisbelow": True,
    "grid.color": "0.9",
    "grid.linestyle": "-",
    "grid.linewidth": 0.4,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.color": "0.15",
    "ytick.color": "0.15",
    "xtick.major.size": 3.0,
    "ytick.major.size": 3.0,
    "xtick.minor.size": 1.6,
    "ytick.minor.size": 1.6,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.minor.width": 0.4,
    "ytick.minor.width": 0.4,
    "lines.linewidth": 1.2,
    "lines.markersize": 4.0,
    "legend.frameon": False,
    "legend.handlelength": 1.8,
    "legend.handletextpad": 0.5,
    "legend.columnspacing": 1.2,
    "legend.borderaxespad": 0.3,
    "figure.dpi": 150,
    "savefig.dpi": 600,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.02,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
}


def apply_style():
    """Set the journal rcParams (call before drawing any figure)."""
    plt.rcParams.update(STYLE)


# Okabe-Ito colors (colorblind-safe). Each method has ONE color in every figure.
METHODS = ("LEGENDRE_a", "LEGENDRE_b", "PINES", "METRIS", "CLENSHAW")
METHOD_COLOR = {
    "LEGENDRE_a": "#0072B2",  # blue
    "LEGENDRE_b": "#56B4E9",  # sky blue (same algorithm, interpreted implementation)
    "PINES": "#009E73",  # bluish green
    "METRIS": "#D55E00",  # vermillion
    "CLENSHAW": "#CC79A7",  # reddish purple
}
METHOD_MARKER = {
    "LEGENDRE_a": "o",
    "LEGENDRE_b": "s",
    "PINES": "^",
    "METRIS": "D",
    "CLENSHAW": "v",
}
C_PAPER = "0.55"  # values transcribed from the paper (neutral grey)
C_REF = "k"  # reference / truth
C_SPAN = "0.93"  # background spans
C_FCR = "#0072B2"  # forward column recursion
C_MSR = "#D55E00"  # modified sub-sectorial recursion (Pines' original)
C_HFCR = "#009E73"  # Helmholtz polynomials via FCR
QUANT_COLOR = ("#0072B2", "#D55E00", "#009E73")  # V / grad V / grad grad V


def tex_escape(s):
    """Escape the characters usetex treats specially in plain text (_ and %)."""
    return s.replace("_", r"\_").replace("%", r"\%") if USE_TEX else s


def method_label(name):
    """Display name of a method, e.g. LEGENDRE_a (underscore escaped for usetex)."""
    return tex_escape(name)


def _tint(color, a):
    """Mix `color` with white: a = 0 -> white, a = 1 -> color."""
    c = np.asarray(mcolors.to_rgb(color))
    return tuple((1.0 - a) + a * c)


# ============================================================
# FIGURE OUTPUT
# ============================================================
_FIG_DIR = Path("results")


def set_figure_dir(path):
    """Set (and create) the folder where save_figure() writes the PDFs."""
    global _FIG_DIR
    _FIG_DIR = Path(path)
    _FIG_DIR.mkdir(parents=True, exist_ok=True)
    return _FIG_DIR


def figure_dir():
    """Current figure folder."""
    return _FIG_DIR


def _write_bytes_retry(fname, data, tries=6, wait=2.0):
    """Write `data` to `fname`, retrying with backoff (cloud-synced folders can time out)."""
    for k in range(tries):
        try:
            fname.parent.mkdir(parents=True, exist_ok=True)
            tmp = fname.with_name(f".{fname.name}.part")
            tmp.write_bytes(data)
            os.replace(tmp, fname)
            return True
        except OSError as err:
            print(f"[Figure] write of {fname.name} failed ({err}); retry {k + 1}/{tries}")
            time.sleep(wait * 2**k)
    return False


def write_output(fname, data, tag="Figure"):
    """
    Write bytes (or str) to a file, robust to a slow cloud-synced (OneDrive) folder.

    Atomic (fname.part, then rename), up to 6 tries with exponential backoff; on failure the
    bytes go to FALLBACK_DIR/<folder>/<file>. Nothing is raised, so a long run is not lost.
    """
    fname = Path(fname)
    if isinstance(data, str):
        data = data.encode("utf-8")
    if _write_bytes_retry(fname, data):
        return fname
    local = FALLBACK_DIR / fname.parent.name / fname.name
    if _write_bytes_retry(local, data, tries=2):
        print(f"[{tag}] WARNING: could not write to {fname.parent}; saved {local} instead")
        return local
    print(f"[{tag}] ERROR: could not write {fname.name}; NOT saved")
    return None


def save_figure(fig, name):
    """Save `fig` as <figure dir>/<name>.pdf, show it if SHOW_FIGS, then close it."""
    buf = io.BytesIO()
    fig.savefig(buf, format="pdf")
    fname = write_output(_FIG_DIR / f"{name}.pdf", buf.getvalue())
    if fname is not None:
        print(f"[Figure] saved {fname}")
    interactive = plt.get_backend().lower() not in {
        "agg", "pdf", "svg", "ps", "cairo", "pgf", "template",
    }
    if SHOW_FIGS and interactive:
        plt.show()
    plt.close(fig)
    return fname


def panel_title(ax, letter, text):
    """MCMCPlots convention: '(a) text' on the left."""
    ax.set_title(f"({letter}) {text}", loc="left")
