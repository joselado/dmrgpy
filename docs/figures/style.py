"""Shared matplotlib style of the user-guide figures.

Every docs/figures/<id>/plot.py uses it the same way:

    import os, sys
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    import style

    def draw(theme):
        c = style.colors(theme)            # c["ink"], c["ink2"], c["series"][k]
        fig, ax = style.figure(1)          # 1, 2 or 3 panels
        ax.plot(x, y, color=c["series"][0], lw=style.LW, label="DMRG")
        ax.plot(x, yref, **style.reference(theme), label="exact")
        ax.legend()
        return fig

    style.main(draw, "<id>")   # parses --theme, writes <id>.png and <id>-dark.png

The PNG files go to docs/guide/img/. The plot is transparent, so the page
color shows through; the ink colors switch with the theme. This module
never imports dmrgpy.
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
IMG = os.path.join(os.path.dirname(HERE), "guide", "img")

INK = {"light": "#0b0b0b", "dark": "#f0efec"}
INK2 = {"light": "#52514e", "dark": "#c3c2b7"}

# blue, orange, aqua, yellow, magenta, green, violet, red
SERIES = {
    "light": ["#2a78d6", "#eb6834", "#1baf7a", "#eda100",
              "#e87ba4", "#008300", "#6250d6", "#e34948"],
    "dark": ["#3987e5", "#d95926", "#199e70", "#c98500",
             "#d55181", "#008300", "#9085e9", "#e66767"],
}

LW = 1.6          # computed data, lines
MS = 5            # computed data, markers
LW_REF = 1.2      # reference curves (exact, ED, free fermions)
CMAP = "magma"    # heat maps, both themes
SIZES = {1: (6.4, 4.0), 2: (9.6, 3.8), 3: (11.0, 3.6)}
DPI = 160


def colors(theme):
    """Ink, secondary ink and the series colors of a theme."""
    return {"ink": INK[theme], "ink2": INK2[theme], "series": SERIES[theme]}


def reference(theme):
    """Line keywords of a reference curve: ink, dashed, 1.2 pt."""
    return {"color": INK[theme], "ls": "--", "lw": LW_REF}


def heatmap_line():
    """Line keywords of a curve drawn on a heat map: white, dashed."""
    return {"color": "white", "ls": "--", "lw": LW_REF}


def apply(theme):
    """Set the rcParams of a theme. Called by figure()."""
    ink, ink2 = INK[theme], INK2[theme]
    plt.rcdefaults()
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 10,
        "axes.labelsize": 11,
        "legend.fontsize": 9,
        "legend.frameon": False,
        "axes.grid": False,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.edgecolor": ink,
        "axes.labelcolor": ink,
        "axes.facecolor": "none",
        "figure.facecolor": "none",
        "savefig.facecolor": "none",
        "text.color": ink,
        "xtick.color": ink,
        "ytick.color": ink,
        "xtick.labelcolor": ink2,
        "ytick.labelcolor": ink2,
        "legend.labelcolor": ink,
        "axes.prop_cycle": plt.cycler(color=SERIES[theme]),
        "image.cmap": CMAP,
        "lines.linewidth": LW,
        "lines.markersize": MS,
    })


def figure(npanels=1, theme=None, **kwargs):
    """A figure with 1, 2 or 3 panels side by side, in the plan's sizes."""
    if theme is not None:
        apply(theme)
    fig, axs = plt.subplots(1, npanels, figsize=SIZES[npanels], **kwargs)
    return fig, axs


def colorbar(fig, mappable, ax, label, theme):
    """A colorbar with its outline, ticks and label in the ink color."""
    cb = fig.colorbar(mappable, ax=ax)
    cb.set_label(label, color=INK[theme])
    cb.outline.set_edgecolor(INK[theme])
    cb.ax.tick_params(colors=INK[theme], labelcolor=INK2[theme])
    return cb


def save(fig, figid, theme):
    """Write docs/guide/img/<id>.png (light) or <id>-dark.png (dark)."""
    os.makedirs(IMG, exist_ok=True)
    name = figid + ("-dark" if theme == "dark" else "") + ".png"
    path = os.path.join(IMG, name)
    fig.savefig(path, dpi=DPI, bbox_inches="tight", transparent=True)
    plt.close(fig)
    print("wrote", path)
    return path


def main(draw, figid):
    """Parse --theme light|dark (default: both), call draw(theme), save."""
    p = argparse.ArgumentParser()
    p.add_argument("--theme", choices=["light", "dark"], default=None)
    args = p.parse_args()
    for theme in ([args.theme] if args.theme else ["light", "dark"]):
        apply(theme)
        save(draw(theme), figid, theme)
