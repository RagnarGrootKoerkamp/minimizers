#!/usr/bin/env python3
# Plot in sus-anchors paper.
# based on thesis.py.

import header as h
from header import plot
import matplotlib.pyplot as plt

n = 10000000
w = 0
k = 1

ws = list(range(2, 64)) + list(2**i for i in range(6, 11))

fig, axes = plt.subplots(1, 3, figsize=(12, 4.8))

sigmas = [2, 4, 32]
for ax, sigma in zip(axes, sigmas):
    h.gen(n, sigma)
    print(f"Sigma={sigma}, k={k}")

    fs = [
        ("BdAnchor", {"r": 0, "lc": "#00ff00", "label": "Bd-anchor r=0"}),
        ("BdAnchor", {"r": 1, "lc": "#00dd00", "label": "Bd-anchor r=1"}),
        ("BdAnchor", {"r": 2, "lc": "#00bb00", "label": "Bd-anchor r=2"}),
        ("BdAnchor", {"r": 3, "lc": "#009900", "label": "Bd-anchor r=3"}),
        ("BdAnchor", {"r": 4, "lc": "#006600", "label": "Bd-anchor r=4"}),
        ("BdAnchor", {"r": 6, "lc": "#003300", "label": "Bd-anchor r=6"}),
        ("SusLex", {"lc": "#a00080", "label": "Sus-anchor Lex"}),
        # ("SusABB", {"lc": "#6010b0", "label": "Sus-anchor ABB"}),
        # ("SusABB2", {"lc": "#6010b0", "label": "Sus-anchor ABB+"}),
        ("SusAntiLex", {"lc": "#0000ff", "label": "Sus-anchor AntiLex"}),
    ]
    fs = [x for x in fs if x[1].get("r", 1000) >= k - 1]

    ymax = {2: 4.0, 4: 2.7, 32: 2.2}[sigma]
    plot(
        f"tmp-s{sigma}-k{k}",
        sigma,
        w,
        fs,
        ax=ax,
        save=False,
        show_legend=False,
        k=k,
        plot_w=True,
        loose=False,
        trivial=False,
        tight=True,
        ncols=3,
        add=0,
        ws=ws,
        ymax=ymax,
        ymin=2.0 - (ymax - 2.0) * 0.05,
        df=True,
        split=64,
        lw=0.3,
        height=4.3,
    )
    ax.set_title(f"$\\sigma={sigma}$", fontsize=11)


# Shared legend from the first subplot's handles
handles, labels = axes[0].get_legend_handles_labels()
fig.legend(
    handles,
    labels,
    loc="lower center",
    bbox_to_anchor=(0.5, -0.12),
    ncols=3,
    fontsize=9,
)

plt.savefig(f"sus-anchors-1.png", bbox_inches="tight", dpi=400)
plt.savefig(f"sus-anchors-1.svg", bbox_inches="tight")
plt.close()
