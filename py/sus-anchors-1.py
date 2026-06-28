#!/usr/bin/env python3
# Plot in sus-anchors paper.
# based on thesis.py.

import header as h
from header import plot
import matplotlib.pyplot as plt

n = 1000000
w = 0
k = 1

ws = list(range(2, 64, 2)) + list(2**i for i in range(6, 11))

fig, axes = plt.subplots(1, 3, figsize=(12, 4.8))

sigmas = [2, 4, 32]
for ax, sigma in zip(axes, sigmas):
    h.gen(n, sigma)
    print(f"Sigma={sigma}, k={k}")

    fs = [
        # ("BdAnchor", {"r": 0, "lc": "#00ff00", "label": "Bd-anchor r=0"}),
        # ("BdAnchor", {"r": 1, "lc": "#00dd00", "label": "Bd-anchor r=1"}),
        # ("BdAnchor", {"r": 2, "lc": "#00bb00", "label": "Bd-anchor r=2"}),
        # ("BdAnchor", {"r": 3, "lc": "#009900", "label": "Bd-anchor r=3"}),
        # ("BdAnchor", {"r": 4, "lc": "#006600", "label": "Bd-anchor r=4"}),
        # ("BdAnchor", {"r": 6, "lc": "#003300", "label": "Bd-anchor r=6"}),
        ("SusLex", {"lc": "#a00080", "label": "Sus-anchor Lex"}),
        # FIXME
        # ("SusABB", {"lc": "#6010b0", "label": "Sus-anchor ABB"}),
        # ("SusABB2", {"lc": "#6010b0", "label": "Sus-anchor ABB+"}),
        ("SusAntiLex", {"lc": "#0000ff", "label": "Sus-anchor AntiLex"}),
        # Minimizers
        ("Lex", {"k": 2, "label": "Lex", "lc": "#0000ff"}),
        ("Lex", {"k": 3, "label": "Lex", "lc": "#0000dd"}),
        ("Lex", {"k": 4, "label": "Lex", "lc": "#0000bb"}),
        ("Lex", {"k": 5, "label": "Lex", "lc": "#000099"}),
        ("Lex", {"k": 6, "label": "Lex", "lc": "#000066"}),
        ("Random", {"k": 2, "lc": "orange", "label": "Random"}),
        ("Random", {"k": 3, "lc": "orange", "label": "Random"}),
        ("Random", {"k": 4, "lc": "orange", "label": "Random"}),
        ("Random", {"k": 5, "lc": "orange", "label": "Random"}),
        ("Random", {"k": 6, "lc": "orange", "label": "Random"}),
        ("Alternating", {"k": 2, "label": "Alternating", "lc": "#ff0000"}),
        ("Alternating", {"k": 3, "label": "Alternating", "lc": "#dd0000"}),
        ("Alternating", {"k": 4, "label": "Alternating", "lc": "#bb0000"}),
        ("Alternating", {"k": 5, "label": "Alternating", "lc": "#990000"}),
        ("Alternating", {"k": 6, "label": "Alternating", "lc": "#660000"}),
        ("AntiLex", {"k": 2, "label": "AntiLex", "lc": "#880088"}),
        ("AntiLex", {"k": 3, "label": "AntiLex", "lc": "#770077"}),
        ("AntiLex", {"k": 4, "label": "AntiLex", "lc": "#660066"}),
        ("AntiLex", {"k": 5, "label": "AntiLex", "lc": "#550055"}),
        ("AntiLex", {"k": 6, "label": "AntiLex", "lc": "#440044"}),
        # ("ABB", {"label": "ABB"}),
        # ("ABB2", {"lc": "teal", "label": "ABB+"}),
        ("DoubleDecycling", {"k": 2, "lc": "#dddddd", "label": "Double Decycling"}),
        ("DoubleDecycling", {"k": 3, "lc": "#bbbbbb", "label": "Double Decycling"}),
        ("DoubleDecycling", {"k": 4, "lc": "#999999", "label": "Double Decycling"}),
        ("DoubleDecycling", {"k": 5, "lc": "#777777", "label": "Double Decycling"}),
        ("DoubleDecycling", {"k": 6, "lc": "#555555", "label": "Double Decycling"}),
        # (
        #     "OpenClosed",
        #     {"k": 2, "r": t, "closed": 1, "lc": "#00dd00", "label": "C=Miniception"},
        # ),
        # (
        #     "OpenClosed",
        #     {"k": 3, "r": t, "closed": 1, "lc": "#00dd00", "label": "C=Miniception"},
        # ),
        # (
        #     "OpenClosed",
        #     {"k": 4, "r": t, "closed": 1, "lc": "#00dd00", "label": "C=Miniception"},
        # ),
        # (
        #     "OpenClosed",
        #     {"k": 5, "r": t, "closed": 1, "lc": "#00dd00", "label": "C=Miniception"},
        # ),
        # (
        #     "OpenClosed",
        #     {"k": 6, "r": t, "closed": 1, "lc": "#00dd00", "label": "C=Miniception"},
        # ),
        # ("Greedy", {"k": 2, "label": "Greedy", "lc": "grey"}),
        # ("Greedy", {"k": 3, "label": "Greedy", "lc": "grey"}),
        # ("Greedy", {"k": 4, "label": "Greedy", "lc": "grey"}),
        # ("Greedy", {"k": 5, "label": "Greedy", "lc": "grey"}),
        # ("Greedy", {"k": 6, "label": "Greedy", "lc": "grey"}),
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

plt.savefig("sus-anchors-1.png", bbox_inches="tight", dpi=400)
plt.savefig("sus-anchors-1.svg", bbox_inches="tight")
plt.close()
