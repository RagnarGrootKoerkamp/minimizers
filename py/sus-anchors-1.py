#!/usr/bin/env python3
# Plot in sus-anchors paper.
# based on thesis.py.

import header as h
from header import plot
import matplotlib.pyplot as plt

n = 10000000
w = 0
k = 1

ws = list(range(2, 64, 2)) + list(2**i for i in range(6, 11))

fig, axes = plt.subplots(1, 3, figsize=(12, 4.8))

sigmas = [2, 4, 32]
for ax, sigma in zip(axes, sigmas):
    h.gen(n, sigma)
    print(f"Sigma={sigma}, k={k}")

    fs = [
        ("BdAnchor", {"r": 0, "lc": "#00aa00", "label": "Bd-anchor", "k": 0}),
        ("BdAnchor", {"r": 1, "lc": "#00dd00", "label": "Bd-anchor", "k": 1}),
        ("BdAnchor", {"r": 2, "lc": "#00bb00", "label": "Bd-anchor", "k": 2}),
        ("BdAnchor", {"r": 3, "lc": "#009900", "label": "Bd-anchor", "k": 3}),
        ("BdAnchor", {"r": 4, "lc": "#006600", "label": "Bd-anchor", "k": 4}),
        ("BdAnchor", {"r": 5, "lc": "#006600", "label": "Bd-anchor", "k": 5}),
        ("BdAnchor", {"r": 6, "lc": "#003300", "label": "Bd-anchor", "k": 6}),
        ("BdAnchor", {"r": 7, "lc": "#003300", "label": "Bd-anchor", "k": 7}),
        ("BdAnchor", {"r": 8, "lc": "#003300", "label": "Bd-anchor", "k": 8}),
        ("BdAnchor", {"r": 9, "lc": "#003300", "label": "Bd-anchor", "k": 9}),
        ("BdAnchor", {"r": 10, "lc": "#003300", "label": "Bd-anchor", "k": 10}),
        ("DoubleDecycling", {"k": 2, "lc": "black", "label": "Double Decycling"}),
        ("DoubleDecycling", {"k": 3, "lc": "black", "label": "Double Decycling"}),
        ("DoubleDecycling", {"k": 4, "lc": "black", "label": "Double Decycling"}),
        ("DoubleDecycling", {"k": 5, "lc": "black", "label": "Double Decycling"}),
        ("DoubleDecycling", {"k": 6, "lc": "black", "label": "Double Decycling"}),
        ("DoubleDecycling", {"k": 7, "lc": "black", "label": "Double Decycling"}),
        ("DoubleDecycling", {"k": 8, "lc": "black", "label": "Double Decycling"}),
        ("DoubleDecycling", {"k": 9, "lc": "black", "label": "Double Decycling"}),
        ("DoubleDecycling", {"k": 10, "lc": "black", "label": "Double Decycling"}),
        ("Greedy", {"k": 2, "label": "Greedy Miminizer", "lc": "red"}),
        ("Greedy", {"k": 3, "label": "Greedy Miminizer", "lc": "grey"}),
        ("Greedy", {"k": 4, "label": "Greedy Miminizer", "lc": "grey"}),
        ("Greedy", {"k": 5, "label": "Greedy Miminizer", "lc": "grey"}),
        ("Greedy", {"k": 6, "label": "Greedy Miminizer", "lc": "grey"}),
        # Minimizers
        ("Lex", {"k": 2, "label": "Lex Minimizer", "lc": "#660066"}),
        ("Lex", {"k": 3, "label": "Lex Minimizer", "lc": "#0000dd"}),
        ("Lex", {"k": 4, "label": "Lex Minimizer", "lc": "#0000bb"}),
        ("Lex", {"k": 5, "label": "Lex Minimizer", "lc": "#000099"}),
        ("Lex", {"k": 6, "label": "Lex Minimizer", "lc": "#000066"}),
        ("Lex", {"k": 7, "label": "Lex Minimizer", "lc": "#000066"}),
        ("Lex", {"k": 8, "label": "Lex Minimizer", "lc": "#000066"}),
        ("Lex", {"k": 9, "label": "Lex Minimizer", "lc": "#000066"}),
        ("Lex", {"k": 10, "label": "Lex Minimizer", "lc": "#000066"}),
        ("Random", {"k": 2, "lc": "orange", "label": "Random Minimizer"}),
        ("Random", {"k": 3, "lc": "orange", "label": "Random Minimizer"}),
        ("Random", {"k": 4, "lc": "orange", "label": "Random Minimizer"}),
        ("Random", {"k": 5, "lc": "orange", "label": "Random Minimizer"}),
        ("Random", {"k": 6, "lc": "orange", "label": "Random Minimizer"}),
        ("Random", {"k": 7, "lc": "orange", "label": "Random Minimizer"}),
        ("Random", {"k": 8, "lc": "orange", "label": "Random Minimizer"}),
        ("Random", {"k": 9, "lc": "orange", "label": "Random Minimizer"}),
        ("Random", {"k": 10, "lc": "orange", "label": "Random Minimizer"}),
        ("AntiLex", {"k": 2, "label": "Antilex Minimizer", "lc": "#000099"}),
        ("AntiLex", {"k": 3, "label": "Antilex Minimizer", "lc": "#770077"}),
        ("AntiLex", {"k": 4, "label": "Antilex Minimizer", "lc": "#660066"}),
        ("AntiLex", {"k": 5, "label": "Antilex Minimizer", "lc": "#550055"}),
        ("AntiLex", {"k": 6, "label": "Antilex Minimizer", "lc": "#440044"}),
        ("AntiLex", {"k": 7, "label": "Antilex Minimizer", "lc": "#440044"}),
        ("AntiLex", {"k": 8, "label": "Antilex Minimizer", "lc": "#440044"}),
        ("AntiLex", {"k": 9, "label": "Antilex Minimizer", "lc": "#440044"}),
        ("AntiLex", {"k": 10, "label": "Antilex Minimizer", "lc": "#440044"}),
        ("SusLex", {"lc": "#a00080", "label": "Lex SUS-anchor"}),
        # FIXME
        # ("SusABB", {"lc": "#6010b0", "label": "Sus-anchor ABB"}),
        # ("SusABB2", {"lc": "#6010b0", "label": "Sus-anchor ABB+"}),
        # ("SusAlternating", {"lc": "#0000ff", "label": "Sus-anchor Alternating"}),
        ("SusAntiLex", {"lc": "#0000ff", "label": "Antilex SUS-anchor"}),
        # ("Alternating", {"k": 2, "label": "Alternating", "lc": "#ff0000"}),
        # ("Alternating", {"k": 3, "label": "Alternating", "lc": "#dd0000"}),
        # ("Alternating", {"k": 4, "label": "Alternating", "lc": "#bb0000"}),
        # ("Alternating", {"k": 5, "label": "Alternating", "lc": "#990000"}),
        # ("Alternating", {"k": 6, "label": "Alternating", "lc": "#660000"}),
        # ("Alternating", {"k": 7, "label": "Alternating", "lc": "#660000"}),
        # ("Alternating", {"k": 8, "label": "Alternating", "lc": "#660000"}),
        # ("Alternating", {"k": 9, "label": "Alternating", "lc": "#660000"}),
        # ("Alternating", {"k": 10, "label": "Alternating", "lc": "#660000"}),
        # ("ABB", {"k": 2, "label": "ABB", "lc": "#880088"}),
        # ("ABB", {"k": 3, "label": "ABB", "lc": "#770077"}),
        # ("ABB", {"k": 4, "label": "ABB", "lc": "#660066"}),
        # ("ABB", {"k": 5, "label": "ABB", "lc": "#550055"}),
        # ("ABB", {"k": 6, "label": "ABB", "lc": "#440044"}),
        # ("ABB", {"k": 7, "label": "ABB", "lc": "#440044"}),
        # ("ABB", {"k": 8, "label": "ABB", "lc": "#440044"}),
        # ("ABB", {"k": 9, "label": "ABB", "lc": "#440044"}),
        # ("ABB", {"k": 10, "label": "ABB+", "lc": "#440044"}),
        # ("ABB2", {"k": 2, "label": "ABB+", "lc": "#880088"}),
        # ("ABB2", {"k": 3, "label": "ABB+", "lc": "#770077"}),
        # ("ABB2", {"k": 4, "label": "ABB+", "lc": "#660066"}),
        # ("ABB2", {"k": 5, "label": "ABB+", "lc": "#550055"}),
        # ("ABB2", {"k": 6, "label": "ABB+", "lc": "#440044"}),
        # ("ABB2", {"k": 7, "label": "ABB+", "lc": "#440044"}),
        # ("ABB2", {"k": 8, "label": "ABB+", "lc": "#440044"}),
        # ("ABB2", {"k": 9, "label": "ABB+", "lc": "#440044"}),
        # ("ABB2", {"k": 10, "label": "ABB+", "lc": "#440044"}),
        # ("ABB", {"label": "ABB"}),
        # ("ABB2", {"lc": "teal", "label": "ABB+"}),
        # # MORE
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
        best_only=True,
    )
    ax.set_title(f"$\\sigma={sigma}$", fontsize=11)
    if sigma > 2:
        ax.set_ylabel(None)


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
