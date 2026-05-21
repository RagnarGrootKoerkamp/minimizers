#!/usr/bin/env python3
# Experiments for selection schemes.
# Initially copies from thesis.py.

import header as h
from header import plot
from sympy import Symbol, Lambda, Max, Min

n = 1000009
w = 0

# for sigma in [2, 4, 32, 256]:
for sigma in [4]:
    h.gen(n, sigma)
    # for k in [1, 2, 3, 4]:
    for k in [2, 3, 4]:
        print(f"Sigma={sigma}, k={k}")
        # ws = range(2, 64)
        ws = list(range(2, 64)) + list(2**i for i in range(6, 11))
        fs = [
            ("BdAnchor", {"r": 0, "lc": "#00ff00", "label": "Bd-anchor r=0"}),
            ("BdAnchor", {"r": 1, "lc": "#00dd00", "label": "Bd-anchor r=1"}),
            ("BdAnchor", {"r": 2, "lc": "#00bb00", "label": "Bd-anchor r=2"}),
            ("BdAnchor", {"r": 3, "lc": "#009900", "label": "Bd-anchor r=3"}),
            ("BdAnchor", {"r": 4, "lc": "#006600", "label": "Bd-anchor r=4"}),
            ("BdAnchor", {"r": 6, "lc": "#003300", "label": "Bd-anchor r=6"}),
            ("SusLex", {"lc": "#a00080", "label": "Sus-anchor Lex"}),
            ("SusABB", {"lc": "#6010b0", "label": "Sus-anchor ABB", "sampling": 0}),
            # ("SusABB2", {"lc": "#2005e0", "label": "Sus-anchor ABB+"}),
            ("SusAntiLex", {"lc": "#0000ff", "label": "Sus-anchor AntiLex"}),
        ]
        plot(
            f"tmp-s{sigma}-k{k}",
            # f"tmp",
            sigma,
            w,
            fs,
            k=k,
            plot_w=True,
            loose=False,
            tight=True,
            ncols=3,
            add=3,
            ws=ws,
            ymin=1.951,
            ymax=2.90,
            df=True,
            split=64,
            lw=0.3,
        )
