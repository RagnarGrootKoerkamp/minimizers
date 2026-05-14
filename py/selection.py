#!/usr/bin/env python3
# Experiments for selection schemes.
# Initially copies from thesis.py.

import header as h
from header import plot
from sympy import Symbol, Lambda, Max, Min

sigma = 2
w = 24
k = Symbol("k")
ks = None
t = 4

n = 100000
h.gen(n, sigma)

k = 1
ws = range(2, 80)
fs = [
    # ("BdAnchor", {"r": 0, "lc": "#00ff00", "label": "Bd-anchor r=0"}),
    # ("BdAnchor", {"r": 1, "lc": "#00dd00", "label": "Bd-anchor r=1"}),
    # ("BdAnchor", {"r": 2, "lc": "#00bb00", "label": "Bd-anchor r=2"}),
    # ("BdAnchor", {"r": 3, "lc": "#009900", "label": "Bd-anchor r=3"}),
    # ("SusLex", {"lc": "#a00080", "label": "Sus-anchor Lex"}),
    # ("SusABB2", {"lc": "#6010b0", "label": "Sus-anchor ABB+"}),
    ("SusAntiLex", {"lc": "#0000ff", "label": "Sus-anchor AntiLex"}),
    (
        "SusRandomLex",
        {"lc": "red", "label": "SUS", "seed": 576843},
    ),
    (
        "SusRandomLex",
        {"lc": "orange", "label": "SUS", "seed": 2332244},
    ),
    (
        "SusRandomLex",
        {"lc": "pink", "label": "SUS", "seed": 23408432},
    ),
    (
        "SusRandomLex",
        {"lc": "green", "label": "SUS", "seed": 823424},
    ),
    (
        "SusRandomLex",
        {"lc": "magenta", "label": "SUS", "seed": 9234},
    ),
    (
        "SusRandomLex",
        {"lc": "magenta", "label": "SUS", "seed": 923492840},
    ),
    (
        "SusRandomLex",
        {"lc": "magenta", "label": "SUS", "seed": 92342834},
    ),
]
plot(
    "tmp",
    sigma,
    w,
    fs,
    plot_w=True,
    loose=True,
    tight=True,
    ymin=1.9,
    ymax=2.5,
    ncols=3,
    add=1,
    bold_last=True,
    ws=ws,
    df=True,
)

plot(
    "tmp-diff",
    sigma,
    w,
    fs,
    # diff="tight",
    diff="loose",
    plot_w=True,
    # loose=True,
    # tight=False,
    trivial=False,
    ymin=-0.0001,
    ymax=0.1,
    ncols=3,
    add=1,
    # bold_last=True,
    ws=ws,
    rol=False,
)
