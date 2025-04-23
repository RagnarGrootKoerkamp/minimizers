#!/usr/bin/env python3

# Plots for greedy mini analysis: https://curiouscoding.nl/posts

import header as h
from sympy import Symbol, Lambda, Max, Min

w = 24
ks = range(3, 16)

# n = 10000000
n = 100000

hl = "#0000cc"
gr = "#00cc00"


def fst(t):
    return [
        # only for w=5
        ("Random", {"lc": hl, "label": "Random", "thin": True}),
        ("Random", {"mod": 1, "r": t, "lc": hl, "label": "Mod-mini"}),
        # always
        ("DoubleDecycling", {"lc": gr, "label": "Double Decycling", "mod": 1}),
        ("Greedy", {"label": "Greedy", "lc": "grey"}),
        # only for sigma=4
        ("SusAntiLex", {"lc": "brown", "label": "Sus-anchor", "mod": 1}),
        ("ABB2", {"lc": "teal", "label": "ABB+", "mod": 1}),
        # only for sigma=256
        ("ThresholdABB2", {"label": "Theshold+AntiLex", "thr": 100, "mod": 1}),
    ]


def plot(fs, *, add):
    global w, sigma, ks
    h.plot(
        f"greedymini/w{w}-s{sigma}",
        sigma,
        w,
        fs,
        add=add,
        ncols=3,
        trivial=True,
        ymin=0.99 / w,
        ymax=2.02 / (w + 1),
        ks=ks,
        title=f"w={w}, σ={sigma}",
        ilp=True,
    )


sigma = 2
h.gen(n, sigma)
fs = fst(1)

for w in [3, 4]:
    plot(fs[:6], add=2)

sigma = 4
h.gen(n, sigma)
w = 5
fs = fst(2)

plot(fs[:6], add=2)

sigma = 256
h.gen(n, sigma)
fs = fst(1)
fs.pop(4)
fs.pop(4)

for w in [3, 4, 5, 6, 8, 10, 12]:
    plot(fs[2 * (w != 5) :], add=3 + (w != 5))
