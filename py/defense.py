#!/usr/bin/env python3

# Plots for defense slides: https://curiouscoding.nl/slides/defense

import header as h
from header import plot
from sympy import Symbol, Lambda, Max, Min

sigma = 4
w = 24
k = Symbol("k")
ks = None
t = 4

n = 10000000
h.gen(n, sigma)


hl = "#0000cc"
gr = "#00cc00"
ks = range(1, 83)

# Plot 1: Baseline
fs = [
    ("DoubleDecycling", {"lc": gr, "label": "Double Decycling"}),
    ("Random", {"lc": hl, "label": "Random"}),
]
plot(
    "defense/defense-1-before",
    sigma,
    w,
    fs,
    ymin=0.040,
    add=4,
    ncols=4,
    trivial=True,
    tight=False,
    marcais=True,
    ks=ks,
)

# Plot 2: Mod-mini
fs[1][1]["thin"] = True
fs.append(("Random", {"mod": 1, "r": t, "lc": hl, "label": "Mod-mini"}))
plot(
    "defense/defense-2-mod",
    sigma,
    w,
    fs,
    ymin=0.040,
    add=4,
    ncols=4,
    trivial=True,
    tight=False,
    marcais=True,
    ks=ks,
)
# Plot 3: Lower bound
plot(
    "defense/defense-3-lb",
    sigma,
    w,
    fs,
    ymin=0.040,
    add=3,
    ncols=4,
    trivial=True,
    tight=True,
    marcais=True,
    ks=ks,
)

# Plot 4: Small-k schemes
fs.extend(
    [
        ("SusAntiLex", {"lc": "brown", "label": "Sus-anchor"}),
        ("ABB2", {"lc": "teal", "label": "ABB+"}),
    ]
)


plot(
    "defense/defense-4-small-k",
    sigma,
    w,
    fs,
    add=1,
    ncols=4,
    trivial=True,
    ymin=0.04,
    marcais=True,
    ks=ks,
)

# Plot 5: Extended-mod-mini
fs = [
    ("DoubleDecycling", {"lc": gr, "label": "Double Decycling", "thin": True}),
    ("DoubleDecycling", {"lc": gr, "label": "Double Decycling", "mod": 1}),
    ("Random", {"lc": hl, "label": "Random", "thin": True}),
    ("Random", {"mod": 1, "r": t, "lc": hl, "label": "Mod-mini"}),
    ("SusAntiLex", {"lc": "brown", "label": "Sus-anchor", "thin": True}),
    ("SusAntiLex", {"lc": "brown", "label": "Sus-anchor", "mod": 1}),
    ("ABB2", {"lc": "teal", "label": "ABB+", "thin": True}),
    ("ABB2", {"lc": "teal", "label": "ABB+", "mod": 1}),
]

plot(
    "defense/defense-5-ext-mod",
    sigma,
    w,
    fs,
    add=1,
    ncols=4,
    trivial=True,
    ymin=0.04,
    marcais=True,
    ks=ks,
)

# Large-sigma

# Plot 5: Extended-mod-mini
fs = [
    ("DoubleDecycling", {"lc": gr, "label": "Double Decycling", "thin": True}),
    ("DoubleDecycling", {"lc": gr, "label": "Double Decycling", "mod": 1}),
    ("Random", {"lc": hl, "label": "Random", "thin": True}),
    ("Random", {"mod": 1, "r": t, "lc": hl, "label": "Mod-mini"}),
    ("SusAntiLex", {"lc": "brown", "label": "Sus-anchor", "thin": True}),
    ("SusAntiLex", {"lc": "brown", "label": "Sus-anchor", "mod": 1}),
    ("ABB2", {"lc": "teal", "label": "ABB+", "thin": True}),
    ("ABB2", {"lc": "teal", "label": "ABB+", "mod": 1}),
]

plot(
    "defense/defense-5-ext-mod",
    sigma,
    w,
    fs,
    add=1,
    ncols=4,
    trivial=True,
    ymin=0.04,
    marcais=True,
    ks=ks,
)


sigma = 256
t = 1
h.gen(n, sigma)

# Plot 6: large-sigma
fs = [
    ("DoubleDecycling", {"lc": gr, "label": "Double Decycling", "thin": True}),
    ("DoubleDecycling", {"lc": gr, "label": "Double Decycling", "mod": 1}),
    ("Random", {"lc": hl, "label": "Random", "thin": True}),
    ("Random", {"mod": 1, "r": t, "lc": hl, "label": "Mod-mini"}),
    ("SusAntiLex", {"lc": "brown", "label": "Sus-anchor", "thin": True}),
    ("SusAntiLex", {"lc": "brown", "label": "Sus-anchor", "mod": 1}),
    ("ABB2", {"lc": "teal", "label": "ABB+", "thin": True}),
    ("ABB2", {"lc": "teal", "label": "ABB+", "mod": 1}),
]

plot(
    "defense/defense-6-large-sigma",
    sigma,
    w,
    fs,
    add=1,
    ncols=4,
    trivial=True,
    ymin=0.04,
    marcais=True,
    ks=ks,
)
